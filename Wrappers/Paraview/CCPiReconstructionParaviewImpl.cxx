#include "CCPiReconstructionParaviewImpl.h"

#include "CCPiParaviewUserInterface.h"

//maybe need these? check
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkIntArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <iostream>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/multi_array.hpp>
#include "algorithms.hpp"
#include "instruments.hpp"
#include "cgls.hpp"
#include "mlem.hpp"
#include "sirt.hpp"
#include "tv_reg.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "mpi.hpp"
#include "base_types.hpp"
#include "utils.hpp"
#include "math.h"

vtkStandardNewMacro(CCPiReconstructionParaviewImpl);

CCPiReconstructionParaviewImpl::CCPiReconstructionParaviewImpl()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

CCPiReconstructionParaviewImpl::~CCPiReconstructionParaviewImpl()
{


}

void CCPiReconstructionParaviewImpl::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}


//initialise the variables: converts the vtkImageData inputs to boost multiarrays, 
//and create instrument and algorithm - used by both RequestInformation and
//RequestData so extracted here to reduce code repetition.
int CCPiReconstructionParaviewImpl::initialise_variables(vtkSmartPointer<vtkImageData> pixels, vtkSmartPointer<vtkImageData> angles)
{
    int pixel_dims[3];
    pixels->GetDimensions(pixel_dims);

    int angles_dim[3];
    angles->GetDimensions(angles_dim);

    //check dimensionality of angle array and either x or z axis
    if(pixel_dims[0] == angles_dim[0]) {
        //angle array matches x - don't have to reverse
        pixels_arr.resize(boost::extents[pixel_dims[0]][pixel_dims[1]][pixel_dims[2]]);

        for (int z = 0; z < pixel_dims[2]; z++) {
            for (int y = 0; y < pixel_dims[1]; y++) {
                for (int x = 0; x < pixel_dims[0]; x++) {
                    float value = pixels->GetScalarComponentAsFloat(x,y,z,0);
                    pixels_arr[x][y][z] = value;
                }
            }
        }

    } else if (pixel_dims[2] == angles_dim[0]) {
        //angle array mayches z - have to reverse
        pixels_arr.resize(boost::extents[pixel_dims[2]][pixel_dims[1]][pixel_dims[0]]);
        
        for (int x = 0; x < pixel_dims[2]; x++) {
            for (int y = 0; y < pixel_dims[1]; y++) {
                for (int z = 0; z < pixel_dims[0]; z++) {
                    float value = pixels->GetScalarComponentAsFloat(z,y,x,0);
                    pixels_arr[x][y][z] = value;
                }
            }
        }

    } else {
        //angle array does not match either the x or z dimension - error
        vtkErrorMacro("The length angle array does not match the input");
        return 0;
    }

    angles_arr.resize(boost::extents[angles_dim[0]]);

    for (int a = 0; a < angles_dim[0]; a++) {
        float value = angles->GetScalarComponentAsFloat(a,0,0,0);
        angles_arr[a] = value;
    }

    blocking_factor = 0;
    instrument = new CCPi::Diamond();

    algorithm = 0;
    //TODO implement all?
    switch(Algorithm) {
        case 0:
            algorithm = new CCPi::cgls_3d(Iterations);
            break;
        case 1:
            algorithm = new CCPi::sirt(Iterations);
            break;
        case 2:
            algorithm = new CCPi::mlem(Iterations);
            break;
    }

    return 1;
}


int CCPiReconstructionParaviewImpl::RequestInformation(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    std::cout << "RequestInformation" << std::endl;

    vtkImageData *pixels = vtkImageData::GetData(inputVector[0]);
    vtkImageData *angles = vtkImageData::GetData(inputVector[1]);

    int error = initialise_variables(pixels, angles);

    //need to propogate error in case initialise_variables fails
    if(error == 0) {
        return 0;
    }

    //need to calculate dimensions in advance to set the extent properly, so call 
    //calculate_dimensions to calculate them for us
    machine::initialise(0);
    int *dimensions = calculate_dimensions(instrument, algorithm, pixels_arr, angles_arr,
                                           RotationCentre, Resolution, blocking_factor);
    machine::exit();

    int outWholeExt[6];
    outWholeExt[0] = 0;
    outWholeExt[1] = dimensions[2];
    outWholeExt[2] = 0;
    outWholeExt[3] = dimensions[1];
    outWholeExt[4] = 0;
    outWholeExt[5] = dimensions[0];

    double outOrigin[3];
    double outSpacing[3];

    outOrigin[0] = 0.0;
    outOrigin[1] = 0.0;
    outOrigin[2] = 0.0;

    outSpacing[0] = 1.0;
    outSpacing[1] = 1.0;
    outSpacing[2] = 1.0;

    //set some stuff for the extent
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outWholeExt, 6);

    outInfo->Set(vtkDataObject::ORIGIN(), outOrigin, 3);
    outInfo->Set(vtkDataObject::SPACING(), outSpacing, 3);

    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);

    //reset execution count
    execution_count = 0;

    return 1;
}

int CCPiReconstructionParaviewImpl::RequestUpdateExtent(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    int inExt[6];


    vtkInformation *pixelsInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *anglesInfo = inputVector[1]->GetInformationObject(0);

    pixelsInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
    pixelsInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

    anglesInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inExt);
    anglesInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), inExt, 6);

    return 1;
}

int CCPiReconstructionParaviewImpl::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    //count executions - currently it repeats for no reason whilst being applied
    //and also when views are changed. TODO: find out what is causing the repeat.
    //UPDATE: removing the outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), outWholeExt, 6); line
    //in RequestInformation fixes the repeat problem, but then the extents are wrong... 
    std::cout << "Execution Count: "  << execution_count << std::endl;

    //! Update Paraview UI with Progress bar
    CCPiParaviewUserInterface userInterface(this);
    this->UpdateProgress(0.1);

    vtkImageData *pixels = vtkImageData::GetData(inputVector[0]);
    vtkImageData *angles = vtkImageData::GetData(inputVector[1]);

    vtkImageData *output = vtkImageData::GetData(outputVector->GetInformationObject(0));

    //if execution count is not 0, that means we have repeated RequestData without
    //calling RequestInformation, so the variables have been deleted so we need
    //to reinitialise. Check if this still needed after repeat execution bug fixed.
    if(execution_count != 0) {
        initialise_variables(pixels, angles);
    }

    this->UpdateProgress(0.5);

    bool beam_hardening = false;

    //apply the algorithm
    machine::initialise(0);

    voxel_data *voxels = reconstruct(instrument, algorithm, pixels_arr, angles_arr,
                                     RotationCentre, Resolution, blocking_factor,
                                     beam_hardening, false);
    machine::exit();

    //delete to free up memory
    delete instrument;
    delete algorithm;

    //resize boost arrays to 0,0,0 to deallocate memory
    pixels_arr.resize(boost::extents[0][0][0]);
    angles_arr.resize(boost::extents[0]);

    this->UpdateProgress(0.8);

    //need to set extent (again!) so get dimensions from voxels
    int dims[3];
    if (voxels == 0) {
        dims[0] = 1;
        dims[1] = 1;
        dims[2] = 1;
    } else {
        // Todo - remove buffered region
        dims[0] = voxels->shape()[0];
        dims[1] = voxels->shape()[1];
        dims[2] = voxels->shape()[2];
    }
    int outWholeExt[6];
    outWholeExt[0] = 0;
    outWholeExt[1] = dims[2]-1;
    outWholeExt[2] = 0;
    outWholeExt[3] = dims[1]-1;
    outWholeExt[4] = 0;
    outWholeExt[5] = dims[0]-1;

    output->SetOrigin(0,0,0);
    output->SetSpacing(1,1,1);
    output->SetExtent(outWholeExt);
    output->AllocateScalars(VTK_FLOAT,1);

    if (voxels == 0) {
        output->SetScalarComponentFromFloat(0,0,0,0,0);
    } else {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {
                    //reverse x and z as converting from numpy back to vtk
                    output->SetScalarComponentFromFloat(k,j,i,0,(*voxels)[i][j][k]);
                }
            }
        }
    }

    delete voxels;

    this->SetProgress(0.99);

    execution_count = execution_count + 1;

    return 1;
}

int CCPiReconstructionParaviewImpl::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  else
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
return 1;
}

int CCPiReconstructionParaviewImpl::FillOutputPortInformation(int port, vtkInformation* info)
{
    if (port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
    return 1;
}