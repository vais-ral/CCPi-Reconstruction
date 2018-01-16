#include "CCPiReconstructionParaviewImpl.h"

#include "CCPiParaviewUserInterface.h"

//maybe need these? check
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkIntArray.h"

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

int CCPiReconstructionParaviewImpl::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    //! Update Paraview UI with Progress bar
    CCPiParaviewUserInterface userInterface(this);
    this->UpdateProgress(0.1);

    vtkImageData *pixels = vtkImageData::GetData(inputVector[0]);
    vtkImageData *angles = vtkImageData::GetData(inputVector[1]);
    vtkImageData *output = vtkImageData::GetData(outputVector->GetInformationObject(0));

    int pixel_dims[3];
    pixels->GetDimensions(pixel_dims);
    std::cout << "Dimensions: " << pixel_dims[0] << ", " << pixel_dims[1] << ", " <<  pixel_dims[2] << std::endl;

    float* pixel_vals;

    /*boost::multi_array<float, 3> pixels_arr(boost::extents[pixel_dims[2]][pixel_dims[1]][pixel_dims[0]]);

    std::cout << "test 3" << std::endl;

    for (int x = 0; x < pixel_dims[2]; x++) {
        for (int y = 0; y < pixel_dims[1]; y++) {
            for (int z = 0; z < pixel_dims[0]; z++) {
                float value = pixels->GetScalarComponentAsFloat(z,y,x,0);
                pixels_arr[x][y][z] = value;
            }
        }
    }*/

    boost::multi_array<float, 3> pixels_arr(boost::extents[pixel_dims[0]][pixel_dims[1]][pixel_dims[2]]);

    for (int z = 0; z < pixel_dims[2]; z++) {
        for (int y = 0; y < pixel_dims[1]; y++) {
            for (int x = 0; x < pixel_dims[0]; x++) {
                float value = pixels->GetScalarComponentAsFloat(x,y,z,0);
                pixels_arr[x][y][z] = value;
            }
        }
    }

    this->UpdateProgress(0.5);
    std::cout << "New Dimensions: " << pixels_arr.shape()[0] << ", " << pixels_arr.shape()[1] << ", " <<  pixels_arr.shape()[2] << std::endl;

    int angles_dim[3];
    angles->GetDimensions(angles_dim);
    std::cout << "Angles length: " << angles_dim[0] << std::endl;

    boost::multi_array<float, 1> angles_arr(boost::extents[angles_dim[0]]);

    for (int a = 0; a < angles_dim[0]; a++) {
        float value = angles->GetScalarComponentAsFloat(a,0,0,0);
        angles_arr[a] = value;
    }

    CCPi::reconstruction_alg *algorithm = 0;
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

    const int blocking_factor = 0;
    bool beam_hardening = false;
    CCPi::instrument *instrument = new CCPi::Diamond();

    //do the thing
    machine::initialise(0);

    voxel_data *voxels = reconstruct(instrument, algorithm, pixels_arr, angles_arr,
                                     RotationCentre, Resolution, blocking_factor,
                                     beam_hardening, false);
    machine::exit();

    delete instrument;
    delete algorithm;

    this->UpdateProgress(0.8);

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
    std::cout << "Output Dimensions: " << dims[0] << ", " << dims[1] << ", " <<  dims[2] << std::endl;

    output->SetDimensions(dims);
    output->AllocateScalars(VTK_FLOAT,1);

    std::cout << (*voxels)[0][0][0] << std::endl;

    std::cout << "New extents: " << output->GetExtent()[0] << " " << output->GetExtent()[1] << " " << output->GetExtent()[2] << " " << output->GetExtent()[3] << " " << output->GetExtent()[4] << " " << output->GetExtent()[5]  << std::endl;
    
    if (voxels == 0) {
        output->SetScalarComponentFromFloat(0,0,0,0,0);
    } else {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {
                    output->SetScalarComponentFromFloat(i,j,k,0,(*voxels)[i][j][k]);
                }
            }
        }
    }

    delete voxels;

    this->SetProgress(0.99);

    return 1;
}

int CCPiReconstructionParaviewImpl::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 1)
    {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
        //info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
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