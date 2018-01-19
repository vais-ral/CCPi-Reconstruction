#ifndef CCPIRECONSTRUCTIONPARAVIEWIMPL_H
#define CCPIRECONSTRUCTIONPARAVIEWIMPL_H

#include "vtkImageAlgorithm.h"
#include "base_types.hpp"
#include "algorithms.hpp"
#include "instruments.hpp"
#include "vtkSmartPointer.h"

#include <string>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bp = boost::python;
namespace np = boost::python::numpy;

class VTK_EXPORT CCPiReconstructionParaviewImpl : public vtkImageAlgorithm 
{
public:
  static CCPiReconstructionParaviewImpl* New();
  vtkTypeMacro(CCPiReconstructionParaviewImpl, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  //TODO: algorithm?
  vtkGetMacro(Algorithm, int);
  vtkSetMacro(Algorithm, int);

  vtkGetMacro(Iterations, int);
  vtkSetMacro(Iterations, int);

  vtkGetMacro(Resolution, int);
  vtkSetMacro(Resolution, int);

  //vtkSetVector3Macro(RotationCentre, float);
  //vtkGetVector3Macro(RotationCentre, float);

  vtkGetMacro(RotationCentre, double);
  vtkSetMacro(RotationCentre, double);

protected:
  CCPiReconstructionParaviewImpl();
  ~CCPiReconstructionParaviewImpl();

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);
  int RequestInformation(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *);

  int Algorithm;
  int Iterations;
  int Resolution;
  double RotationCentre;

  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);

private:
  CCPiReconstructionParaviewImpl(const CCPiReconstructionParaviewImpl&);  // Not implemented.
  void operator=(const CCPiReconstructionParaviewImpl&);  // Not implemented.
  int intialise_variables(vtkSmartPointer<vtkImageData> pixels, vtkSmartPointer<vtkImageData> angles);

  //to save recalculation, initialise in RequestInformation, reuse from here in
  //Request Data, then delete in RequestData once finished.
  CCPi::reconstruction_alg *algorithm;
  int blocking_factor;
  CCPi::instrument *instrument;
  boost::multi_array<float, 3> pixels_arr;
  boost::multi_array<float, 1> angles_arr;
  bool deleted_vars = false;

  int execution_count = 0;
};

#endif