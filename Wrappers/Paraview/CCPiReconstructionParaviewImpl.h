#ifndef CCPIRECONSTRUCTIONPARAVIEWIMPL_H
#define CCPIRECONSTRUCTIONPARAVIEWIMPL_H

#include "vtkImageAlgorithm.h"
#include "base_types.hpp"
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

  int Algorithm;
  int Iterations;
  int Resolution;
  double RotationCentre;

  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);

private:
  CCPiReconstructionParaviewImpl(const CCPiReconstructionParaviewImpl&);  // Not implemented.
  void operator=(const CCPiReconstructionParaviewImpl&);  // Not implemented.
};

#endif