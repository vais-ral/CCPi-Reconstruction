
/*
numpy_boost<float, 2> test(const numpy_boost<float, 3> &pixels,
			   const numpy_boost<float, 1> &angles,
			   const int niterations)
{
  std::cout << "Address " << (long(pixels.data()) % 64) << '\n';
  std::cout << "Size " << pixels.shape()[0] << ' ' << pixels.shape()[1]
	    << ' ' << pixels.shape()[2] << '\n';
  int dims[] = {5, 8};
  //dims[0] = 5;
  //dims[1] = 4;
  //dims[2] = 8;
  //numpy_boost<float, 1> voxels(dims);
  //PyObject *v = voxels.py_ptr();
  //Py_INCREF(v);
  numpy_boost<float, 2> v(dims);
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 8; j++)
      v[i][j] = i + j;
  return v;
}
*/

numpy_boost<float, 3> reconstruct_cgls(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
				       double rotation_centre, int resolution,
				       int niterations)
{
  // Todo ring artefacts choice etc.

  bool beam_harden = false;
  // vertical size to break data up into for processing
  const int blocking_factor = 0;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  CCPi::instrument *instrument = new CCPi::Diamond(true, 0.05, 2, true, 0.000,
						   0.0050, 1,
						   CCPi::ring_artefacts_aml);
  CCPi::reconstruction_alg *algorithm = new CCPi::cgls_3d(niterations);
  //if (blocking_factor > 0 and instrument->supports_blocks())
  //  recon_algorithm = new CCPi::cgls_2d(niterations, pixels_per_voxel);
  machine::initialise();
  // instrument setup from pixels/angles will probably copy
  voxel_data *voxels = reconstruct(instrument, algorithm, pixels, angles,
				   rotation_centre, resolution,
				   blocking_factor, beam_harden);
  machine::exit();
  delete algorithm;
  delete instrument;
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
  numpy_boost<float, 3> varray(dims);
  if (voxels == 0)
    varray[0][0][0] = 0.0;
  else {
    // Todo - vector? and remove buffered region
    for (int i = 0; i < dims[0]; i++)
      for (int j = 0; j < dims[1]; j++)
	for (int k = 0; k < dims[2]; k++)
	  varray[i][j][k] = (*voxels)[i][j][k];
    delete voxels;
  }
  return varray;
}

void reconstruct_tvreg()
{
}
