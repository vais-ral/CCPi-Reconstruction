
extern numpy_boost<double, 3>
ring_artefacts_aml(const numpy_boost<double, 3> &pixels,
		   const float param_n, const float param_r,
		   const int num_series);
extern numpy_boost<float, 3>
reconstruct_cgls(const numpy_boost<float, 3> &pixels,
		 const numpy_boost<float, 1> &angles,
		 double rotation_centre, int resolution,
		 int niterations, int nthreads);
extern numpy_boost<float, 3>
reconstruct_sirt(const numpy_boost<float, 3> &pixels,
		 const numpy_boost<float, 1> &angles,
		 double rotation_centre, int resolution,
		 int niterations, int nthreads);
extern numpy_boost<float, 3>
reconstruct_mlem(const numpy_boost<float, 3> &pixels,
		 const numpy_boost<float, 1> &angles,
		 double rotation_centre, int resolution,
		 int niterations, int nthreads);
extern numpy_boost<float, 3>
reconstruct_cgls_tikhonov(const numpy_boost<float, 3> &pixels,
			  const numpy_boost<float, 1> &angles,
			  double rotation_centre, int resolution,
			  int niterations, int nthreads, double regularize);
extern numpy_boost<float, 3>
reconstruct_cgls_tvreg(const numpy_boost<float, 3> &pixels,
		       const numpy_boost<float, 1> &angles,
		       double rotation_centre, int resolution,
		       int niterations, int nthreads, double regularize);
extern void reconstruct_tvreg();
