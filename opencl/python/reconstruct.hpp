
extern numpy_boost<float, 3>
reconstruct_cgls(const numpy_boost<float, 3> &pixels,
		 const numpy_boost<float, 1> &angles,
		 double rotation_centre, int resolution,
		 int niterations);
extern void reconstruct_tvreg();
