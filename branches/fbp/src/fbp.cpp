
#include <complex>
#include "base_types.hpp"
#include "fbp.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "fft.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

// Todo - not sure 1d.data() is legal
// - move these to header as will be needed elsewhere.
typedef std::vector<real> real_1d;
typedef boost::multi_array<real, 2> real_2d;
typedef boost::multi_array<real, 3> real_3d;
typedef std::vector<std::complex<real> > complex_1d;
typedef boost::multi_array<std::complex<real>, 2> complex_2d;
typedef std::vector<pixel_type> pixel_1d;
typedef boost::multi_array<pixel_type, 2> pixel_2d;
typedef boost::multi_array<pixel_type, 3> pixel_3d;
typedef boost::multi_array<std::complex<pixel_type>, 3> complex_pixel_3d;
typedef std::vector<voxel_type> voxel_1d;
typedef boost::multi_array<voxel_type, 2> voxel_2d;
typedef boost::multi_array<voxel_type, 3> voxel_3d;

static bool fbp_axial(complex_1d &vf, complex_1d &vi, real_1d &vo,
		      const real angle_step, const int n, const int nangles);
static void create_filter(real_1d &filter, const filter_name_t name,
			  const filter_window_t window,
			  const filter_norm_t norm, const int n,
			  const real bandwidth, const real pixel_size,
			  const real angle_step, const int nangles);
static void find_params(const int nangles, const int nx, const int ny,
			const int nh, /*int &mod_x, int &mod_y,*/ int &x_width,
			int &x_power, const real centre, real &roi_r,
			real &roi_a, int &required_width);
static void transform_sinogram_sol(const int nx, const int ny,
				   const int nangles, const int nh,
				   const int x_width, const real roi_r,
				   const real roi_a, const real centre,
				   real_2d &mapx, real_2d &mapy,
				   const real angle_step,
				   const int required_width, int &ta_nx,
				   int &ta_ny, int &tb_nx, int &tb_ny, int &nz,
				   const real alpha);
static void transform_sinogram(real_3d &mapping, const real *pixels,
			       const int z, const int nangles, const int nv,
			       const int nh, const real_2d &mapx,
			       const real_2d &mapy, int &ta_nx, int &ta_ny,
			       int &tb_nx, int &tb_ny, int &nz);

/*
void fbp_axial(real* dev_out, float da, unsigned int nx, unsigned int ny, float fnr){
	unsigned int ixx, i;
	float r, s, sum;
		
	ixx = blockIdx.x*blockDim.x + threadIdx.x;
	r = float(ixx)-fnr;
	sum = 0.0f;
			
	for(i = 0; i<ny; i++){
		s = fnr + r*cosf(da*float(i));
		sum += tex1Dfetch(texA,s);
	}

	dev_out[ixx] = sum;
}
*/

// parallel_gpu.cu
bool fbp_axial(complex_1d &vf, complex_1d &vi, real_1d &vo,
	       const real angle_step, const int n, const int nangles)
{
  // angle step looks like 1 degree for total of 180 in pi, vs 192 nangles?
  real fnr = real(0.5) * real(n - 1);
  complex_1d &dev_filter = vf;
  complex_1d &dev_data = vi;
  fft_1d_forward(dev_data.data(), n, 1);
  for (int i = 0; i < n; i++)
    dev_data[i] *= dev_filter[i];
  fft_1d_inverse(dev_data.data(), n, 1);
  real_1d dev_in(n);
  for (int i = 0; i < n; i++)
    dev_in[i] = dev_data[i].real();
  // bind dev_in to 1d linear clamped texture
  real_1d &dev_out = vo;
  // CUDA fbp_axial call
  for (int i = 0; i < n; i++) {
    // Check that nx/ny are the same, or was ta_nx a bit bigger to avoid issues
    // at end/start of array?
    // fnr = (n-1)/2 so i from 0 to n-1 -> r from -(n-1)/2 to (n-1)/2
    real r = real(i) - fnr;
    real sum = 0.0;
    for (int j = 0; j < nangles; j++) {
      // cosine ranges from -1/+1 so (n-1)/2 +/- (n-1)/2 is 0 to n-1
      // so should be a safe mapping into the array
      real s = fnr + r * std::cos(angle_step * float(j));
      int si = int(s);
      // Todo - improve this, do we handle the end points correctly?
      if (si >= 0 and si < n - 1) {
	real diff = s - real(si);
	sum += (real(1.0) - diff) * dev_in[si] + diff * dev_in[si + 1];
      } else if (si < 0)
	sum += dev_in[0];
      else
	sum += dev_in[n - 1];
    }
    dev_out[i] = sum;
  }
  return true;
}

// diamond_cpp - form_filter
void create_filter(real_1d &filter, const filter_name_t name,
		   const filter_window_t window,
		   const filter_norm_t norm, const int n,
		   const real bandwidth, const real pixel_size,
		   const real angle_step, const int nangles)
{
  real pif = M_PI;
  {
  real_1d v1(n);
  real_1d v2(n);
  real_1d v3(n);
  real_1d v4(n);
  /*
  if (v1 == 0 or v2 == 0 or v3 == 0 or v4 == 0 or filter == 0) {
    std::cerr << "can not allocate memory (filter)\n";
    return 0;
  }
  */
  real omega = bandwidth;
  if (omega > 1.0 or omega < 0.00001) {
    std::cerr << "Bandwidth wrong (should be in (0, 1])\n";
    return;
  }
  for (int i = 0; i < n; i++)
    v1[i] = real(i);
  for (int i = 0; i < n; i++) {
    if (v1[i] > real(n) * omega)
      v4[i] = 0.0;
    else
      v4[i] = v1[i];
  }
  //ippsThreshold_GTVal_32f(v1,v4,n,float(n)*omega,0.0);
  switch (name) {
  case Ram_Lak_filter:
    for (int i = 0; i < n; i++)
      v1[i] = v4[i];
    break;
  case Shepp_Logan_filter:
    for (int i = 0; i < n; i++)
      v1[i] *= pif / (real(n) * omega);
    for (int i = 0; i < n; i++)
      if (v1[i] > real(0.5) * pif)
	v1[i] = 0.0;
    // from gen_functions
    for (int i = 0; i < n; i++) {
      if (std::abs(v1[i]) < 1.e-5) {
	real y = v1[i] * v1[i];
	v2[i] = (real(1.0) - y / real(6.0) + y * y / real(120.0));
      } else
	v2[i] = std::sin(v1[i]) / v1[i];
    }
    for (int i = 0; i < n; i++)
      v1[i] = v2[i] * v4[i];
    break;
  case Cosine_filter:
    for (int i = 0; i < n; i++)
      v1[i] *= pif / (real(n) * omega);
    for (int i = 0; i < n; i++)
      if (v1[i] > real(0.5) * pif)
	v1[i] = 0.0;
    for (int i = 0; i < n; i++)
      v2[i] = std::cos(v1[i]);
    for (int i = 0; i < n; i++)
      v1[i] = v2[i] * v4[i];
    break;
  default:
    std::cerr << "Unknown filter - abort\n";
    return;
    break;
  }
  for (int i = 0; i < n; i++)
    v4[i] = v1[i];
  for (int i = 0; i < n; i++)
    v1[i] = real(i);
  for (int i = 0; i < n; i++)
    v1[i] *= real(2.0) * pif / real(n);
  for (int i = 0; i < n; i++)
    v2[i] = v1[i] * real(2.0);
  for (int i = 0; i < n; i++)
    v3[i] = std::cos(v1[i]);
  for (int i = 0; i < n; i++)
    v1[i] = v3[i];
  for (int i = 0; i < n; i++)
    v3[i] = std::cos(v2[i]);
  switch (window) {
  case No_window:
    for (int i = 0; i < n; i++)
      v3[i] = 1.0;
    break;
  case Hann_window:
    for (int i = 0; i < n; i++)
      v1[i] *= real(0.5);
    for (int i = 0; i < n; i++)
      v3[i] = v1[i] + real(0.5);
    break;
  case Hamming_window:
    for (int i = 0; i < n; i++)
      v1[i] *= real(0.46);
    for (int i = 0; i < n; i++)
      v3[i] = v1[i] + real(0.54);
    break;
  case Blackman_window:
    for (int i = 0; i < n; i++)
      v3[i] *= real(0.08);
    for (int i = 0; i < n; i++)
      v3[i] += real(0.42);
    for (int i = 0; i < n; i++)
      v3[i] += v1[i] * real(0.5);
    break;
  default:
    std::cerr << "Unknown window - abort\n";
    return;
    break;
  }
  for (int i = 0; i < n; i++)
    v4[i] *= v3[i];
  int n1 = (n + 1) / 2;
  for (int i = 0; i < n1; i++)
    filter[i] = v4[i];
  int n2 = n - n1;
  for (int i = 0; i < n2; i++)
    filter[n1 + i] = v4[n2 - i];
  //ippsFlip_32f(v4+1,xd->vfilter+n1,n2);
  }
  real_1d vecti(n);
  real_1d vecto(n);
  real_1d vectr(n);
  real fxc = real(n - 1) * real(0.5);
  real frad = real(0.1) * real(n - 1);
  for (int i = 0; i < n; i++)
    vecti[i] = -fxc + real(i);
  for (int i = 0; i < n; i++)
    vecti[i] = std::abs(vecti[i]);
  for (int i = 0; i < n; i++)
    vectr[i] = vecti[i];
  for (int i = 0; i < n; i++)
    vecto[i] = vecti[i];
  for (int i = 0; i < n; i++)
    vecti[i] = vecti[i] * vecti[i];
  for (int i = 0; i < n; i++)
    vecti[i] = frad * frad - vecti[i];
  for (int i = 0; i < n; i++)
    if (vecti[i] < real(0.0))
      vecti[i] = 0.0;
  for (int i = 0; i < n; i++)
    vecti[i] = std::sqrt(vecti[i]);
  switch (norm) {
  case CC_norm:
  case CA_norm:
    for (int i = 0; i < n; i++) {
      if (vectr[i] < frad)
	vectr[i] = 1.0;
      else if (vectr[i] > frad)
	vectr[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
      vecti[i] *= real(2.0);
    break;
  case PC_norm:
  case PA_norm:
    for (int i = 0; i < n; i++)
      vectr[i] = frad - vectr[i];
    for (int i = 0; i < n; i++)
      if (vectr[i] < real(0.0))
	vectr[i] = 0.0;
    for (int i = 0; i < n; i++) {
      real x = vecto[i];
      real gamma = x / frad;
      real sg = real(1.0) - gamma * gamma;
      if (sg <= real(0.0))
	vecti[i] = 0.0;
      else {
	sg = std::sqrt(sg);
	vecti[i] = frad * frad * (sg - gamma * gamma
				  * std::log((real(1.0) + sg) / gamma));
      }
    }
    break;
  case SC_norm:
  case SA_norm:
    for (int i = 0; i < n; i++)
      vectr[i] = vecti[i];
    for (int i = 0; i < n; i++)
      vecti[i] = vecti[i] * vecti[i];
    for (int i = 0; i < n; i++)
      vecti[i] *= real(0.5) * pif;
    break;
  default:
    std::cerr << "Unknown sample for normalisation - abort\n";
    return;
    break;
  }
  {
  complex_1d veccF(n);
  complex_1d veccI(n);
  for (int i = 0; i < n; i++)
    veccF[i] = filter[i];
  for (int i = 0; i < n; i++)
    veccI[i] = vecti[i];
  if (!fbp_axial(veccF, veccI, vecto, angle_step, n, nangles))
    return;
  }
  int ns = int(real(0.47) * real(n));
  int nl = int(real(0.06) * real(n));
  real sr = 0.0;
  real so = 0.0;
  switch (norm) {
  case CC_norm:
  case PC_norm:
  case SC_norm:
    for (int i = 0; i < nl; i++)
      sr = vectr[ns + i];
    sr /= real(nl);
    for (int i = 0; i < nl; i++)
      so = vecto[ns + i];
    so /= real(nl);
    break;
  case CA_norm:
  case PA_norm:
  case SA_norm:
    for (int i = 0; i < n; i++)
      sr = vectr[i];
    sr /= real(n);
    for (int i = 0; i < n; i++)
      so = vecto[i];
    so /= real(n);
    break;
  default:
    break;
  }
  real su = sr / (so * pixel_size);
  for (int i = 0; i < n; i++)
    filter[i] *= su;
}

// diamond_cpp
void find_params(const int nangles, const int nx, const int ny, const int nh,
		 /*int &mod_x, int &mod_y,*/ int &x_width, int &x_power,
		 const real centre, real &roi_r, real &roi_a,
		 int &required_width)
{
  // block sizes from GPU, CPU can probably do differently
  // Todo - infinite blocks would be new_x = nx, new_y = ny
  const int x_block = 16;
  const int y_block = 12;
  int new_x = nx / x_block;
  if (nx % x_block != 0)
    new_x++;
  new_x *= x_block;
  int new_y = ny / y_block;
  if (ny % y_block != 0)
    new_y++;
  new_y *= y_block;
  // Todo?
  real out_pixel_size = 1.0; // pixels_per_voxel * voxel_size[0];
  //real width = new_x * out_pixel_size;
  //real height = new_y * out_pixel_size;
  //real radius = 0.5 * std::sqrt(width * width + height * height);
  real image_centre = centre; //real(nh - 1) / 2.0; // was xc
  real roi_orig_x = out_pixel_size * real(0.5) * new_x;
  real roi_orig_y = out_pixel_size * real(0.5) * new_y;
  real x1 = roi_orig_x - image_centre;
  real y1 = roi_orig_y - image_centre;
  // would use angle_roi here if had one
  roi_orig_x = x1 + image_centre;
  roi_orig_y = y1 + image_centre;
  real new_width = new_x * out_pixel_size;
  real new_height = new_y * out_pixel_size;
  real new_radius = real(0.5) * std::sqrt(new_width * new_width
					  + new_height * new_height);
  roi_r = std::sqrt(x1 * x1 + y1 * y1);
  roi_a = std::atan2(y1, x1);
  //real x_off = roi_orig_x - image_centre;
  //real y_off = roi_orig_y - image_centre;
  int x_min1 = (int)std::floor(image_centre - new_radius) - 1;
  int x_max1 = (int)std::ceil(image_centre + new_radius) + 1;
  // No cropping and shape is point, as input, so
  int x_min2 = 0;
  int x_max2 = nx;
  // GPU code has error here in assigning these to unsigned!
  int x_min = std::min(x_min1, x_min2);
  int x_max = std::max(x_max1, x_max2);
  required_width = x_max - x_min;
  // find power of 2 to build this
  //real exponent = std::log2(real(required_width));
  int extpower = 0;
  int next = 0;
  {
    int u = 0;
    int w = 1;
    while (w < required_width) {
      u++;
      w *= 2;
    }
    next = w;
    extpower = u;
  }
  //mod_x = new_x;
  //mod_y = new_y;
  x_width = next;
  x_power = extpower;
}

// transform_sinogram.cpp
void transform_sinogram_sol(const int nx, const int ny, const int nangles,
			    const int nh, const int x_width, const real roi_r,
			    const real roi_a, const real centre, real_2d &mapx,
			    real_2d &mapy, const real angle_step,
			    const int required_width, int &ta_nx, int &ta_ny,
			    int &tb_nx, int &tb_ny, int &nz, const real alpha)
{
  if (alpha < 1e-6)
    std::cerr << "PixelParam too small\n";
  // basic check_roi values
  //real xmin = 0.0;
  //real xmax = real(nx - 1);
  //real ymin = 0.0;
  //real ymax = real(ny - 1);
  // no missed proj so
  ta_ny = nangles; // with angle step as expected
  ta_nx = x_width; // from find params.
  tb_nx = nh + 6; // or new_x?
  tb_ny = nangles + 3;
  // pixel shape
  //tb_ny++;
  // tb_nx/tb_ny stored
  //int ta_ny_13 = ta_ny;
  // ta_ny changed to value required for num y threads
  // SinCosChunk nom(192, warpsize) -> 192 -> /2 -> 96 (gi->vertH)
  const int y_block = 96;
  int new_y = ta_ny / y_block;
  if (ta_ny % y_block != 0)
    new_y++;
  ta_ny = new_y * y_block;
  // ta_nx/ta_ny -> gi->nxo, nyo and aD_.ta_nx, ta_ny
  nz = 4;
  // allocate vtb tb_nx * tb_ny, mapx+mapy ta_nx * ta_ny
  // and vta ta_nx * ta_ny * nz
  // position old data in new region
  int extra = x_width - required_width;
  int extra_left = extra / 2;
  int extra_right = extra - extra_left;
  real xc = centre;
  real fxc = real(extra_left) + xc;
  int warpsize = 32;
  int ixc = (int)std::ceil(fxc / warpsize) * warpsize;
  real dxc = real(ixc) - fxc;
  if (dxc > real(extra_right))
    dxc = 0.0;
  real new_xc = real(extra_left) + xc + dxc; //124.5 -> 224 warp=32?
  dxc = new_xc - xc;

  mapx.resize(boost::extents[ta_nx][ta_ny]);
  for (int i = 0; i < ta_nx; i++)
    mapx[i][0] = real(i);
  for (int j = 1; j < ta_ny; j++)
    for (int i = 0; i < ta_nx; i++)
      mapx[i][j] = mapx[i][0];
  for (int j = 0; j < ta_ny; j++) {
    real shift = dxc - roi_r * std::cos(roi_a - angle_step * real(j))
      - real(3.0);
    for (int i = 0; i < ta_nx; i++)
      mapx[i][j] -= shift;
  }
  mapy.resize(boost::extents[ta_nx][ta_ny]);
  for (int j = 0; j < ta_ny; j++)
    for (int i = 0; i < ta_nx; i++)
      mapy[i][j] = real(j);
  float cropbottom = 0.0; // ?
  for (int j = 0; j < ta_ny; j++)
    for (int i = 0; i < ta_nx; i++)
      mapy[i][j] += cropbottom;
  int croptop = 0;
  int mp = 0; // missing angles
  for (int j = 0; j < ta_ny; j++) {
    for (int i = 0; i < ta_nx; i++) {
      if (mapy[i][j] > real(nangles + mp - 1 - croptop))
	mapy[i][j] = real(tb_ny - 1);
    }
  }
  int cropleft = 0;
  for (int j = 0; j < ta_ny; j++) {
    for (int i = 0; i < ta_nx; i++) {
      if (mapx[i][j] < real(cropleft + 3))
	mapx[i][j] = 1.0;
    }
  }
  int cropright = 0;
  for (int j = 0; j < ta_ny; j++) {
    for (int i = 0; i < ta_nx; i++) {
      if (mapx[i][j] > real(nh + 2 - cropright))
	mapx[i][j] = real(tb_nx) - real(2.0);
    }
  }
}

void transform_sinogram(real_3d &mapping, const real *pixels, const int z,
			const int nangles, const int nv, const int nh,
			const real_2d &mapx, const real_2d &mapy, int &ta_nx,
			int &ta_ny, int &tb_nx, int &tb_ny, int &nz)
{
  // mapping is vta
  pixel_2d vtb(boost::extents[tb_nx][tb_ny], boost::fortran_storage_order());
  for (int i = 0; i < nz; i++)
    for (int j = 0; j < ta_ny; j++)
      for (int k = 0; k < ta_nx; k++)
	mapping[k][j][i] = 0.0;
  for (int j = 0; j < tb_ny; j++)
    for (int i = 0; i < tb_nx; i++)
      vtb[i][j] = 0.0;
  // This is for shape point with no missing angles
  for (int j = 0; j < nangles; j++) {
    for (int i = 0; i < nh; i++)
      vtb[i + 3][j] = pixels[j * nh * nv + z * nh + i];
  }
  for (int j = 0; j < nangles; j++) {
    std::cerr << "loop " << j << '\n';
    for (int i = 0; i < nh; i++)
      std::cerr << j * nh + i << ' ' << pixels[j * nh * nv + z * nh + i] << '\n';
  }
  // no missed stuff so skip next if
  int extrapolation_pixels = 10; // Todo
  int crop_left = 0;
  int crop_right = 0;
  for (int i = 0; i < nangles; i++) {
    real sv = 0.0;
    for (int j = 0; j < extrapolation_pixels; j++)
      sv += vtb[3 + crop_left + j][i];
    sv /= real(extrapolation_pixels);
    for (int j = 0; j < 3; j++)
      vtb[j][i] = sv;
    sv = 0.0;
    for (int j = 0; j < extrapolation_pixels; j++)
      sv += vtb[nh - extrapolation_pixels + 3 - crop_right + j][i];
    sv /= real(extrapolation_pixels);
    for (int j = 0; j < 3; j++)
      vtb[tb_nx - 3 + j][i] = sv;
  }
  for (int j = 0; j < tb_ny; j++)
    for (int i = 0; i < tb_nx; i++)
      std::cerr << j * tb_nx + i << ' ' << vtb[i][j] << '\n';
  // remap - pixel order, uses interpolation type from option in transform
  // section of xml - Nearest-Neighbour, Linear, or Cubic. This is linear
  const real tol = 1e-8;
  // map arrays seem to be [ny][nx], vtb is [y][x] also
  for (int j = 0; j < ta_ny; j++) {
    for (int i = 0; i < ta_nx; i++) {
      //mapping[j][i] = vtb[mapy[j][i]][mapx[j][i]];
      // Todo - we can effectively pre-map all this surely since ta_nx etc
      // are known when we generate mapx/mapy
      if (mapx[i][j] < real(0.0) or mapx[i][j] > real(tb_nx - 1) or
	  mapy[i][j] < real(0.0) or mapy[i][j] > real(tb_ny - 1))
	; // no edge smoothing was requested
      else {
	int ix = int(mapx[i][j]);
	int iy = int(mapy[i][j]);
	real diffx = mapx[i][j] - real(ix);
	real diffy = mapy[i][j] - real(iy);
	// Todo - the 1D interp here assumes integer steps of 1.0
	if (diffx < tol) {
	  if (diffy < tol)
	    mapping[i][j][0] = vtb[ix][iy];
	  else {
	    real f1 = vtb[ix][iy];
	    real f2 = vtb[ix][iy + 1];
	    real y = mapy[i][j];
	    real y1 = real(iy);
	    real y2 = real(iy + 1);
	    mapping[i][j][0] = (f1 * (y2 - y) + f2  * (y - y1)) / (y2 - y1);
	  }
	} else if (diffy < tol) {
	  real f1 = vtb[ix][iy];
	  real f2 = vtb[ix + 1][iy];
	  real x = mapx[i][j];
	  real x1 = real(ix);
	  // if x is 3.25 then we are interpolating vtb[3] to vtb[4]
	  // so the step is 1.0 surely or do we have a separate set of coords
	  // for vtb somewhere?
	  real x2 = real(ix + 1);
	  mapping[i][j][0] = (f1 * (x2 - x) + f2  * (x - x1)) / (x2 - x1);
	} else {
	  // 2D interpolation required - is this right?
	  real x = mapx[i][j];
	  real x1 = real(ix);
	  real x2 = real(ix + 1);
	  real y = mapy[i][j];
	  real y1 = real(iy);
	  real y2 = real(iy + 1);
	  real f11 = vtb[ix][iy];
	  real f12 = vtb[ix][iy + 1];
	  real f21 = vtb[ix + 1][iy];
	  real f22 = vtb[ix + 1][iy + 1];
	  mapping[i][j][0] = (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1)
			      * (y2 - y) + f12 * (x2 - x) * (y - y1)
			      + f22 * (x - x1) * (y - y1)) /
	    ((x2 - x1) * (y2 - y1));
	}
      }
    }
  }
  // replicate across nz
  for (int i = 1; i < nz; i++) {
    for (int j = 0; j < ta_ny; j++)
      for (int k = 0; k < ta_nx; k++)
	mapping[k][j][i] = mapping[k][j][0];
  }
}

bool CCPi::fbp_reconstruction(const instrument *device, voxel_data &voxels,
			      const real origin[3], const real voxel_size[3],
			      const filter_name_t name,
			      const filter_window_t window,
			      const filter_norm_t norm, const real bandwidth)
{
  bool ok = device->filtered_back_project(voxels, origin, voxel_size,
					  name, window, norm, bandwidth);
  return ok;
}

bool CCPi::parallel_beam::filtered_back_project(voxel_data &voxels,
						const real origin[3],
						const real voxel_size[3],
						const filter_name_t name,
						const filter_window_t window,
						const filter_norm_t norm,
						const real bandwidth) const
{
  // Todo?
  real alpha = 1.0; // pixelParam
  const voxel_data::size_type *sz = voxels.shape();
  //sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  //voxel_type *const data = voxels.data();
  pixel_type *const pixels = get_pixel_data();
  //sl_int n_rays = get_data_size();

  int nv = get_num_v_pixels();
  int npixels_per_voxel = nv / sz[2];
  if (npixels_per_voxel != 1) {
    std::cerr << "FBP error\n";
    return false;
  }
  // scan z slices separately
  int z_vox = 0;
  int z_count = 0;
  int nangles = get_num_angles();
  real *angles = get_phi();
  real angle_step = angles[1] - angles[0];
  int nh = get_num_h_pixels();
  real *h_pixels = get_h_pixels();
  real h_step = (h_pixels[1] - h_pixels[0]);
  real inv_h_step = real(1.0) / h_step;
  //int new_x = 0;
  //int new_y = 0;
  int x_width = 0;
  int x_power = 0;
  int required_width = 0;
  real roi_r;
  real roi_a;
  real centre = real(nh - 1) / real(2.0); // Todo - read this or find?
  find_params(nangles, int(sz[0]), int(sz[1]), nh, /*new_x, new_y,*/
	      x_width, x_power, centre, roi_r, roi_a, required_width);
  timer fbptime(USE_TIMER);
  // Todo - save filter between block calls?
  //std::cerr << "width " << x_width << ' ' << new_x << ' ' << new_y << '\n';
  real_1d filter(x_width);
  create_filter(filter, name, window, norm, x_width, bandwidth,
		0.5 /*h_step*/, angle_step, get_num_angles());
  fbptime.accumulate();
  fbptime.output("create filter");
  // centres of voxels
  std::vector<real> xdata(sz[0]);
  for (int i = 0; i < (int)sz[0]; i++)
    xdata[i] = origin[0] + (real(i) + real(0.5)) * voxel_size[0];
  std::vector<real> ydata(sz[1]);
  for (int i = 0; i < (int)sz[1]; i++)
    ydata[i] = origin[1] + (real(i) + real(0.5)) * voxel_size[1];
  fbptime.reset();
  std::cerr << nangles << ' ' << nv << ' ' << nh << '\n';
  int ta_nx = nh;
  int ta_ny = nangles;
  int tb_nx = nh;
  int tb_ny = nangles;
  int nz = 0;
  // Yuck, pity resize can't specify storage order
  real_2d mapx(boost::extents[1][1], boost::fortran_storage_order());
  real_2d mapy(boost::extents[1][1], boost::fortran_storage_order());
  bool first = true;
  //for (int z = 0; z < nv; z++) {
  for (int z = 49; z < 51; z++) {
    std::cerr << "whats going on with slice[z]?\n";
    if (first)
      transform_sinogram_sol(sz[0], sz[1], nangles, nh, x_width, roi_r, roi_a,
			     centre, mapx, mapy, angle_step, required_width,
			     ta_nx, ta_ny, tb_nx, tb_ny, nz, alpha);
    first = false;
    complex_pixel_3d slice(boost::extents[ta_nx][ta_ny][nz],
			   boost::fortran_storage_order());
    {
      pixel_3d mapping(boost::extents[ta_nx][ta_ny][nz],
		       boost::fortran_storage_order());
      transform_sinogram(mapping, pixels, z, nangles, nv, nh, mapx, mapy,
			 ta_nx, ta_ny, tb_nx, tb_ny, nz);
      for (int ak = 0; ak < nz; ak++)
	for (int aj = 0; aj < ta_ny; aj++)
	  for (int ai = 0; ai < ta_nx; ai++)
	    slice[ai][aj][ak] = mapping[ai][aj][ak];
    }
    if (z == 49) {
      std::cerr << "slice\n";
      for (int iy = 0; iy < ta_ny; iy++) {
	for (int ix = 0; ix < ta_nx; ix++) {
	  std::cerr << ix << ' ' << slice[ix][iy][0] << '\n';
	}
      }
    }
    // Todo - is this really complex?
    fft_1d_forward(slice.data(), ta_nx, ta_ny);
    // apply filter
    for (int a = 0; a < ta_ny; a++) {
      for (int h = 0; h < ta_nx; h++)
	slice[h][a][0] *= filter[h];
    }
    fft_1d_inverse(slice.data(), ta_nx, ta_ny);
    if (z == 49) {
      std::cerr << "transform\n";
      for (int iy = 0; iy < ta_ny; iy++) {
	for (int ix = 0; ix < ta_nx; ix++) {
	  std::cerr << ix << ' ' << slice[ix][iy][0] << '\n';
	}
      }
    }
    // Only use real part (or convert?)
    // project into voxel slice - based on fbp_std2
    // Todo - should these be centred in the voxel?, probably although
    // related to projection if multiple pixels
    std::cerr << "Should these be new_x/new_y?\n";
    for (int y = 0; y < (int)sz[1]; y++) {
      real ypos = ydata[y];
      for (int x = 0; x < (int)sz[0]; x++) {
	real xpos = xdata[x];
	// Todo loop number of pixels per voxel? in h <=> x
	real sum = 0.0;
	for (int a = 0; a < ta_ny; a++) {
	  // find rotated position in slice, gpu interp might help here.
	  // rotates into theta = 0 plane to find position in horizontal
	  // pixels, interpolate to relate to centre of voxel.
	  // ?? y' should be in h_pixel, x' map into filter? Todo
	  // this is x, so seems wrong
	  //std::cerr << "Check this mapping\n";
	  //real p = (xpos * cos(angles[a]) + ypos * sin(angles[a]));
	  real p = (- xpos * sin(angles[a]) + ypos * cos(angles[a]));
	  // y' = -x sin + y cos
	  // Can I do a lot of this mapping once?
	  // Todo - expand h_pixels to guarantee its in the array and avoid if
	  // Todo - what is ta_nx_pixels[0]?
	  int i = int((p - h_pixels[0]) * inv_h_step);
	  if (i >= 0 and i < ta_nx - 1) { // Todo - end point < nh? or skip
	    real diff = (p - h_pixels[i]) * inv_h_step;
	    sum += ((real(1.0) - diff) * slice[i][a][0].real()
		    + diff * slice[i + 1][a][0].real());
	  }
	}
	voxels[x][y][z] += sum;
      }
    }
    z_count++;
    if (z_count == npixels_per_voxel) {
      z_vox++;
      z_count = 0;
    }
  }
  fbptime.accumulate();
  fbptime.output("filtered back projection");
  return true;
}

bool CCPi::cone_beam::filtered_back_project(voxel_data &voxels,
					    const real origin[3],
					    const real voxel_size[3],
					    const filter_name_t name,
					    const filter_window_t window,
					    const filter_norm_t norm,
					    const real bandwidth) const
{
  std::cerr << "FBP not implemented for cone beam\n";
  return false;
}
