
#ifndef OCL_WRAPPER
#define OCL_WRAPPER

#ifdef USE_OPENCL
#  include "cl.hpp"
typedef cl::Buffer dev_ptr;
typedef cl::Event event_t;
#else
typedef void *dev_ptr;
typedef void event_t;
#endif // OPENCL

namespace machine {

  void init_accelerator();
  void close_accelerator();
  bool has_accelerator();
  int number_of_accelerators();
  sl_int available_mem(const int device);
  sl_int largest_alloc(const int device);
  dev_ptr device_allocate(const sl_int size, const bool read_only,
			  const int device);
  void copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
		      const int device);
  void copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
		      const sl_int size, const int device);
  void copy_from_device(dev_ptr buffer, void *ptr, const sl_int size,
			const int device);
  void copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			const sl_int size, const int device);
  void copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
		      const int device, event_t *event);
  void copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
		      const sl_int size, const int device, event_t *event);
  void copy_from_device(dev_ptr buffer, void *ptr, const sl_int size,
			const int device, event_t *event);
  void copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			const sl_int size, const int device, event_t *event);
  void copy_to_device(void *ptr, dev_ptr buffer, const sl_int y_host,
		      const sl_int z_host, const sl_int x_size,
		      const sl_int y_size, const int z_size,
		      const int device, event_t *event);
  void copy_from_device(void *ptr, dev_ptr buffer, const sl_int y_host,
			const sl_int z_host, const sl_int x_size,
			const sl_int y_size, const int z_size,
			const int device);
  void device_free(dev_ptr data, const int device);
  bool check_accelerator_kernel(const char name[], const int device);
  void run_parallel_ah(const char name[], dev_ptr pix_buf,
		       dev_ptr vox_buf, const int vox_offset,
		       dev_ptr xy_buff, dev_ptr xy_offsets,
		       dev_ptr nlengths, const int nv,
		       const int nz, const int xy_size, const int ix,
		       const int size, const int dim3, const int device,
		       std::vector<event_t> *events);
  void run_parallel_xy(const char name[], dev_ptr pix_buf, dev_ptr vox_buf,
		       dev_ptr xy_buff, dev_ptr xy_offsets, dev_ptr h,
		       dev_ptr lengths, const int nv, const int nz,
		       const int start, const int ah_size, const int size,
		       const int dim3, const int device,
		       std::vector<event_t> *events);
  void run_cone_xy(const char name[], dev_ptr pix_buf, dev_ptr vox_buf,
		   dev_ptr xy_buff, dev_ptr xy_offsets, dev_ptr h,
		   dev_ptr lengths, const int nv, const int nz, const int start,
		   const int ah_size, const float pzbz, const int midp,
		   dev_ptr delta_z, dev_ptr inv_delz, dev_ptr vox_z,
		   const int size, const int dim3, const int device,
		   std::vector<event_t> *events);
  void run_cone_ah(const char name[], dev_ptr pix_buf, dev_ptr vox_buf,
		   dev_ptr xy0_buff, dev_ptr xy1_buff, dev_ptr xy_offsets,
		   dev_ptr h, dev_ptr lengths, const int nv, const int nz,
		   const int start, const int xy_size, const float pzbz,
		   const int midp, dev_ptr delta_z, dev_ptr inv_delz,
		   dev_ptr vox_z, const int size, const int dim3,
		   const int device, std::vector<event_t> *events);
  void accelerator_barrier(const int device);
  void accelerator_flush(const int device);
  void accelerator_complete(const int device);

}

#endif // OCL_WRAPPER
