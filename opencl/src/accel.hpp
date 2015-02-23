
#ifndef OCL_WRAPPER
#define OCL_WRAPPER

typedef void *dev_ptr;

namespace machine {

  void init_accelerator();
  void close_accelerator();
  bool has_accelerator();
  int number_of_accelerators();
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
  void device_free(dev_ptr data, const int device);
  bool check_accelerator_kernel(const char name[], const int device);
  void run_parallel_ah(const char name[], dev_ptr pix_buf, dev_ptr vox_buf,
		       const int vox_offset, dev_ptr xy_buff,
		       dev_ptr xy_offsets, const int n, const int nv,
		       const int nz, const int size, const int device);
  void run_parallel_xy(const char name[], dev_ptr pix_buf, const int offset,
		       dev_ptr vox_buf, dev_ptr xy_buff,
		       dev_ptr xy_offsets, const int n, const int nv,
		       const int nz, const int size, const int device);
  void accelerator_barrier(const int device);

}

#endif // OCL_WRAPPER
