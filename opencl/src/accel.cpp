
#ifdef MATLAB_MEX_FILE
#  ifdef WIN32
#    define snprintf _snprintf
#  endif // WIN32
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // MEX_FILE
#include "accel.hpp"
#include "ui_calls.hpp"
#include <iostream>
#include <fstream>

#ifdef USE_OPENCL

#include "cl.hpp"

namespace machine {

  bool accelerator = false;
  int num_devices = 0;
  int max_work0 = 1000000;
  int max_work01 = 10000000;
  int max_work2 = 1000000;
  cl::Context context;
  cl::Program program;
  cl::Kernel kernel;
  std::vector<cl::CommandQueue *> queue;

  sl_int max_alloc = 1000000000000;
  sl_int max_mem = 1000000000000;

}

static char *get_file(const char name[])
{
  // Open file stream
  std::fstream f(name, (std::fstream::in | std::fstream::binary));

  // Check if we have opened file stream
  if (f.is_open()) {
    size_t size;
    // Find the stream size
    f.seekg(0, std::fstream::end);
    size = f.tellg();
    f.seekg(0, std::fstream::beg);

    char *str = new char[size + 1];
    if (!str) {
      f.close();
      return 0;
    }

    // Read file
    f.read(str, size);
    f.close();
    str[size] = '\0';

    return str;
  }
  return 0;
}

void machine::init_accelerator()
{
  bool ok = true;
  cl_int err;

  std::vector<cl::Platform> platforms;
  //std::cout<<"Me!\nGetting Platform Information\n";
  err = cl::Platform::get(&platforms);
  if (err != CL_SUCCESS) {
    //std::cerr << "Platform::get() failed (" << err << ")" << std::endl;
    ok = false;
  }

  if (ok) {
    int plat_num = -1;
    if (platforms.size() > 0) {
      int j = 0;
      std::vector<cl::Platform>::iterator i;
      for (i = platforms.begin(); i != platforms.end(); ++i) {
	//std::cout << "Platform - ";
	//std::cout << (*i).getInfo<CL_PLATFORM_NAME>(&err).c_str() << '\n';
	if (strncmp((*i).getInfo<CL_PLATFORM_VERSION>(&err).c_str(),
		    "OpenCL 1.1", 10) == 0 and plat_num == -1)
	  plat_num = j;
	j++;
      }
      if (plat_num == -1)
	plat_num = 0;
    } else {
      report_error("No platforms");
      ok = false;
    }

    cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM,
		    (cl_context_properties)(platforms[plat_num])(),
				     0 };

    if (ok) {
      //std::cout<<"Creating a context\n";
      context = cl::Context(CL_DEVICE_TYPE_GPU, cps, NULL, NULL, &err);
      if (err != CL_SUCCESS) {
	//std::cerr << "Context::Context() failed (" << err << ")\n";
	ok = false;
      }
    }

    if (ok) {
      //std::cout<<"Getting device info\n";
      std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
      if (err != CL_SUCCESS) {
	report_error("Context::getInfo() failed - ", err);
	ok = false;
      }
      if (ok) {
	if (devices.size() == 0) {
	  report_error("No device available");
	  ok = false;
	} else {
	  num_devices = (int)devices.size();
	  queue.resize(num_devices);
	  int k = 0;
	  for (int i = 0; i < num_devices; i++) {
	    queue[k] = 0;
	    //std::cout << devices.size() << " devices available\n";
	    cl_uint val;
	    err = devices[i].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &val);
	    if (err != CL_SUCCESS)
	      report_error("Compute units failed ", err);
	    else {
	      add_output("Max compute units = ");
	      add_output(val);
	      send_output();
	    }
	    err = devices[i].getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &val);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Work item dims failed " << err << '\n'*/;
	    else {
	      if (val < 2)
		report_error("Invalid assumption about work dims >= 2");
	      else {
		add_output("Max work item dims = ");
		add_output(val);
		send_output();
	      }
	    }
	    size_t wval[3];
	    err = clGetDeviceInfo(devices[i](), CL_DEVICE_MAX_WORK_ITEM_SIZES,
				  sizeof(wval), wval, NULL);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Work item sizes failed " << err << '\n'*/;
	    else {
	      add_output("Max work item sizes =");
	      for (int z = 0; z < (int)val; z++) {
		add_output(" ");
		add_output(wval[z]);
	      }
	      send_output();
	      if ((int)wval[0] < max_work0)
		max_work0 = wval[0];
	      if (val == 3) {
		if (int(wval[0] * wval[2]) < max_work01)
		  max_work01 = wval[0] * wval[2];
	      } else {
		max_work01 = max_work0;
	      }
	      if (int(wval[1]) < max_work2)
		max_work2 = wval[1];
	    }
	    if (val < 3)
	      report_error("Only 2 work dims supported");
	    size_t tval;
	    err = devices[i].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &tval);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Work group size failed " << err << '\n'*/;
	    else {
	      add_output("Max work group size = ");
	      add_output(tval);
	      send_output();
	    }
	    cl_ulong lval;
	    err = devices[i].getInfo(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, &lval);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Constant buffer size failed " << err << '\n'*/;
	    else {
	      add_output("Max constant buffer size = ");
	      add_output(lval);
	      send_output();
	    }
	    err = devices[i].getInfo(CL_DEVICE_MAX_CONSTANT_ARGS, &val);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Constant args failed " << err << '\n'*/;
	    else {
	      add_output("Max constant args = ");
	      add_output(val);
	      send_output();
	    }
	    err = devices[i].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &lval);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Constant args failed " << err << '\n'*/;
	    else {
	      if ((sl_int)lval < max_mem)
		max_mem = (sl_int)lval;
	      add_output("Max memory size = ");
	      add_output(lval);
	      send_output();
	    }
	    err = devices[i].getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &lval);
	    if (err != CL_SUCCESS)
	      /*std::cerr << "Constant args failed " << err << '\n'*/;
	    else {
	      if ((sl_int)lval < max_alloc)
		max_alloc = (sl_int)lval;
	      add_output("Max block allocation size = ");
	      add_output(lval);
	      send_output();
	    }

	    queue[k] = new cl::CommandQueue(context, devices[i],
				       CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE,
					    &err);
	    if (err == CL_SUCCESS)
	      k++;
	    else
	      report_error("CommandQueue::CommandQueue() failed -", err);
	  }
	  if (k == 0) {
	    ok = false;
	    num_devices = 0;
	    report_error("No queues for available devices");
	  } else {
	    if ((unsigned long)max_alloc > 0x7fffffff * sizeof(int))
	      report_error("Possible issue with buffer size and int offsets");
	  }
	}

	if (ok) {
	  //std::cout<<"Loading and compiling CL source\n";
	  const char name[] = "kernels.cl";
	  if (access(name, R_OK) != 0) {
	    report_error("We couldn't load CL source code");
	    ok = false;
	  } else {
	    cl::Program::Sources sources(1);
	    char *file_data = get_file(name);
	    sources[0] = std::make_pair(file_data, strlen(file_data));

	    program = cl::Program(context, sources, &err);
	    if (err != CL_SUCCESS) {
	      report_error("Program::Program() failed (", err);
	      ok = false;
	    } else {
	      err = program.build(devices);
	      if (err != CL_SUCCESS) {
		if (err == CL_BUILD_PROGRAM_FAILURE) {
		  std::string str = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
		  report_error("\nBUILD LOG");
		  report_error(str);
		} else
		  report_error("Program::build() failed (", err);
		ok = false;
	      }
	    }
	    // Todo - delete file data?
	  }
	}
      }
    }
  }
  if (!ok)
    report_error("OpenCL accelerator setup failure");
  accelerator = ok;
}

void machine::close_accelerator()
{
  for (int i = 0; i < num_devices; i++)
    delete queue[i];
}

bool machine::has_accelerator()
{
  return accelerator;
}

int machine::number_of_accelerators()
{
  return num_devices;
}

sl_int machine::available_mem(const int device)
{
  // force all to be the same size for simplicity in the rest of the code
  return max_mem;
}

sl_int machine::largest_alloc(const int device)
{
  // force all to be the same size for simplicity in the rest of the code
  return max_alloc;
}

dev_ptr machine::device_allocate(const sl_int size, const bool read_only,
				 const int device)
{
  cl_int status;
  cl_mem_flags flags = CL_MEM_READ_WRITE;
  if (read_only)
    flags = CL_MEM_READ_ONLY;
  cl::Buffer data(context, flags, size, 0, &status);
  if (status == CL_SUCCESS)
    return data;
  else {
    report_error("Accelerator allocation failed");
    //return 0;
    return data;
  }
}

// Todo - pass buffers by reference
void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
			     const int device)
{
  cl_int err;
  err = queue[device]->enqueueWriteBuffer(buffer, CL_TRUE, 0, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy to accelerator failed ", err);
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
			     const sl_int size, const int device)
{
  cl_int err;
  err = queue[device]->enqueueWriteBuffer(buffer, CL_TRUE, offset, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy to accelerator failed ", err);
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int size,
			       const int device)
{
  cl_int err;
  err = queue[device]->enqueueReadBuffer(buffer, CL_TRUE, 0, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy from accelerator failed ", err);
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			       const sl_int size, const int device)
{
  cl_int err;
  err = queue[device]->enqueueReadBuffer(buffer, CL_TRUE, offset, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy from accelerator failed ", err);
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
			     const int device, event_t *event)
{
  cl_int err;
  err = queue[device]->enqueueWriteBuffer(buffer, CL_FALSE, 0, size, ptr,
					  NULL, event);
  if (err != CL_SUCCESS)
    report_error("Copy to accelerator failed ", err);
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
			     const sl_int size, const int device,
			     event_t *event)
{
  cl_int err;
  err = queue[device]->enqueueWriteBuffer(buffer, CL_FALSE, offset, size, ptr,
					  NULL, event);
  if (err != CL_SUCCESS)
    report_error("Copy to accelerator failed ", err);
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int size,
			       const int device, event_t *event)
{
  cl_int err;
  err = queue[device]->enqueueReadBuffer(buffer, CL_FALSE, 0, size, ptr,
					 NULL, event);
  if (err != CL_SUCCESS)
    report_error("Copy from accelerator failed ", err);
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			       const sl_int size, const int device,
			       event_t *event)
{
  cl_int err;
  err = queue[device]->enqueueReadBuffer(buffer, CL_FALSE, offset, size, ptr,
					 0, event);
  if (err != CL_SUCCESS)
    report_error("Copy from accelerator failed ", err);
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int y_host,
			     const sl_int z_host, const sl_int x_size,
			     const sl_int y_size, const int z_size,
			     const int device, event_t *event)
{
  cl_int err;
  cl::size_t<3> dev_origin;
  dev_origin[0] = 0;
  dev_origin[1] = 0;
  dev_origin[2] = 0;
  cl::size_t<3> host_origin;
  host_origin[0] = 0;
  host_origin[1] = 0;
  host_origin[2] = 0;
  cl::size_t<3> region;
  region[0] = z_size;
  region[1] = y_size;
  region[2] = x_size;
  size_t host_1D = z_host;
  size_t host_2D = y_host * z_host;
  size_t dev_1D = z_size;
  size_t dev_2D = y_size * z_size;
  err = queue[device]->enqueueWriteBufferRect(buffer, CL_FALSE, dev_origin,
					      host_origin, region,  dev_1D,
					      dev_2D, host_1D, host_2D, ptr,
					      NULL, event);
  if (err != CL_SUCCESS)
    report_error("Copy rect to accelerator failed ", err);
}

void machine::copy_from_device(void *ptr, dev_ptr buffer, const sl_int y_host,
			     const sl_int z_host, const sl_int x_size,
			     const sl_int y_size, const int z_size,
			     const int device)
{
  cl_int err;
  cl::size_t<3> dev_origin;
  dev_origin[0] = 0;
  dev_origin[1] = 0;
  dev_origin[2] = 0;
  cl::size_t<3> host_origin;
  host_origin[0] = 0;
  host_origin[1] = 0;
  host_origin[2] = 0;
  cl::size_t<3> region;
  region[0] = z_size;
  region[1] = y_size;
  region[2] = x_size;
  size_t host_1D = z_host;
  size_t host_2D = y_host * z_host;
  size_t dev_1D = z_size;
  size_t dev_2D = y_size * z_size;
  err = queue[device]->enqueueReadBufferRect(buffer, CL_TRUE, dev_origin,
					     host_origin, region,  dev_1D,
					     dev_2D, host_1D, host_2D, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy rect from accelerator failed ", err);
}

void machine::device_free(dev_ptr data, const int device)
{
  // destructor should tidy
  //cl::Buffer *buf = (cl::Buffer *)data;
  //delete buf;
}

bool machine::check_accelerator_kernel(const char name[], const int device)
{
  cl_int err;
  cl::Kernel kernel(program, name, &err);
  return (err == CL_SUCCESS);
}

void machine::run_parallel_ah(const char name[], dev_ptr pix_buf,
			      dev_ptr vox_buf, const int vox_offset,
			      dev_ptr xy_buff, dev_ptr xy_offsets,
			      dev_ptr nlengths, const int nv,
			      const int nz, const int xy_size, const int ix,
			      const int size, const int dim3, const int device,
			      std::vector<event_t> *events)
{
  if (size > max_work01 or dim3 > max_work2)
    report_error("Too much work for device");
  else {
    cl_int err;
    cl::Kernel kernel(program, name, &err);
    if (err != CL_SUCCESS)
      report_error("Failed to setup OpenCL kernel");
    else {
      // It should be a multiple of 32 in the allocator, so group
      // based on that and use the 3rd dim if we want to submit
      // multiple jobs? not sure its ideal usage of the available space
      cl::NDRange globalThreads(size, dim3);
      cl_int status = CL_SUCCESS;
      status = kernel.setArg(0, pix_buf);
      if (status == CL_SUCCESS)
	status = kernel.setArg(1, vox_buf);
      if (status == CL_SUCCESS)
	status = kernel.setArg(2, vox_offset);
      if (status == CL_SUCCESS)
	status = kernel.setArg(3, xy_buff);
      if (status == CL_SUCCESS)
	status = kernel.setArg(4, xy_offsets);
      if (status == CL_SUCCESS)
	status = kernel.setArg(5, nlengths);
      if (status == CL_SUCCESS)
	status = kernel.setArg(6, nv);
      if (status == CL_SUCCESS)
	status = kernel.setArg(7, nz);
      if (status == CL_SUCCESS)
	status = kernel.setArg(8, xy_size);
      if (status == CL_SUCCESS)
	status = kernel.setArg(9, ix);
      if (status == CL_SUCCESS) {
	if (queue[device]->enqueueNDRangeKernel(kernel, cl::NullRange,
						globalThreads, cl::NullRange,
						events) != CL_SUCCESS) {
	  report_error("Failed to queue job");
	}
      }
      if (status != CL_SUCCESS)
	report_error("Failed to setup job args");
    }
  }
}

void machine::run_parallel_xy(const char name[], dev_ptr pix_buf,
			      dev_ptr vox_buf, dev_ptr xy_buff,
			      dev_ptr xy_offsets, dev_ptr h, dev_ptr lengths,
			      const int nv, const int nz, const int start,
			      const int ah_size, const int size,
			      const int dim3, const int device,
			      std::vector<event_t> *events)
{
  if (size > max_work01 or dim3 > max_work2)
    report_error("Too much work for device");
  else {
    cl_int err;
    cl::Kernel kernel(program, name, &err);
    if (err != CL_SUCCESS)
      report_error("Failed to setup OpenCL kernel");
    else {
      // It should be a multiple of 32 in the allocator, so group
      // based on that and use the 3rd dim if we want to submit
      // multiple jobs? not sure its ideal usage of the available space
      cl::NDRange globalThreads(size, dim3);
      cl_int status = CL_SUCCESS;
      status = kernel.setArg(0, pix_buf);
      if (status == CL_SUCCESS)
	status = kernel.setArg(1, vox_buf);
      if (status == CL_SUCCESS)
	status = kernel.setArg(2, xy_buff);
      if (status == CL_SUCCESS)
	status = kernel.setArg(3, xy_offsets);
      if (status == CL_SUCCESS)
	status = kernel.setArg(4, h);
      if (status == CL_SUCCESS)
	status = kernel.setArg(5, lengths);
      if (status == CL_SUCCESS)
	status = kernel.setArg(6, nv);
      if (status == CL_SUCCESS)
	status = kernel.setArg(7, nz);
      if (status == CL_SUCCESS)
	status = kernel.setArg(8, start);
      if (status == CL_SUCCESS)
	status = kernel.setArg(9, ah_size);
      if (status == CL_SUCCESS) {
	if (queue[device]->enqueueNDRangeKernel(kernel, cl::NullRange,
						globalThreads, cl::NullRange,
						events) != CL_SUCCESS) {
	  report_error("Failed to queue job");
	}
      }
      if (status != CL_SUCCESS)
	report_error("Failed to setup job args");
    }
  }
}

void machine::accelerator_barrier(const int device)
{
  queue[device]->enqueueBarrier();
}
 
#else

void machine::init_accelerator()
{
}

void machine::close_accelerator()
{
}

bool machine::has_accelerator()
{
  return false;
}

int machine::number_of_accelerators()
{
  return 0;
}

sl_int machine::available_mem(const int device)
{
  return -1;
}

sl_int machine::largest_alloc(const int device)
{
  return -1;
}

dev_ptr machine::device_allocate(const sl_int size, const bool read_only,
				 const int device)
{
  return 0;
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
			     const int device)
{
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
			     const sl_int size, const int device)
{
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			       const sl_int size, const int device)
{
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
			     const int device, event_t *event)
{
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
			     const sl_int size, const int device,
			     event_t *event)
{
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int size,
			       const int device, event_t *event)
{
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			       const sl_int size, const int device,
			       event_t *event)
{
}

void machine::device_free(dev_ptr data, const int device)
{
}

bool machine::check_accelerator_kernel(const char name[], const int device)
{
  return false;
}

void machine::run_parallel_ah(const char name[], dev_ptr pix_buf,
			      dev_ptr vox_buf, const int vox_offset,
			      dev_ptr xy_buff, dev_ptr xy_offsets,
			      dev_ptr nlengths, const int nv,
			      const int nz, const int xy_size, const int ix,
			      const int size, const int dim3, const int device,
			      std::vector<event_t> *events)
{
}

void machine::run_parallel_xy(const char name[], dev_ptr pix_buf,
			      dev_ptr vox_buf, dev_ptr xy_buff,
			      dev_ptr xy_offsets, dev_ptr h, dev_ptr lengths,
			      const int nv, const int nz, const int start,
			      const int ah_size, const int size,
			      const int dim3, const int device,
			      std::vector<event_t> *events)
{
}

void machine::accelerator_barrier(const int device)
{
}

#endif // USE_OPENCL
