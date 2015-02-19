
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
  cl::Context context;
  cl::Program program;
  cl::Kernel kernel;
  std::vector<cl::CommandQueue *> queue;

  sl_int max_alloc = 1000000000000;

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
	      if (int(wval[0] * wval[1]) < max_work01)
		max_work01 = wval[0] * wval[1];
	    }
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

	    queue[k] = new cl::CommandQueue(context, devices[i], 0, &err);
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
  cl::Buffer *data = new cl::Buffer(context, flags, size, 0, &status);
  if (status == CL_SUCCESS)
    return (void *)data;
  else {
    report_error("Accelerator allocation failed");
    delete data;
    return 0;
  }
}

// Todo - switch to CL_FALSE non-blocking transfers + events.
void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int size,
			     const int device)
{
  cl_int err;
  cl::Buffer *buf = (cl::Buffer *)buffer;
  err = queue[device]->enqueueWriteBuffer(*buf, CL_TRUE, 0, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy to accelerator failed");
}

void machine::copy_to_device(void *ptr, dev_ptr buffer, const sl_int offset,
			     const sl_int size, const int device)
{
  cl_int err;
  cl::Buffer *buf = (cl::Buffer *)buffer;
  err = queue[device]->enqueueWriteBuffer(*buf, CL_TRUE, offset, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy to accelerator failed");
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int size,
			       const int device)
{
  cl_int err;
  cl::Buffer *buf = (cl::Buffer *)buffer;
  err = queue[device]->enqueueReadBuffer(*buf, CL_TRUE, 0, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy from accelerator failed");
}

void machine::copy_from_device(dev_ptr buffer, void *ptr, const sl_int offset,
			       const sl_int size, const int device)
{
  cl_int err;
  cl::Buffer *buf = (cl::Buffer *)buffer;
  err = queue[device]->enqueueReadBuffer(*buf, CL_TRUE, offset, size, ptr);
  if (err != CL_SUCCESS)
    report_error("Copy from accelerator failed");
}

void machine::device_free(dev_ptr data, const int device)
{
  // destructor should tidy
  cl::Buffer *buf = (cl::Buffer *)data;
  delete buf;
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
			      const int n, const int nv,
			      const int nz, const int size, const int device)
{
  if (n > max_work01)
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
      cl::NDRange globalThreads(32, size / 32, 1);
      cl_int status = CL_SUCCESS;
      status = kernel.setArg(0, *((cl::Buffer *)pix_buf));
      if (status == CL_SUCCESS)
	status = kernel.setArg(1, *((cl::Buffer *)vox_buf));
      if (status == CL_SUCCESS)
	status = kernel.setArg(2, vox_offset);
      if (status == CL_SUCCESS)
	status = kernel.setArg(3, *((cl::Buffer *)xy_buff));
      if (status == CL_SUCCESS)
	status = kernel.setArg(4, *((cl::Buffer *)xy_offsets));
      if (status == CL_SUCCESS)
	status = kernel.setArg(5, n);
      if (status == CL_SUCCESS)
	status = kernel.setArg(6, nv);
      if (status == CL_SUCCESS)
	status = kernel.setArg(7, nz);
      if (status == CL_SUCCESS) {
	if (queue[device]->enqueueNDRangeKernel(kernel, cl::NullRange,
						globalThreads, 0)
	    != CL_SUCCESS) {
	  report_error("Failed to queue job");
	}
	// Todo get event and wait on it to make synchronous?
      }
      if (status != CL_SUCCESS)
	report_error("Failed to setup job args");
    }
  }
}

void machine::run_parallel_xy(const char name[], dev_ptr pix_buf,
			      const int offset, dev_ptr vox_buf,
			      dev_ptr xy_buff, dev_ptr xy_offsets, const int n,
			      const int nv, const int nz, const int size,
			      const int device)
{
  if (n > max_work01)
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
      cl::NDRange globalThreads(32, size / 32, 1);
      cl_int status = CL_SUCCESS;
      status = kernel.setArg(0, *((cl::Buffer *)pix_buf));
      if (status == CL_SUCCESS)
	status = kernel.setArg(1, offset);
      if (status == CL_SUCCESS)
	status = kernel.setArg(2, *((cl::Buffer *)vox_buf));
      if (status == CL_SUCCESS)
	status = kernel.setArg(3, *((cl::Buffer *)xy_buff));
      if (status == CL_SUCCESS)
	status = kernel.setArg(4, *((cl::Buffer *)xy_offsets));
      if (status == CL_SUCCESS)
	status = kernel.setArg(5, n);
      if (status == CL_SUCCESS)
	status = kernel.setArg(6, nv);
      if (status == CL_SUCCESS)
	status = kernel.setArg(7, nz);
      if (status == CL_SUCCESS) {
	if (queue[device]->enqueueNDRangeKernel(kernel, cl::NullRange,
						globalThreads, 0)
	    != CL_SUCCESS) {
	  report_error("Failed to queue job");
	}
      }
      if (status != CL_SUCCESS)
	report_error("Failed to setup job args");
    }
  }
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
			      const int n, const int nv,
			      const int nz, const int size, const int device)
{
}

void machine::run_parallel_xy(const char name[], dev_ptr pix_buf,
			      const int offset, dev_ptr vox_buf,
			      dev_ptr xy_buff, dev_ptr xy_offsets, const int n,
			      const int nv, const int nz, const int size,
			      const int device)
{
}

#endif // USE_OPENCL