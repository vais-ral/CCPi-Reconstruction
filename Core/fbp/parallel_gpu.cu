#include "defs_gpu.cuh"

#include "time_stamp.h"

texture<float, 1, cudaReadModeElementType> texC;
texture<float, 1, cudaReadModeElementType> texA;
texture<float, 1, cudaReadModeElementType> texX;
texture<float, 1, cudaReadModeElementType> texY;

texture<float, 1, cudaReadModeElementType> texCos;
texture<float, 2, cudaReadModeElementType> texAR;

__global__ void find_cart_new(float * pc_res_d, float *pc_a_d, float *pc_r_d, int cart_len){
	int start = blockIdx.x*cart_len;
	for(int i=start;i<(start+cart_len);i++){
		//pc_res_d[i]= tex2D(texPC,pc_a_d[i],pc_r_d[i]);
		pc_res_d[i]= tex2D(texAR,pc_a_d[i]+0.5f,pc_r_d[i]+0.5f);
	}
}



static __device__ __host__ inline Complex ComplexMul(Complex a, Complex b)
{
    Complex c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}


__global__ void cuda_data_r2c(Complex *dev_inc, float *dev_in, unsigned int sp){
	unsigned int i, j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = i + sp;
	dev_inc[i].x = dev_in[j];
	dev_inc[i].y = 0.0f;
}


__global__ void cuda_data_c2r_3(float *dev_in, Complex *dev_inc, unsigned int shift, unsigned int nc, unsigned int nr, unsigned int nrow){
	unsigned int i, i2, j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	i2 = i + shift;
	j =  blockIdx.y * blockDim.y + threadIdx.y;
	dev_in[(nrow+j)*nr+i] = dev_inc[j*nc+i2].x;
}



__global__ void cuda_data_c2r(float *dev_in, Complex *dev_inc, unsigned int sp){
	unsigned int i, j;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = i + sp;
	dev_in[j] = dev_inc[i].x;
}

__global__ void cuda_mul_c(const Complex *dev_fc, Complex *dev_inc){
	unsigned int i;
	Complex c;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	c = ComplexMul(dev_fc[i],dev_inc[i]);
	dev_inc[i] = c;
}


__global__ void cprod(Complex* a, const Complex* b)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	a[i] = ComplexMul(a[i], b[i]);
}

__global__ void creal(const Complex* a, float* b)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
	b[i] = a[i].x;
} 


__global__ void cuda_filter_r2c(Complex *dev_fc, float *dev_fr, unsigned int nx){
	unsigned int i, j, k;
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = i%nx;
	k = blockIdx.y * nx + i;
	dev_fc[k].x = dev_fr[j];
	dev_fc[k].y = 0.0f;
}



__global__ void fbp_axial(float* dev_out, float da, unsigned int nx, unsigned int ny, float fnr){
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


__global__ void set_zero(float *veco, unsigned int nxo){	
	unsigned int tig;
	tig = (blockIdx.y*blockDim.y + threadIdx.y)*nxo + blockIdx.x*blockDim.x + threadIdx.x;
	veco[tig] = 0.0f;
}

bool get_gpu_info(allData *aD){
	

/*	
	aD->gi->major = 1;
	aD->gi->minor = 1;
	aD->gi->multiProcessorCount = 16;
	aD->gi->regsPerBlock = 8*1024;
	aD->gi->warpSize = 32;
	aD->gi->sharedMemPerBlock = 16*1024;
//	aD->gi->maxResidentThreads = prop.maxThreadsPerMultiProcessor;//
	aD->gi->maxThreadsPerBlock = 512;

	aD->gi->maxResidentBlocks = 8;
	if(aD->gi->major == 2){
		aD->gi->maxResidentThreads = 1536;//new
		aD->gi->sharedMemBanks = 32;
		aD->gi->maxResidentWarps = 48;
	}else if(aD->gi->minor > 1){
		aD->gi->maxResidentThreads = 1024;//new
		aD->gi->sharedMemBanks = 16;
		aD->gi->maxResidentWarps = 32;
	}else{
		aD->gi->maxResidentThreads = 768;//new
		aD->gi->sharedMemBanks = 16;
		aD->gi->maxResidentWarps = 24;
	}
*/
	/*aD->gi->major = 2;
	aD->gi->minor = 0;
	aD->gi->multiProcessorCount = 30;
	aD->gi->regsPerBlock = 32*1024;
	aD->gi->warpSize = 32;
	aD->gi->sharedMemPerBlock = 48*1024;
//	aD->gi->maxResidentThreads = prop.maxThreadsPerMultiProcessor;//
	aD->gi->maxThreadsPerBlock = 1024;

	aD->gi->maxResidentBlocks = 8;
	if(aD->gi->major == 2){
		aD->gi->maxResidentThreads = 1536;//new
		aD->gi->sharedMemBanks = 32;
		aD->gi->maxResidentWarps = 48;
	}else if(aD->gi->minor > 1){
		aD->gi->maxResidentThreads = 1024;//new
		aD->gi->sharedMemBanks = 16;
		aD->gi->maxResidentWarps = 32;
	}else{
		aD->gi->maxResidentThreads = 768;//new
		aD->gi->sharedMemBanks = 16;
		aD->gi->maxResidentWarps = 24;
	}*/
	
	cudaError_t ce;
	int deviceCount;
        int MyDevice;

	printTagStart(aD,"GPUDevices");
	
	ce = cudaGetDeviceCount(&deviceCount);
	if(ce != cudaSuccess){
		sprintf(aD->message, "can not count GPUs: \"%s\"", cudaGetErrorString(ce));
		printError(aD);
		return false;
	}else{
		printTag(aD,"NumberOfDevices",deviceCount);
	}

	if(deviceCount < 1 || deviceCount > 10){
		printError(aD,"device count - wrong number");
		return false;
	}
	
	printTag(aD,"RequestedIndexOfDevice",aD->hms->fbp->GPUDeviceNumber);
        cudaGetDevice(&MyDevice);
	printTag(aD,"MyCudaDevice",MyDevice);
	
	cudaDeviceProp prop;



	cudaGetDeviceProperties(&prop,MyDevice);


	aD->gi->major = (unsigned int)prop.major;
	aD->gi->minor = (unsigned int)prop.minor;
	aD->gi->multiProcessorCount = (unsigned int)prop.multiProcessorCount;
	aD->gi->regsPerBlock = (unsigned int)prop.regsPerBlock;
	aD->gi->warpSize = (unsigned int)prop.warpSize;
	aD->gi->sharedMemPerBlock = (unsigned int)(prop.sharedMemPerBlock);
//	aD->gi->maxResidentThreads = prop.maxThreadsPerMultiProcessor;//
	aD->gi->maxThreadsPerBlock = (unsigned int)prop.maxThreadsPerBlock;

	aD->gi->maxResidentBlocks = 8;
	if(aD->gi->major == 2){
		aD->gi->maxResidentThreads = 1536;//new
		aD->gi->sharedMemBanks = 32;
		aD->gi->maxResidentWarps = 48;
	}else if(aD->gi->minor > 1){
		aD->gi->maxResidentThreads = 1024;//new
		aD->gi->sharedMemBanks = 16;
		aD->gi->maxResidentWarps = 32;
	}else{
		aD->gi->maxResidentThreads = 768;//new
		aD->gi->sharedMemBanks = 16;
		aD->gi->maxResidentWarps = 24;
	}

	printTagStart(aD,"GPUDevice");
	printTag(aD,"Name",prop.name);
	printTag(aD,"MultiProcessorCount",(unsigned int)prop.multiProcessorCount);
	printTag(aD,"ClockRate",prop.clockRate/1000,"in MHz");
	printTag(aD,"TotalGlobalMemory",(unsigned int)(prop.totalGlobalMem),"bytes");
	printTag(aD,"TotalGlobalMemory",(float)(float(prop.totalGlobalMem)/1073741824.0),"GB");
	printTag(aD,"RevisionMajor",(unsigned int)prop.major);
	printTag(aD,"RevisionMinor",(unsigned int)prop.minor);

	sprintf(aD->message, "%i x %i x %i", (unsigned int)prop.maxGridSize[0], (unsigned int)prop.maxGridSize[1], (unsigned int)prop.maxGridSize[2]);
	printTag(aD,"MaximumGridSize",aD->message,"maximum size of a grid of thread blocks");

	sprintf(aD->message, "%i x %i x %i", (unsigned int)prop.maxThreadsDim[0], (unsigned int)prop.maxThreadsDim[1], (unsigned int)prop.maxThreadsDim[2]);
	printTag(aD,"MaximumThreadSize",aD->message,"the maximum size of each dimension of a block");

	printTag(aD,"WarpSize",(unsigned int)(prop.warpSize));
	printTag(aD,"MaxWarpsPerMultiProcessor",(unsigned int)(aD->gi->maxResidentWarps),"maximum number of resident warps per multiprocessor");
	
	printTag(aD,"MaxThreadsPerBlock",(unsigned int)(prop.maxThreadsPerBlock));
//	printTag(aD,"MaxThreadsPerMultiProcessor",prop.maxThreadsPerMultiProcessor);

	printTag(aD,"MaxBlocksPerMultiProcessor",(unsigned int)(aD->gi->maxResidentBlocks),"maximum number of resident blocks per multiprocessor");

	printTag(aD,"RegistersPerBlock",(unsigned int)(prop.regsPerBlock),"maximum number of 32-bit registers available to a thread block");
	printTag(aD,"SharedMemoryPerBlock",(unsigned int)(prop.sharedMemPerBlock),"maximum amount of shared memory available to a thread block in bytes");
	printTag(aD,"SharedMemoryBanks",(unsigned int)(aD->gi->sharedMemBanks),"number of shared memory banks");

//	printTag(aD,"L2cache",prop.l2CacheSize,"bytes");

	printTagEnd(aD);//GPUDevice
	printTagEnd(aD);//GPUDevices
	
	return true;
}




bool fbp_axial_cuda(allData *aD){
	gpu_info *gi;
	xData *xd;
	
	unsigned int mb, mt, nx, ny;
	unsigned int uu, ub, maxt;
	float fnr, da;
	size_t mem_size, mem_size_c;
		
	float *dev_in, *dev_out;
	cufftComplex *dev_filter, *dev_data;

	cufftHandle plan;
        timestamp("Starting fbp_axial_cuda",4);

	xd = aD->data;
	gi = aD->gi;

	nx = aD->ta_nx;
	ny = aD->ny;

	mb = gi->maxResidentBlocks;
	mt = gi->maxResidentThreads;
	
	maxt = mt/mb;
	
	uu = 1;
	while(maxt>1){
		uu*=2;
		maxt/=2;	
	}
	maxt = uu;
	ub = nx/maxt;

	da = aD->gi->rotAngleStep;

	fnr = 0.5f*float(nx-1);

	dim3 grids_p(ub,1,1);
	dim3 threads_p(maxt,1,1);
          

	mem_size = nx*sizeof(float);
	mem_size_c = nx*sizeof(cufftComplex);
	
        timestamp("calling cuda to automatically get a device",4);
        cudaFree(NULL);
	// cudaSetDevice(aD->hms->fbp->GPUDeviceNumber);
         {
		int MyDevice;
		char  devmsg[128];
                cudaGetDevice(&MyDevice);
		snprintf(devmsg,128,"my device %i",MyDevice);
		timestamp(devmsg,4);
          }

        timestamp("calling cudaMalloc",4);
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_filter,mem_size_c));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_data,mem_size_c));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_in,mem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_out,mem_size));

        timestamp("finished cudaMalloc",4);
	
	CUDA_SAFE_CALL(cudaMemcpy(dev_filter,xd->veccF,mem_size_c,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dev_data,xd->veccI,mem_size_c,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaThreadSynchronize());

	HMFFT_SAFE_CALL(cufftPlan1d(&plan, nx, CUFFT_C2C, 1));

	HMFFT_SAFE_CALL(cufftExecC2C(plan, dev_data, dev_data, CUFFT_FORWARD));
	CUDA_SAFE_CALL(cudaThreadSynchronize());

	cprod<<<grids_p,threads_p>>>(dev_data,dev_filter);
	CUDA_SAFE_CALL(cudaThreadSynchronize());

	HMFFT_SAFE_CALL(cufftExecC2C(plan, dev_data, dev_data, CUFFT_INVERSE));
	CUDA_SAFE_CALL(cudaThreadSynchronize());

	creal<<<grids_p,threads_p>>>(dev_data,dev_in);
	CUDA_SAFE_CALL(cudaThreadSynchronize());

	cufftDestroy(plan);

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	texA.normalized = false;
   	texA.filterMode = cudaFilterModeLinear;
   	texA.addressMode[0] = cudaAddressModeClamp;

	CUDA_SAFE_CALL(cudaBindTexture(NULL, texA, dev_in, channelDesc, mem_size));
	
	fbp_axial<<<grids_p, threads_p>>>(dev_out, da, nx, ny, fnr);
	CUDA_SAFE_CALL(cudaThreadSynchronize());

	CUDA_SAFE_CALL(cudaMemcpy(xd->vecto, dev_out, mem_size, cudaMemcpyDeviceToHost));

	cudaFree(dev_in); dev_in = NULL;
	cudaFree(dev_out); dev_out = NULL;
	cudaFree(dev_data); dev_data = NULL;
	cudaFree(dev_filter); dev_filter = NULL;
        timestamp("Finishing fbp_axial_cuda",4);
				
	return true;
}





__global__ void fbp_std2(float* temp_dev, unsigned int nxo, unsigned int mw, float xc, float *dev_CS, unsigned int vH, unsigned int chnum){
	extern __shared__ float shm[];
	float2 *vcs = (float2 *)&shm[0];
	float sum, p, fitt;
	float x, y;
	unsigned int ixx, iyy, itt, j, tib, tpb;

	tib = threadIdx.y * blockDim.x + threadIdx.x;
	tpb = blockDim.x * blockDim.y;
	ixx = blockIdx.x*blockDim.x + threadIdx.x;
	iyy = blockIdx.y*blockDim.y + threadIdx.y;
	itt = ixx+iyy*nxo;
	fitt = float(itt);
	
	x = tex1Dfetch(texX,fitt);
	y = tex1Dfetch(texY,fitt);
		
	sum = 0.0f;
	
	shm[tib] = dev_CS[tib+chnum*tpb];
	__syncthreads();
	for(j=0; j<vH; j++){
		p = xc + x*vcs[j].x+ y*vcs[j].y + float(j*mw);
		sum += (tex1Dfetch(texA,p));	
	}
	
	temp_dev[itt] += sum;
}


bool fbp_cuda(allData *aD){
	xData *xd = aD->data;

	size_t mem_CS;
	size_t mem_out, mem_in;
	size_t mem_chunk;
	size_t mem_comp;
	size_t mem_shared;

	float *dev_in, *dev_out;
	float *dev_X, *dev_Y;
	float *dev_CS;
		
	unsigned int hTimer;
	unsigned int nx, ny, nxo, nyo, wo, ho;
	unsigned int cl, cw;
	unsigned int vH, sp;
		
	cufftComplex *dev_fc, *dev_inc;
	cufftHandle plan;
	double gpuTime;
	
	xFBP_gpu *fg;
        timestamp("starting fbp_cuda",4);
	fg = aD->fg;

	nx = fg->nx;
	ny = fg->ny;
	nxo = fg->nxo;
	nyo = fg->nyo;
	wo = fg->blockWidth;
	ho = fg->blockHeight;
	cl = fg->chunkLeft;
	cw = fg->chunkWidth;
	vH = fg->vH;
	
	mem_shared = wo*ho*sizeof(float);		
	mem_comp = nx*vH*sizeof(cufftComplex);
	mem_chunk = cw*vH*sizeof(float);
	mem_CS = 2*ny*sizeof(float);
	mem_out = nxo*nyo*sizeof(float);
	mem_in = nx*ny*sizeof(float);

	dim3 grids_in(nx*vH/(wo*ho),1);
	dim3 threads_in(wo*ho,1);

	dim3 grids_ch(cw/wo,vH/ho,1);
	dim3 threads_ch(wo,ho,1);
	
	dim3 grids_bp(nxo/wo,nyo/ho,1);	
	dim3 threads_bp(wo,ho,1);
	
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	
	CUT_SAFE_CALL( cutCreateTimer(&hTimer) );
	CUT_SAFE_CALL( cutResetTimer(hTimer));
    CUT_SAFE_CALL( cutStartTimer(hTimer));

        timestamp("calling cudaMalloc",4);
         {
		int MyDevice;
		char  devmsg[128];
                cudaGetDevice(&MyDevice);
		snprintf(devmsg,128,"my device %i",MyDevice);
		timestamp(devmsg,4);
          }
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_fc,mem_comp));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_inc,mem_comp));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_in,mem_in));
        timestamp("finished cudaMalloc",4);

	CUDA_SAFE_CALL(cudaMemcpy(dev_fc,xd->veccF,mem_comp,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dev_in,xd->vta,mem_in,cudaMemcpyHostToDevice));

	CUFFT_SAFE_CALL(cufftPlan1d(&plan, nx, CUFFT_C2C, vH));

	texA.normalized = false;
   	texA.filterMode = cudaFilterModeLinear;
   	texA.addressMode[0] = cudaAddressModeClamp;

	texX.normalized = false;
   	texX.filterMode = cudaFilterModeLinear;
   	texX.addressMode[0] = cudaAddressModeClamp;

	texY.normalized = false;
   	texY.filterMode = cudaFilterModeLinear;
   	texY.addressMode[0] = cudaAddressModeClamp;


        timestamp("calling cudaMalloc",4);
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_CS,mem_CS));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_out,mem_out));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_X,mem_out));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_Y,mem_out));
        timestamp("finished cudaMalloc",4);

	CUDA_SAFE_CALL(cudaMemcpy(dev_X,xd->vecX,mem_out,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dev_Y,xd->vecY,mem_out,cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dev_CS,xd->vecCS,mem_CS,cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaBindTexture(NULL, texA, dev_in, channelDesc, mem_chunk));	
	CUDA_SAFE_CALL(cudaBindTexture(NULL, texX, dev_X, channelDesc, mem_out));
	CUDA_SAFE_CALL(cudaBindTexture(NULL, texY, dev_Y, channelDesc, mem_out));



	set_zero<<<grids_bp, threads_bp>>>(dev_out, nxo);
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	for(unsigned int i=0;i<fg->numChunks; i++){
		sp = i * nx * vH;
		cuda_data_r2c<<<grids_in,threads_in>>>(dev_inc,dev_in,sp);
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
		CUFFT_SAFE_CALL(cufftExecC2C(plan, dev_inc, dev_inc, CUFFT_FORWARD));
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
		
		cuda_mul_c<<<grids_in,threads_in>>>(dev_fc,dev_inc);
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
		CUFFT_SAFE_CALL(cufftExecC2C(plan, dev_inc, dev_inc, CUFFT_INVERSE));
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
		cuda_data_c2r_3<<<grids_ch,threads_ch>>>(dev_in, dev_inc, cl, nx, cw, 0);
		
		CUDA_SAFE_CALL( cudaThreadSynchronize() );	
		
		fbp_std2<<<grids_bp, threads_bp, mem_shared>>>(dev_out, nxo, cw, fg->xc, dev_CS, vH, i);
		CUDA_SAFE_CALL( cudaThreadSynchronize() );
	}

	cufftDestroy(plan);
	cudaFree(dev_inc); dev_inc = NULL;
	cudaFree(dev_fc); dev_fc = NULL;

	CUDA_SAFE_CALL(cudaMemcpy(xd->vto,dev_out,mem_out,cudaMemcpyDeviceToHost));
	
	cudaFree(dev_in); dev_in = NULL;
	cudaFree(dev_out); dev_out = NULL;
	cudaFree(dev_CS); dev_CS = NULL;
	cudaFree(dev_X); dev_X = NULL;
	cudaFree(dev_Y); dev_Y = NULL;
	
	CUT_SAFE_CALL( cutStopTimer(hTimer) );
    gpuTime = cutGetTimerValue(hTimer);
	printTag(aD,"BackprojectionTime",(float)gpuTime,"in ms");
        timestamp("finished fbp_cuda",4);
    
	return true;
}


__global__ void xb_proj_new(unsigned int indexs, unsigned int starty, float *vif, float *vecxy, float *vec_sc, float xc, unsigned int vH){
	unsigned int ij, ijk, tib, i;
	float s;
	float2 cs, xy;
	
	tib = (blockDim.x*threadIdx.y + threadIdx.x)/2;
	cs = *(float2 *)(&vec_sc[2*(tib+starty)]);
		
	ij = 2*blockIdx.x + indexs;
	for(i = 0; i<2; i++){	
		xy = *(float2 *)(&vecxy[2*ij]);	
		ijk = ij*vH+tib;
		s = __fmaf_rn(xy.x,cs.x,xc);
		vif[ijk] =  __fmaf_rn(xy.y,cs.y,s);
		ij++;
	}

}


__global__ void xb_proj2_new(unsigned int starty, float *vec_x, float *vecxy, float *vec_cs){
	unsigned int ij, ijk, tib, i, tpb;
	float s;
	float2 cs, xy;

	tpb = blockDim.x * blockDim.y;
	tib = (blockDim.x*threadIdx.y + threadIdx.x)/2;
	cs = *(float2 *)(&vec_cs[2*(tib+starty)]);
		
	ij = 2*blockIdx.x;
	for(i = 0; i<2; i++){
		xy = *(float2 *)(&vecxy[2*ij]);	
		ijk = tib*tpb+ij;
		s = __fmul_rn(xy.x,cs.x);
		vec_x[ijk] = __fmaf_rn(xy.y,cs.y,s);
		ij++;
	}
}


__global__ void fbp_std5(float *vb, float *veco, float *vt, unsigned int nx, unsigned int vH, unsigned int nxo, unsigned int nxb){
	float rr, sum;
	
	unsigned int ixq, ivc, ic;
	unsigned int tib, tig, big, tpb;
		
	tib = blockDim.x*threadIdx.y + threadIdx.x;
	big = blockIdx.y*nxb + blockIdx.x;
	tig = (blockIdx.y*blockDim.y + threadIdx.y)*nxo + blockIdx.x*blockDim.x + threadIdx.x;
	tpb = blockDim.x * blockDim.y;
	
	ixq = tib;
	ivc = big*vH;
	sum = 0.0f;	
	
	for(ic = 0; ic<vH; ic++){		
		rr = __fadd_rn(vt[ixq],vb[ivc+ic]);
		sum += tex1Dfetch(texA,rr+float(ic*nx));	
		ixq += tpb;
	}
	
	veco[tig] += sum;
}


bool fbp_cuda_20(allData *aD){
	xData *xd = aD->data;

	unsigned int nx, ny, nxo, nyo;
	unsigned int wo, ho;
	unsigned int nxb, nyb;
	
	unsigned int mt, mb;
	unsigned int uu, ub;
	unsigned int maxt;
	unsigned int vH;
	unsigned int nf, sp;
	unsigned int starty, indexs;
	unsigned int nxyb, nlen, nxyb2;

	unsigned int hTimer;

	size_t mem_fc, mem_chunk;
	size_t mem_in, mem_out;
	size_t mem_CS, mem_block_xy;
	size_t mem_thread_x, mem_thread_xy, mem_block;
	
	float xc;

	float *dev_thread_x, *dev_thread_xy, *dev_block;
	float *dev_in, *dev_out;
	float *dev_block_xy;
	float *dev_CS;
	
	gpu_info *gi;
	
	double gpuTime;
			
	cufftComplex *dev_fc, *dev_inc;

	cufftHandle plan;
        timestamp("starting fbp_cuda_20",4);

	gi = aD->gi;
	
	xc = aD->new_xc;
	vH = gi->vertH;

	wo = gi->wo;
	ho = gi->ho;
	
	nx = aD->ta_nx;
	ny = aD->ta_ny;
	nxo = gi->mxo;
	nyo = gi->myo;	
	
	nxb = nxo/wo;
	nyb = nyo/ho;

	mem_block = nxb*nyb*vH*sizeof(float); 
	mem_thread_x = wo*ho*vH*sizeof(float);
	mem_thread_xy = 2*wo*ho*sizeof(float);
	mem_block_xy = 2 * nxb * nyb * sizeof(float);

	mem_in = nx * ny *sizeof(float);
	mem_out = nxo * nyo *sizeof(float);
	mem_CS = 2*ny*sizeof(float);
	mem_fc = nx * vH * sizeof(cufftComplex);
	mem_chunk = nx * vH * sizeof(float);
	
	mb = gi->maxResidentBlocks;
	mt = gi->maxResidentThreads;
	
	maxt = mt/mb;
	
	uu = 1;
	while(maxt>1){
		uu*=2;
		maxt/=2;	
	}
	maxt = uu;
	ub = nx/maxt;

	nf = ny/vH;

	nxyb = nxb*nyb;
	
	nlen = (nxyb/4);
	if(nxyb%4 > 0) nlen++;
	
	nxyb2 = nxyb-2*nlen;

	dim3 threads_xb(1,1,1);
	dim3 threads_1(maxt,1,1);
	dim3 threads_bp(wo,ho,1);

	dim3 grid_xb2(vH,1,1);	
	dim3 grid_xb_new(nlen,1,1);	
	dim3 grids_2(ub*vH,1,1);
	dim3 grid_bp(nxb,nyb,1);	
	
	HMCUT_SAFE_CALL( cutCreateTimer(&hTimer) );
	HMCUT_SAFE_CALL( cutResetTimer(hTimer));
        HMCUT_SAFE_CALL( cutStartTimer(hTimer));
	
        timestamp("calling cudaMalloc ",4);
         {
		int MyDevice;
		char  devmsg[128];
                cudaGetDevice(&MyDevice);
		snprintf(devmsg,128,"my device %i",MyDevice);
		timestamp(devmsg,4);
          }
	HM_SAFE_CALL(cudaMalloc((void**)&dev_block,mem_block));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_thread_x,mem_thread_x));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_thread_xy,mem_thread_xy));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_fc,mem_fc));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_inc,mem_fc));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_in,mem_in));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_out,mem_out));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_block_xy,mem_block_xy));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_CS,mem_CS));
        timestamp("finished  cudaMalloc ",4);

	HM_SAFE_CALL(cudaMemcpy(dev_fc,xd->veccF,mem_fc,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(dev_in,xd->vta,mem_in,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(dev_block_xy,xd->vecbXY,mem_block_xy,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(dev_thread_xy,xd->vecXY_block,mem_thread_xy,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(dev_CS,xd->vecCS,mem_CS,cudaMemcpyHostToDevice));

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	texA.normalized = false;
   	texA.filterMode = cudaFilterModeLinear;
   	texA.addressMode[0] = cudaAddressModeClamp;

	CUDA_SAFE_CALL(cudaBindTexture(NULL, texA, dev_in, channelDesc, mem_chunk));	

        /*rca: AARGH -- plan is done for each slice -- this defeats the purpose of planning the FFT */
        timestamp("before plan creation",4);
	HMFFT_SAFE_CALL(cufftPlan1d(&plan, nx, CUFFT_C2C, vH));
        timestamp("after plan creation",4);

	set_zero<<<grid_bp, threads_bp>>>(dev_out, nxo);
	HM_SAFE_CALL( cudaThreadSynchronize() );
        timestamp("before backprojection loop ",4);
	
	for(unsigned int i = 0;i<nf;i++){
		sp = i*nx*vH;
		cuda_data_r2c<<<grids_2,threads_1>>>(dev_inc,dev_in,sp);
		HM_SAFE_CALL( cudaThreadSynchronize() );
		HMFFT_SAFE_CALL(cufftExecC2C(plan, dev_inc, dev_inc, CUFFT_FORWARD));
		HM_SAFE_CALL( cudaThreadSynchronize() );
		cuda_mul_c<<<grids_2,threads_1>>>(dev_fc,dev_inc);
		
		HM_SAFE_CALL( cudaThreadSynchronize() );
		HMFFT_SAFE_CALL(cufftExecC2C(plan, dev_inc, dev_inc, CUFFT_INVERSE));
		HM_SAFE_CALL( cudaThreadSynchronize() );
		sp = 0;
		cuda_data_c2r<<<grids_2,threads_1>>>(dev_in,dev_inc,sp);
		HM_SAFE_CALL( cudaThreadSynchronize() );

		starty = i*vH;
		indexs = 0;
		xb_proj_new<<<grid_xb_new, threads_bp>>>(indexs, starty, dev_block, dev_block_xy, dev_CS, xc, vH);
		HM_SAFE_CALL( cudaThreadSynchronize() );
		indexs = nxyb2;
		xb_proj_new<<<grid_xb_new, threads_bp>>>(indexs, starty, dev_block, dev_block_xy, dev_CS, xc, vH);
		HM_SAFE_CALL( cudaThreadSynchronize() );

		xb_proj2_new<<<grid_xb2, threads_bp>>>(starty, dev_thread_x, dev_thread_xy, dev_CS);
		HM_SAFE_CALL( cudaThreadSynchronize() );
		fbp_std5<<<grid_bp, threads_bp>>>(dev_block, dev_out, dev_thread_x, nx, vH, nxo, nxb);
		HM_SAFE_CALL( cudaThreadSynchronize() );
	}
        timestamp("after backprojection loop ",4);
	
	cufftDestroy(plan);
        timestamp("after plan destruction",4);
		
	HM_SAFE_CALL(cudaMemcpy(xd->vto, dev_out, mem_out, cudaMemcpyDeviceToHost));

	cudaFree(dev_fc); dev_fc = NULL;
	cudaFree(dev_inc); dev_inc = NULL;
	cudaFree(dev_in); dev_in = NULL;
	cudaFree(dev_block); dev_block = NULL;
	cudaFree(dev_thread_x); dev_thread_x = NULL;
	cudaFree(dev_thread_xy); dev_thread_xy = NULL;
	cudaFree(dev_out); dev_out = NULL;
	cudaFree(dev_CS); dev_CS = NULL;
	cudaFree(dev_block_xy); dev_block_xy = NULL;

	HMCUT_SAFE_CALL( cutStopTimer(hTimer) );

        gpuTime = cutGetTimerValue(hTimer);

	printTag(aD,"TimeBackprojection",float(gpuTime),"time (in ms)");
   	
        timestamp("finished fbp_cuda_20",4);
	return true;
}


__global__ void fbp_cpu2(float *dev_in, float *dev_pol, float *dev_Cos, unsigned int starty, unsigned int nx, unsigned int na, unsigned int vH, unsigned int pola, float ps, float xc_new){
	float x, xx, sum;
	
	unsigned int ic, irr;
	unsigned int tib, tig, tpb;
	int iaa, aa, ia, dr, dm;
		
	tib = blockDim.x*threadIdx.y + threadIdx.x;
	iaa = blockIdx.x*blockDim.x + threadIdx.x;
	irr = blockIdx.y*blockDim.y + threadIdx.y;
	tig = irr*pola + iaa;
	tpb = blockDim.x * blockDim.y;

	x = ps * float(irr);
		
	sum = 0.0f;	
	float t;

	ia = starty + iaa;
	dm = ia/na;
	dr = ia - dm * na;
	if(dm%2 == 0){
		t= 1.0f;
	}else{
		t= -1.0f;
	}
	
	for(ic = 0; ic<vH; ic++){
		if(dr == na){
			dr = 0;
			t*= -1.0f;
		}
		
		xx = xc_new + t* x * tex1Dfetch(texCos, dr);
		sum += tex1Dfetch(texA,xx + float(ic*nx));	
		dr++;
		
	}
	
	dev_pol[tig] += sum;
}



bool fbp_cuda_cpu2(allData *aD){
	xData *xd = aD->data;

	unsigned int nx, ny, polr, pola, nba, nbr;
	unsigned int wo, ho;
		
	unsigned int mt, mb;
	unsigned int uu, ub;
	unsigned int maxt;
	unsigned int vH;
	unsigned int nf, sp, na;
	unsigned int starty, indexs;
	
	unsigned int hTimer;

	size_t mem_fc, mem_chunk;
	size_t mem_in;
	size_t mem_Cos, mem_pol;
		
	float xc, ps;
	
	float *dev_in, *dev_pol;
	float *dev_Cos;
	
	gpu_info *gi;
	
	double gpuTime;
        timestamp("calling fbp_cuda_cpu2",4);

	gi = aD->gi;

	pola = gi->pol_a;
	polr = gi->pol_r;

	ps = gi->outputPixelSize;
			
	cufftComplex *dev_fc, *dev_inc;

	cufftHandle plan;

	
	
	xc = aD->new_xc;
	vH = gi->vertH;

	wo = gi->wo;
	ho = gi->ho;
	
	nx = aD->ta_nx;
	ny = aD->ta_ny;

	na = aD->ta_ny_13;

	nba = pola/wo;
	nbr = polr/ho;

	printf("pola: %i, polr: %i\n",pola, polr);
		

	mem_in = nx * ny *sizeof(float);
	mem_pol = pola *polr *sizeof(float);
	mem_Cos = na*sizeof(float);
	mem_fc = nx * vH * sizeof(cufftComplex);
	mem_chunk = nx * vH * sizeof(float);
	
	mb = gi->maxResidentBlocks;
	mt = gi->maxResidentThreads;
	
	maxt = mt/mb;
	
	uu = 1;
	while(maxt>1){
		uu*=2;
		maxt/=2;	
	}
	maxt = uu;
	ub = nx/maxt;

	nf = ny/vH;

	dim3 threads_1(maxt,1,1);
	dim3 threads_bp(wo,ho,1);
	
	dim3 grids_2(ub*vH,1,1);
	dim3 grid_bp(nba,nbr,1);	
	
	HMCUT_SAFE_CALL( cutCreateTimer(&hTimer) );
	HMCUT_SAFE_CALL( cutResetTimer(hTimer));
        HMCUT_SAFE_CALL( cutStartTimer(hTimer));
	
         {
		int MyDevice;
		char  devmsg[128];
                cudaGetDevice(&MyDevice);
		snprintf(devmsg,128,"my device %i",MyDevice);
		timestamp(devmsg,4);
          }
	HM_SAFE_CALL(cudaMalloc((void**)&dev_fc,mem_fc));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_inc,mem_fc));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_in,mem_in));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_pol,mem_pol));
	HM_SAFE_CALL(cudaMalloc((void**)&dev_Cos,mem_Cos));

	HM_SAFE_CALL(cudaMemcpy(dev_fc,xd->veccF,mem_fc,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(dev_in,xd->vta,mem_in,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(dev_Cos,xd->vecCos,mem_Cos,cudaMemcpyHostToDevice));

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	texA.normalized = false;
   	texA.filterMode = cudaFilterModeLinear;
   	texA.addressMode[0] = cudaAddressModeClamp;

	texCos.normalized = false;
   	texCos.filterMode = cudaFilterModeLinear;
   	texCos.addressMode[0] = cudaAddressModeClamp;

	CUDA_SAFE_CALL(cudaBindTexture(NULL, texA, dev_in, channelDesc, mem_chunk));	
	CUDA_SAFE_CALL(cudaBindTexture(NULL, texCos, dev_Cos, channelDesc, mem_Cos));	

	HMFFT_SAFE_CALL(cufftPlan1d(&plan, nx, CUFFT_C2C, vH));

	set_zero<<<grid_bp, threads_bp>>>(dev_pol, nba);
	HM_SAFE_CALL( cudaThreadSynchronize() );
	
	for(unsigned int i = 0;i<nf;i++){
	//for(unsigned int i = 0;i<1;i++){
		sp = i*nx*vH;
		cuda_data_r2c<<<grids_2,threads_1>>>(dev_inc,dev_in,sp);
		HM_SAFE_CALL( cudaThreadSynchronize() );
		HMFFT_SAFE_CALL(cufftExecC2C(plan, dev_inc, dev_inc, CUFFT_FORWARD));
		HM_SAFE_CALL( cudaThreadSynchronize() );
		cuda_mul_c<<<grids_2,threads_1>>>(dev_fc,dev_inc);
		
		HM_SAFE_CALL( cudaThreadSynchronize() );
		HMFFT_SAFE_CALL(cufftExecC2C(plan, dev_inc, dev_inc, CUFFT_INVERSE));
		HM_SAFE_CALL( cudaThreadSynchronize() );
		sp = 0;
		cuda_data_c2r<<<grids_2,threads_1>>>(dev_in,dev_inc,sp);
		HM_SAFE_CALL( cudaThreadSynchronize() );

		starty = i*vH;
		
		fbp_cpu2<<<grid_bp, threads_bp>>>(dev_in, dev_pol, dev_Cos, starty, nx, na, vH, pola, ps, xc);
		HM_SAFE_CALL( cudaThreadSynchronize() );
	}
	
	cufftDestroy(plan);
		
	HM_SAFE_CALL(cudaMemcpy(xd->vecPol, dev_pol, mem_pol, cudaMemcpyDeviceToHost));

	cudaFree(dev_fc); dev_fc = NULL;
	cudaFree(dev_inc); dev_inc = NULL;
	cudaFree(dev_in); dev_in = NULL;
	
	cudaFree(dev_pol); dev_pol = NULL;
	cudaFree(dev_Cos); dev_Cos = NULL;
	
	HMCUT_SAFE_CALL( cutStopTimer(hTimer) );
    gpuTime = cutGetTimerValue(hTimer);
	printTag(aD,"TimeBackprojection",float(gpuTime),"time (in ms)");
   	
	return true;
}

