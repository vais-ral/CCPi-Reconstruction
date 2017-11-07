#include "defs_gpu.cuh"


__global__ void find_cart(float * pc_res_d, float *pc_a_d, float *pc_r_d, int cart_len){
	int start = blockIdx.x*cart_len;
	for(int i=start;i<(start+cart_len);i++){
		//pc_res_d[i]= tex2D(texPC,pc_a_d[i],pc_r_d[i]);
		pc_res_d[i]= tex2D(texPC,pc_a_d[i]+0.5f,pc_r_d[i]+0.5f);
	}
}

__global__ void find_fdk2(float* vecv_d, float *m_x_d, float *m_y_d, int *i_a_d, int *i_r_d, int im, int na_old, int nrn, int aimin, int w_a, int w_r, size_t pitch_pc){
	int rr1, aa1;
		
	rr1 = (i_r_d[blockIdx.x]*blockDim.y+threadIdx.y)*w_r;
	aa1 = (i_a_d[blockIdx.x]*blockDim.x+threadIdx.x)*w_a;


	int rr, aa, a2, m1, m2;
	int p4 = pitch_pc/4;
	
	for(int j=0;j<w_a;j++){	
		aa = aa1+j;
		a2 = (aa+im+aimin)%(na_old);
		for(int i=0;i<w_r;i++){
			rr = rr1+i;
			m1 = rr*p4+aa;
			m2 = rr*na_old+a2;
			//vecv_d[m1]+= tex2D(texI,m_x_d[m2],m_y_d[m2]);
			vecv_d[m1]+= tex2D(texI, m_x_d[m2]+0.5f, m_y_d[m2]+0.5f);
			//vecv_d[m1] = m2;
			//vecv_d[m1]= m1;//aa1;
			//vecv_d[m1]= rr1;//1.2+blockDim.y;//i_r_d[blockIdx.x];//tex2D(texI,m_x_d[m2],m_y_d[m2]);
		}
	}
	
}


__global__ void find_fdk3(float* vecv_d, float *m_x_d, float *m_y_d, int *i_a_d, int *i_r_d, int im, int na_old, int nrn, int aimin, int w_a, int w_r, size_t pitch_pc, float vx, float vy){
	int rr1, aa1;
		
	rr1 = (i_r_d[blockIdx.x]*blockDim.y+threadIdx.y)*w_r;
	aa1 = (i_a_d[blockIdx.x]*blockDim.x+threadIdx.x)*w_a;


	int rr, aa, a2, m1, m2;
	int p4 = pitch_pc/4;
	
	for(int j=0;j<w_a;j++){	
		aa = aa1+j;
		a2 = (aa+im+aimin)%(na_old);
		for(int i=0;i<w_r;i++){
			rr = rr1+i;
			m1 = rr*p4+aa;
			m2 = rr*na_old+a2;
			//vecv_d[m1]+= tex2D(texI,m_x_d[m2],m_y_d[m2]);
			vecv_d[m1]+= tex2D(texI, m_x_d[m2]+0.5f-vx, m_y_d[m2]+0.5f-vy);
			//vecv_d[m1] = m2;
			//vecv_d[m1]= m1;//aa1;
			//vecv_d[m1]= rr1;//1.2+blockDim.y;//i_r_d[blockIdx.x];//tex2D(texI,m_x_d[m2],m_y_d[m2]);
		}
	}
	
}



bool cudaReconstructFDK(MatPar *matp, PolarToCart *PC, allData *aD, OtherParam *other){
	int nx, ny;
	nx = aD->nx;
	ny = matp->y_maxi - matp->y_mini+1;

	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	cudaSetDevice(0);

	int ustep3 = 32;
	dim3 grid_cr(nx/ustep3,ny);
	dim3 grid_cr2(1,ny);

	int volMatrix;
	volMatrix = matp->nr*matp->na;
	printf("volMatrix: %i\n",volMatrix);

	dim3 gfdk2(PC->ntb,1,1);
	dim3 tfdk2(b_a,b_r);
	
	float *veci_d;
	float *m_x_d, *m_y_d;

	int mem_allm = volMatrix*sizeof(float);

	HM_SAFE_CALL(cudaMalloc((void**)&m_x_d,mem_allm));
	HM_SAFE_CALL(cudaMalloc((void**)&m_y_d,mem_allm));
	HM_SAFE_CALL(cudaMemcpy(m_x_d,matp->m_x,mem_allm,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(m_y_d,matp->m_y,mem_allm,cudaMemcpyHostToDevice));

	size_t pitch_pc;
	float *vecv_d;
	HM_SAFE_CALL(cudaMallocPitch((void**)&vecv_d, &pitch_pc, PC->nan * sizeof(float), PC->nrn+1));
	HM_SAFE_CALL(cudaMemset2D(vecv_d, pitch_pc, 0, PC->nan*sizeof(float), PC->nrn+1));

	texPC.normalized = false;
	texPC.filterMode = cudaFilterModeLinear;
	texPC.addressMode[0] = cudaAddressModeClamp;
   	texPC.addressMode[1] = cudaAddressModeClamp;
	cudaChannelFormatDesc channelDesc_pc = cudaCreateChannelDesc<float>();
	cudaBindTexture2D(NULL, &texPC, vecv_d, &channelDesc_pc, PC->nan, PC->nrn+1, pitch_pc);

	printf("PC_nan: %i\n",PC->nan);
	printf("PC_nrn: %i\n",PC->nrn);
	printf("pitch_pc: %i\n",pitch_pc);

	int mem_size_ntb = PC->ntb*sizeof(Ipp32s);

	int *i_a_d, *i_r_d;
	HM_SAFE_CALL(cudaMalloc((void**)&i_a_d,mem_size_ntb));
	HM_SAFE_CALL(cudaMalloc((void**)&i_r_d,mem_size_ntb));
	HM_SAFE_CALL(cudaMemcpy(i_a_d,PC->i_a,mem_size_ntb,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(i_r_d,PC->i_r,mem_size_ntb,cudaMemcpyHostToDevice));

	size_t pitch;
		
	HM_SAFE_CALL(cudaMallocPitch((void**)&veci_d, &pitch, nx * sizeof(float), ny));

	texI.normalized = false;
	texI.filterMode = cudaFilterModeLinear;
   	texI.addressMode[0] = cudaAddressModeClamp;
   	texI.addressMode[1] = cudaAddressModeClamp;

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
	cudaBindTexture2D(NULL, &texI, veci_d, &channelDesc, nx, ny, pitch);
	
	for(int i=0;i<other->images_to_process;i++){
		HM_SAFE_CALL(cudaMemcpy2D(veci_d, pitch, matp->ivol+i*nx*ny,nx*sizeof(float),nx*sizeof(float),ny,cudaMemcpyHostToDevice));		
		find_fdk3<<<gfdk2,tfdk2>>>(vecv_d, m_x_d, m_y_d, i_a_d, i_r_d, i, matp->na, PC->nrn, PC->aimin, PC->w_a, PC->w_r, pitch_pc, aD->data->shift_x[i], aD->data->shift_y[i]);
		//find_fdk2<<<gfdk2,tfdk2>>>(vecv_d, m_x_d, m_y_d, i_a_d, i_r_d, i, matp->na, PC->nrn, PC->aimin, PC->w_a, PC->w_r, pitch_pc);
		HM_SAFE_CALL( cudaThreadSynchronize() );
	}

	int mem_cart = PC->pc_len*sizeof(float);

	cudaUnbindTexture(&texI);
	HM_SAFE_CALL(cudaFree(veci_d)); veci_d = NULL;
	HM_SAFE_CALL(cudaFree(m_x_d)); m_x_d = NULL;
	HM_SAFE_CALL(cudaFree(m_y_d)); m_y_d = NULL;
	HM_SAFE_CALL(cudaFree(i_a_d)); i_a_d = NULL;
	HM_SAFE_CALL(cudaFree(i_r_d)); i_r_d = NULL;


	int cart_len = 256;
	dim3 tcart(1,1,1);
	dim3 gcart(PC->pc_len/cart_len,1,1);
	
	float *pc_r_d, *pc_a_d, *pc_res_d;
	
	HM_SAFE_CALL(cudaMalloc((void**)&pc_r_d,mem_cart));
	HM_SAFE_CALL(cudaMalloc((void**)&pc_a_d,mem_cart));
	HM_SAFE_CALL(cudaMalloc((void**)&pc_res_d,mem_cart));
	HM_SAFE_CALL(cudaMemcpy(pc_r_d,PC->pc_r,mem_cart,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(pc_a_d,PC->pc_a,mem_cart,cudaMemcpyHostToDevice));
	HM_SAFE_CALL(cudaMemcpy(pc_res_d,PC->pc_res,mem_cart,cudaMemcpyHostToDevice));
		

	find_cart<<<gcart,tcart>>>(pc_res_d,pc_a_d,pc_r_d,cart_len);
	HM_SAFE_CALL( cudaThreadSynchronize() );
	HM_SAFE_CALL(cudaMemcpy(PC->pc_res,pc_res_d,mem_cart,cudaMemcpyDeviceToHost));
	HM_SAFE_CALL( cudaThreadSynchronize() );
	cudaUnbindTexture(&texPC);

	ippsZero_32f(matp->rec,PC->OutputHeight * PC->OutputWidth);

	int sp;
	sp = 0;
	
	for(int j=0;j<PC->pc_nc;j++){
		ippsCopy_32f(PC->pc_res+sp,matp->rec+PC->pc_pos[j],PC->pc_size[j]);
		sp+=(PC->pc_size[j]);
	}
	
	HM_SAFE_CALL(cudaFree(vecv_d)); vecv_d = NULL;
	HM_SAFE_CALL(cudaFree(pc_r_d)); pc_r_d = NULL;
	HM_SAFE_CALL(cudaFree(pc_a_d)); pc_a_d = NULL;
	HM_SAFE_CALL(cudaFree(pc_res_d)); pc_res_d = NULL;

	return true;
}






