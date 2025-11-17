// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include <cuda_runtime.h>
#include <cufft.h>
#include <cmath>
#include "bmNakatsha1_CUFFT_omp.h"
#include "bmSparseMat_cC_oBlock_omp.h"
#include "bmPre_fft_nakatsha1_omp.h"
#include "bmPost_fft_nakatsha1_omp.h"

void bmNakatsha1_CUFFT_omp(	float* out_real_ptr_0,
							float* out_imag_ptr_0,
							float* in_real_ptr_0,
							float* in_imag_ptr_0,
							float* KFC_conj_real_ptr_0,
							float* KFC_conj_imag_ptr_0,
							int Nx_u,
							int nCh,
							int r_size,
							int* r_jump_ptr_0,
							int* r_nJump_ptr_0,
							float* m_val_ptr_0,
							int l_size,
							int* l_jump_ptr_0,
							int l_nJump
							)
{
	
	float* in_real_ptr_run = in_real_ptr_0;
	float* in_imag_ptr_run = in_imag_ptr_0;
	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;

	float* in_complex_ptr_0 = new float[2 * Nx_u * nCh];
	float* in_complex_ptr_run = (float*)in_complex_ptr_0;
	float* out_complex_ptr_0 = new float[2 * Nx_u * nCh];
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;

	long long half_size_64   = ((long long)Nx_u)*((long long)nCh); 
	float* interA_real_ptr_0 = (float*)in_complex_ptr_0;
	float* interA_imag_ptr_0 = ((float*)in_complex_ptr_0) + half_size_64;
	float* interB_real_ptr_0 =  (float*)out_complex_ptr_0;
	float* interB_imag_ptr_0 = ((float*)out_complex_ptr_0) + half_size_64;


	cufftHandle plan;
	cufftPlan1d(&plan, Nx_u, CUFFT_C2C, 1);

	cufftComplex* device_in;
	cudaMalloc(&device_in, sizeof(cufftComplex)* Nx_u);
	cufftComplex* device_out;
	cudaMalloc(&device_out, sizeof(cufftComplex)* Nx_u);
	

	// FFT ----------------------------------------------------------------
	


	// gridding ----------------------------------------------------------
	bmSparseMat_cC_oBlock_omp(r_size, r_jump_ptr_0, r_nJump_ptr_0,
		m_val_ptr_0,
		l_size, l_jump_ptr_0, l_nJump,
		in_real_ptr_0, in_imag_ptr_0,
		nCh,
		interB_real_ptr_0, interB_imag_ptr_0);
	// END_gridding ------------------------------------------------------

	

	
	bmPre_fft_nakatsha1_omp((float*)in_complex_ptr_0, interB_real_ptr_0, interB_imag_ptr_0, Nx_u, nCh);
	for (int i = 0; i < nCh; i++)
	{
		cudaMemcpy(device_in, in_complex_ptr_run, Nx_u * sizeof(cufftComplex), cudaMemcpyHostToDevice);
		cufftExecC2C(plan, device_in, device_out, CUFFT_INVERSE);
		cudaDeviceSynchronize();
		cudaMemcpy(out_complex_ptr_run, device_out, Nx_u * sizeof(cufftComplex), cudaMemcpyDeviceToHost);
		in_complex_ptr_run += 2 * Nx_u;
		out_complex_ptr_run += 2 * Nx_u;
	}
	in_complex_ptr_run = in_complex_ptr_0;
	out_complex_ptr_run = out_complex_ptr_0;
	bmPost_fft_nakatsha1_omp((float*)out_complex_ptr_0, out_real_ptr_0, out_imag_ptr_0, interA_real_ptr_0, interA_imag_ptr_0, KFC_conj_real_ptr_0, KFC_conj_imag_ptr_0, Nx_u, nCh);

	// END_FFT -----------------------------------------------------------

	


	// deallocate temp arrays --------------------------------------------------
	cufftDestroy(plan);
	cudaFree(device_in);
	cudaFree(device_out);
	delete[] in_complex_ptr_0;
	delete[] out_complex_ptr_0;
	// END_deallocate temp arrays --------------------------------------------------


} // end function
