// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "bmDFT1_CUFFT.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <cmath>

void pre_fft1( float* in_complex_ptr_0,  float* in_real_ptr_0,  float* in_imag_ptr_0,  int Nx_u, int nCh); 
void post_fft1(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int nCh, float F); 

void bmDFT1_CUFFT(float* out_real_ptr_0, float* out_imag_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, float dKx_u, int nCh)
{

	float* in_real_ptr_run	= in_real_ptr_0;
	float* in_imag_ptr_run	= in_imag_ptr_0; 
	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;

	float* in_complex_ptr_0		= new float[2 * Nx_u * nCh];
	float* in_complex_ptr_run	= (float*)in_complex_ptr_0;
	float* out_complex_ptr_0	= new float[2 * Nx_u * nCh];
	float* out_complex_ptr_run	= (float*)out_complex_ptr_0;
	

	// FFT ----------------------------------------------------------------

	float F = ((float)Nx_u)*dKx_u;

	cufftHandle plan;
	cufftPlan1d(&plan, Nx_u, CUFFT_C2C, 1);

	cufftComplex* device_in;
	cudaMalloc(&device_in, sizeof(cufftComplex)* Nx_u);
	cufftComplex* device_out;
	cudaMalloc(&device_out, sizeof(cufftComplex)* Nx_u);


	pre_fft1(in_complex_ptr_0, in_real_ptr_0, in_imag_ptr_0, Nx_u, nCh);
	for (int i = 0; i < nCh; i++)
	{
		cudaMemcpy(device_in, in_complex_ptr_run, Nx_u * sizeof(cufftComplex), cudaMemcpyHostToDevice);
		cufftExecC2C(plan, device_in, device_out, CUFFT_FORWARD);
		cudaDeviceSynchronize();
		cudaMemcpy(out_complex_ptr_run, device_out, Nx_u * sizeof(cufftComplex), cudaMemcpyDeviceToHost);
		in_complex_ptr_run += 2 * Nx_u;
		out_complex_ptr_run += 2 * Nx_u;
	}
	in_complex_ptr_run = in_complex_ptr_0;
	out_complex_ptr_run = out_complex_ptr_0;
	post_fft1(out_complex_ptr_0, out_real_ptr_0, out_imag_ptr_0, Nx_u, nCh, F);

	// END_FFT -----------------------------------------------------------


	


	// deallocate temp arrays --------------------------------------------------

	cufftDestroy(plan);
	cudaFree(device_in);
	cudaFree(device_out); 
	delete[] in_complex_ptr_0;
	delete[] out_complex_ptr_0;
	
	// END_deallocate temp arrays --------------------------------------------------

}


void pre_fft1(float* in_complex_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int nCh)
{


	float* in_real_ptr_run		= in_real_ptr_0;
	float* in_imag_ptr_run		= in_imag_ptr_0;
	float* in_complex_ptr_run	= (float*)in_complex_ptr_0;


	double Nx_u_double = (double)Nx_u;
	int cx_max;
	int cx;
	int ind_ch;
	int x_incr = (int)floor(Nx_u_double / 2);



	in_real_ptr_run += x_incr;
	in_imag_ptr_run += x_incr;
	for (ind_ch = 0; ind_ch < nCh; ind_ch++)
	{
		cx_max = Nx_u - x_incr;
		for (cx = 0; cx < cx_max; cx++)
		{
			*in_complex_ptr_run++ = *in_real_ptr_run++;
			*in_complex_ptr_run++ = *in_imag_ptr_run++;
		}

		in_real_ptr_run -= Nx_u;
		in_imag_ptr_run -= Nx_u;
		cx_max = Nx_u - cx_max;
		for (cx = 0; cx < cx_max; cx++)
		{
			*in_complex_ptr_run++ = *in_real_ptr_run++;
			*in_complex_ptr_run++ = *in_imag_ptr_run++;
		}

		in_real_ptr_run += Nx_u;
		in_imag_ptr_run += Nx_u;

	}
	in_real_ptr_run = (float*)in_real_ptr_0;
	in_imag_ptr_run = (float*)in_imag_ptr_0;
	in_complex_ptr_run = (float*)in_complex_ptr_0;

}


void post_fft1(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int nCh, float F)
{

	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;


	double Nx_u_double = (double)Nx_u;
	int cx_max;
	int cx;
	int ind_ch;
	int x_incr = (int)ceil(Nx_u_double / 2);



	out_complex_ptr_run += 2*x_incr;
	for (ind_ch = 0; ind_ch < nCh; ind_ch++)
	{
		cx_max = Nx_u - x_incr;
		for (cx = 0; cx < cx_max; cx++)
		{
			*out_real_ptr_run++ = (*out_complex_ptr_run++) / F;
			*out_imag_ptr_run++ = (*out_complex_ptr_run++) / F;
		}

		out_complex_ptr_run -= 2*Nx_u;
		cx_max = Nx_u - cx_max;
		for (cx = 0; cx < cx_max; cx++)
		{
			*out_real_ptr_run++ = (*out_complex_ptr_run++) / F;
			*out_imag_ptr_run++ = (*out_complex_ptr_run++) / F;
		}

		out_complex_ptr_run += 2*Nx_u;

	}
	out_real_ptr_run = (float*)out_real_ptr_0;
	out_imag_ptr_run = (float*)out_imag_ptr_0;
	out_complex_ptr_run = (float*)out_complex_ptr_0;

}

