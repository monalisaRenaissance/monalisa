// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "bmIDF2_CUFFT.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <cmath>

void pre_ifft2(float* in_complex_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, int nCh);
void post_ifft2(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int Ny_u, int nCh, float F);

void bmIDF2_CUFFT(float* out_real_ptr_0, float* out_imag_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, float dKx_u, float dKy_u, int nCh)
{

	float* in_real_ptr_run = in_real_ptr_0;
	float* in_imag_ptr_run = in_imag_ptr_0;
	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;

	float* in_complex_ptr_0 = new float[2 * Nx_u * Ny_u * nCh];
	float* in_complex_ptr_run = (float*)in_complex_ptr_0;
	float* out_complex_ptr_0 = new float[2 * Nx_u * Ny_u * nCh];
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;
	
	cufftHandle plan;
	cufftPlan2d(&plan, Ny_u, Nx_u, CUFFT_C2C);

	cufftComplex* device_in;
	cudaMalloc(&device_in, sizeof(cufftComplex) * Nx_u * Ny_u);
	cufftComplex* device_out;
	cudaMalloc(&device_out, sizeof(cufftComplex) * Nx_u * Ny_u);

	// FFT ----------------------------------------------------------------

	float F = dKx_u*dKy_u;  


	// we copy the data in complex_interleaved order and do ifft shift at the same time
	pre_ifft2((float*)in_complex_ptr_0, in_real_ptr_0, in_imag_ptr_0, Nx_u, Ny_u, nCh);
	for (int i = 0; i < nCh; i++)
	{
		cudaMemcpy(device_in, in_complex_ptr_run, Nx_u * Ny_u * sizeof(cufftComplex), cudaMemcpyHostToDevice);
		cufftExecC2C(plan, device_in, device_out, CUFFT_INVERSE);
		cudaDeviceSynchronize();
		cudaMemcpy(out_complex_ptr_run, device_out, Nx_u * Ny_u * sizeof(cufftComplex), cudaMemcpyDeviceToHost);
		in_complex_ptr_run += 2 * Nx_u * Ny_u;
		out_complex_ptr_run += 2 * Nx_u * Ny_u;
	}
	in_complex_ptr_run = in_complex_ptr_0;
	out_complex_ptr_run = out_complex_ptr_0;
	// we copy the data back in complex_sequential order and do fft shift and F multiply at the same time
	post_ifft2((float*)out_complex_ptr_0, out_real_ptr_0, out_imag_ptr_0, Nx_u, Ny_u, nCh, F);

	// END_FFT -----------------------------------------------------------

	// deallocate temp arrays --------------------------------------------------
	cufftDestroy(plan);
	cudaFree(device_in);
	cudaFree(device_out);
	delete[] in_complex_ptr_0;
	delete[] out_complex_ptr_0;
	// END_deallocate temp arrays --------------------------------------------------

}


void pre_ifft2(float* in_complex_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, int nCh)
{


	float* in_real_ptr_run = in_real_ptr_0;
	float* in_imag_ptr_run = in_imag_ptr_0;
	float* in_complex_ptr_run = in_complex_ptr_0;


	double Nx_u_double = (double)Nx_u;
	double Ny_u_double = (double)Ny_u;

	int Nx_half = (int)floor(Nx_u_double / 2);
	int Ny_half = (int)floor(Ny_u_double / 2);

	int cx_max;
	int cx;
	int cy_max;
	int cy;
	int ch_ind;


	in_real_ptr_run += (Nx_half + Nx_u*Ny_half);
	in_imag_ptr_run += (Nx_half + Nx_u*Ny_half);
	for (ch_ind = 0; ch_ind < nCh; ch_ind++)
	{
		cy_max = Ny_u - Ny_half;
		for (cy = 0; cy < cy_max; cy++)
		{

			cx_max = Nx_u - Nx_half;
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

		in_real_ptr_run -= Nx_u*Ny_u;
		in_imag_ptr_run -= Nx_u*Ny_u;
		cy_max = Ny_u - cy_max;
		for (cy = 0; cy < cy_max; cy++)
		{

			cx_max = Nx_u - Nx_half;
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

		in_real_ptr_run += Nx_u*Ny_u;
		in_imag_ptr_run += Nx_u*Ny_u;

	}
	in_real_ptr_run = in_real_ptr_0;
	in_imag_ptr_run = in_imag_ptr_0;
	in_complex_ptr_run = in_complex_ptr_0;

}


void post_ifft2(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int Ny_u, int nCh, float F)
{

	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;

	double Nx_u_double = (double)Nx_u;
	double Ny_u_double = (double)Ny_u;

	int Nx_half = (int)ceil(Nx_u_double / 2);
	int Ny_half = (int)ceil(Ny_u_double / 2);

	int cx_max;
	int cx;
	int cy_max;
	int cy;
	int ch_ind;


	out_complex_ptr_run += 2 * (Nx_half + Nx_u*Ny_half);
	for (ch_ind = 0; ch_ind < nCh; ch_ind++)
	{
		cy_max = Ny_u - Ny_half;
		for (cy = 0; cy < cy_max; cy++)
		{

			cx_max = Nx_u - Nx_half;
			for (cx = 0; cx < cx_max; cx++)
			{
				*out_real_ptr_run++ = (*out_complex_ptr_run++) * F;
				*out_imag_ptr_run++ = (*out_complex_ptr_run++) * F;
			}

			out_complex_ptr_run -= 2 * Nx_u;
			cx_max = Nx_u - cx_max;
			for (cx = 0; cx < cx_max; cx++)
			{
				*out_real_ptr_run++ = (*out_complex_ptr_run++) * F;
				*out_imag_ptr_run++ = (*out_complex_ptr_run++) * F;
			}

			out_complex_ptr_run += 2 * Nx_u;
		}

		out_complex_ptr_run -= 2 * Nx_u*Ny_u;
		cy_max = Ny_u - cy_max;
		for (cy = 0; cy < cy_max; cy++)
		{

			cx_max = Nx_u - Nx_half;
			for (cx = 0; cx < cx_max; cx++)
			{
				*out_real_ptr_run++ = (*out_complex_ptr_run++) * F;
				*out_imag_ptr_run++ = (*out_complex_ptr_run++) * F;
			}

			out_complex_ptr_run -= 2 * Nx_u;
			cx_max = Nx_u - cx_max;
			for (cx = 0; cx < cx_max; cx++)
			{
				*out_real_ptr_run++ = (*out_complex_ptr_run++) * F;
				*out_imag_ptr_run++ = (*out_complex_ptr_run++) * F;
			}

			out_complex_ptr_run += 2 * Nx_u;
		}

		out_complex_ptr_run += 2 * Nx_u * Ny_u;

	}
	out_real_ptr_run = out_real_ptr_0;
	out_imag_ptr_run = out_imag_ptr_0;
	out_complex_ptr_run = out_complex_ptr_0;

}

