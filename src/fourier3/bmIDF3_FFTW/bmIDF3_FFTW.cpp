// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023


#include "fftw3.h"
#include "bmIDF3_FFTW.h"
#include <cmath>

void pre_ifft3(float* in_complex_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, int Nz_u, int nCh);
void post_ifft3(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int Ny_u, int Nz_u, int nCh, float F);

void bmIDF3_FFTW(float* out_real_ptr_0, float* out_imag_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, int Nz_u, float dKx_u, float dKy_u, float dKz_u, int nCh)
{
	int myNumOfThread = 32; // magic_number

	float* in_real_ptr_run = in_real_ptr_0;
	float* in_imag_ptr_run = in_imag_ptr_0;
	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;


	fftwf_complex* in_complex_ptr_0 = (fftwf_complex*)fftwf_alloc_complex(sizeof(fftwf_complex)* Nx_u * Ny_u * Nz_u * nCh);
	float* in_complex_ptr_run = (float*)in_complex_ptr_0;
	fftwf_complex* out_complex_ptr_0 = (fftwf_complex*)fftwf_alloc_complex(sizeof(fftwf_complex)* Nx_u * Ny_u * Nz_u * nCh);
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;

	float F = dKx_u*dKy_u*dKz_u; // scaling is different between FFTW and MATLAB



	// FFT ----------------------------------------------------------------
	fftwf_init_threads();
	fftwf_plan_with_nthreads(myNumOfThread);
	fftwf_plan p;

	// Create plan before writing data in vector, the fftw planing method may write in the vecor for testing.
	// Last dim must have fastest varying index (largest rapidity). 
	// But I don't see any logic in this order, exept rapididy of index. 
	p = fftwf_plan_dft_3d(Nz_u, Ny_u, Nx_u, in_complex_ptr_0, out_complex_ptr_0, FFTW_BACKWARD, FFTW_ESTIMATE); 


	// we copy the data in complex_interleaved order and do ifft shift at the same time
	pre_ifft3((float*)in_complex_ptr_0, in_real_ptr_0, in_imag_ptr_0, Nx_u, Ny_u, Nz_u, nCh);

	fftwf_execute(p);
	if (nCh > 1)
	{
		for (int i = 0; i < nCh - 1; i++)
		{
			in_complex_ptr_run += 2 * Nx_u*Ny_u*Nz_u;
			out_complex_ptr_run += 2 * Nx_u*Ny_u*Nz_u;
			fftwf_execute_dft(p, (fftwf_complex*)in_complex_ptr_run, (fftwf_complex*)out_complex_ptr_run);
		}
		in_complex_ptr_run = (float*)in_complex_ptr_0;
		out_complex_ptr_run = (float*)out_complex_ptr_0;
	}

	// we copy the data back in complex_sequential order and do fft shift and F multiply at the same time
	post_ifft3((float*)out_complex_ptr_0, out_real_ptr_0, out_imag_ptr_0, Nx_u, Ny_u, Nz_u, nCh, F);

	// END_FFT -----------------------------------------------------------





	// deallocate temp arrays --------------------------------------------------
	fftwf_free(in_complex_ptr_0);
	fftwf_free(out_complex_ptr_0);
	fftwf_destroy_plan(p); 
	fftwf_cleanup_threads(); 
	// END_deallocate temp arrays --------------------------------------------------

}


void pre_ifft3(float* in_complex_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, int Nz_u, int nCh)
{


	float* in_real_ptr_run = in_real_ptr_0;
	float* in_imag_ptr_run = in_imag_ptr_0;
	float* in_complex_ptr_run = in_complex_ptr_0;


	double Nx_u_double = (double)Nx_u;
	double Ny_u_double = (double)Ny_u;
	double Nz_u_double = (double)Nz_u;

	int Nx_half = (int)floor(Nx_u_double / 2);
	int Ny_half = (int)floor(Ny_u_double / 2);
	int Nz_half = (int)floor(Nz_u_double / 2);

	int cx;
	int cx_max;
	int cy;
	int cy_max;
	int cz;
	int cz_max;

	int ch_ind;



	in_real_ptr_run += (Nx_half + Nx_u*Ny_half + Nx_u*Ny_u*Nz_half);
	in_imag_ptr_run += (Nx_half + Nx_u*Ny_half + Nx_u*Ny_u*Nz_half);
	for (ch_ind = 0; ch_ind < nCh; ch_ind++)
	{
		cz_max = Nz_u - Nz_half;
		for (cz = 0; cz < cz_max; cz++)
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


		in_real_ptr_run -= Nx_u*Ny_u*Nz_u;
		in_imag_ptr_run -= Nx_u*Ny_u*Nz_u;
		cz_max = Nz_u - cz_max;
		for (cz = 0; cz < cz_max; cz++)
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

		in_real_ptr_run += Nx_u*Ny_u*Nz_u;
		in_imag_ptr_run += Nx_u*Ny_u*Nz_u;

	}
	in_real_ptr_run = in_real_ptr_0;
	in_imag_ptr_run = in_imag_ptr_0;
	in_complex_ptr_run = in_complex_ptr_0;

}


void post_ifft3(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int Ny_u, int Nz_u, int nCh, float F)
{

	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;

	double Nx_u_double = (double)Nx_u;
	double Ny_u_double = (double)Ny_u;
	double Nz_u_double = (double)Nz_u;

	int Nx_half = (int)ceil(Nx_u_double / 2);
	int Ny_half = (int)ceil(Ny_u_double / 2);
	int Nz_half = (int)ceil(Nz_u_double / 2);

	int cx;
	int cx_max;
	int cy;
	int cy_max;
	int cz;
	int cz_max;
	
	int ch_ind;


	out_complex_ptr_run += 2 * (Nx_half + Nx_u*Ny_half + Nx_u*Ny_u*Nz_half);
	for (ch_ind = 0; ch_ind < nCh; ch_ind++)
	{

		cz_max = Nz_u - Nz_half;
		for (cz = 0; cz < cz_max; cz++)
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

		out_complex_ptr_run -= 2 * Nx_u*Ny_u*Nz_u;
		cz_max = Nz_u - cz_max;
		for (cz = 0; cz < cz_max; cz++)
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
		out_complex_ptr_run += 2 * Nx_u * Ny_u * Nz_u;
	}
	out_real_ptr_run = out_real_ptr_0;
	out_imag_ptr_run = out_imag_ptr_0;
	out_complex_ptr_run = out_complex_ptr_0;
}

