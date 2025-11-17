// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023


#include "fftw3.h"
#include "bmIDF1_FFTW.h"
#include <cmath>

void pre_ifft1( float* in_complex_ptr_0,  float* in_real_ptr_0,  float* in_imag_ptr_0,  int Nx_u, int nCh); 
void post_ifft1(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int nCh, float F); 

void bmIDF1_FFTW(float* out_real_ptr_0, float* out_imag_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, float dKx_u, int nCh)
{
	int myNumOfThread = 32; // magic number

	float* in_real_ptr_run	= in_real_ptr_0;
	float* in_imag_ptr_run	= in_imag_ptr_0; 
	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;


	fftwf_complex* in_complex_ptr_0		= (fftwf_complex*)fftwf_alloc_complex(sizeof(fftwf_complex)* Nx_u * nCh);
	float* in_complex_ptr_run			= (float*)in_complex_ptr_0;
	fftwf_complex* out_complex_ptr_0	= (fftwf_complex*)fftwf_alloc_complex(sizeof(fftwf_complex)* Nx_u * nCh);
	float* out_complex_ptr_run			= (float*)out_complex_ptr_0;

	float F = dKx_u; // FFTW divide the inverse DFT by the array size.  
	

	// FFT ----------------------------------------------------------------
	fftwf_init_threads();
	fftwf_plan_with_nthreads(myNumOfThread);
	fftwf_plan p;

	// create plan before writing data in vector, the fftw method may write in the vecotr for testing. 
	p = fftwf_plan_dft_1d(Nx_u, in_complex_ptr_0, out_complex_ptr_0, FFTW_BACKWARD, FFTW_ESTIMATE); // n and m inverted; row versus column major order 

	pre_ifft1((float*)in_complex_ptr_0, in_real_ptr_0, in_imag_ptr_0, Nx_u, nCh);

	fftwf_execute(p); 
	if (nCh > 1)
	{
		for (int i = 0; i < nCh - 1; i++)
		{
			in_complex_ptr_run  += 2 * Nx_u;
			out_complex_ptr_run += 2 * Nx_u; 
			fftwf_execute_dft(p, (fftwf_complex*)in_complex_ptr_run, (fftwf_complex*)out_complex_ptr_run);
		}
		in_complex_ptr_run  = (float*)in_complex_ptr_0;
		out_complex_ptr_run = (float*)out_complex_ptr_0;
	}
	
	post_ifft1((float*)out_complex_ptr_0, out_real_ptr_0, out_imag_ptr_0, Nx_u, nCh, F);

	// END_FFT -----------------------------------------------------------


	


	// deallocate temp arrays --------------------------------------------------
	fftwf_free(in_complex_ptr_0); 
	fftwf_free(out_complex_ptr_0);
	fftwf_destroy_plan(p); 
	fftwf_cleanup_threads(); 
	
	// END_deallocate temp arrays --------------------------------------------------

}


void pre_ifft1(float* in_complex_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int nCh)
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


void post_ifft1(float* out_complex_ptr_0, float* out_real_ptr_0, float* out_imag_ptr_0, int Nx_u, int nCh, float F)
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
			*out_real_ptr_run++ = (*out_complex_ptr_run++) * F;
			*out_imag_ptr_run++ = (*out_complex_ptr_run++) * F;
		}

		out_complex_ptr_run -= 2*Nx_u;
		cx_max = Nx_u - cx_max;
		for (cx = 0; cx < cx_max; cx++)
		{
			*out_real_ptr_run++ = (*out_complex_ptr_run++) * F;
			*out_imag_ptr_run++ = (*out_complex_ptr_run++) * F;
		}

		out_complex_ptr_run += 2*Nx_u;

	}
	out_real_ptr_run = (float*)out_real_ptr_0;
	out_imag_ptr_run = (float*)out_imag_ptr_0;
	out_complex_ptr_run = (float*)out_complex_ptr_0;

}

