// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include <cmath>
#include <omp.h>
#include "fftw3.h"
#include "bmNakatsha2_FFTW_omp.h"
#include "bmSparseMat_cC_oBlock_omp.h"
#include "bmPre_fft_nakatsha2_omp.h"
#include "bmPost_fft_nakatsha2_omp.h"

void bmNakatsha2_FFTW_omp(	float* out_real_ptr_0,
							float* out_imag_ptr_0,
							float* in_real_ptr_0,
							float* in_imag_ptr_0,
							float* KFC_conj_real_ptr_0,
							float* KFC_conj_imag_ptr_0,
							int Nx_u,
							int Ny_u,
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
	int myNumOfThread = omp_get_max_threads(); 

	float* in_real_ptr_run  = in_real_ptr_0;
	float* in_imag_ptr_run  = in_imag_ptr_0;
	float* out_real_ptr_run = out_real_ptr_0;
	float* out_imag_ptr_run = out_imag_ptr_0;

	fftwf_complex* in_complex_ptr_0 = (fftwf_complex*)fftwf_alloc_complex(sizeof(fftwf_complex)* Nx_u * Ny_u * nCh);
	float* in_complex_ptr_run = (float*)in_complex_ptr_0;
	fftwf_complex* out_complex_ptr_0 = (fftwf_complex*)fftwf_alloc_complex(sizeof(fftwf_complex)* Nx_u * Ny_u * nCh);
	float* out_complex_ptr_run = (float*)out_complex_ptr_0;

	long long half_size_64 = ((long long)Nx_u)*((long long)Ny_u)*((long long)nCh);
	float* interA_real_ptr_0 = (float*)in_complex_ptr_0;
	float* interA_imag_ptr_0 = ((float*)in_complex_ptr_0) + half_size_64;
	float* interB_real_ptr_0 =  (float*)out_complex_ptr_0;
	float* interB_imag_ptr_0 = ((float*)out_complex_ptr_0) + half_size_64;

	fftwf_init_threads();
	fftwf_plan_with_nthreads(myNumOfThread);
	fftwf_plan p;

	// create plan before writing data in vector, the fftw method may write in the vecotr for testing. 
	p = fftwf_plan_dft_2d(Ny_u, Nx_u, in_complex_ptr_0, out_complex_ptr_0, FFTW_BACKWARD, FFTW_ESTIMATE); // n and m inverted; row versus column major order


	// FFT ----------------------------------------------------------------
	


	// gridding ----------------------------------------------------------
	bmSparseMat_cC_oBlock_omp(r_size, r_jump_ptr_0, r_nJump_ptr_0,
		m_val_ptr_0,
		l_size, l_jump_ptr_0, l_nJump,
		in_real_ptr_0, in_imag_ptr_0,
		nCh,
		interB_real_ptr_0, interB_imag_ptr_0);
	// END_gridding ------------------------------------------------------

	

	
	bmPre_fft_nakatsha2_omp((float*)in_complex_ptr_0, interB_real_ptr_0, interB_imag_ptr_0, Nx_u, Ny_u, nCh);
	fftwf_execute(p);
	if (nCh > 1)
	{
		for (int i = 0; i < nCh - 1; i++)
		{
			in_complex_ptr_run += 2 * Nx_u * Ny_u;
			out_complex_ptr_run += 2 * Nx_u * Ny_u;
			fftwf_execute_dft(p, (fftwf_complex*)in_complex_ptr_run, (fftwf_complex*)out_complex_ptr_run);
		}
		in_complex_ptr_run = (float*)in_complex_ptr_0;
		out_complex_ptr_run = (float*)out_complex_ptr_0;
	}
	bmPost_fft_nakatsha2_omp((float*)out_complex_ptr_0, out_real_ptr_0, out_imag_ptr_0, interA_real_ptr_0, interA_imag_ptr_0, KFC_conj_real_ptr_0, KFC_conj_imag_ptr_0, Nx_u, Ny_u, nCh);
	// END_FFT -----------------------------------------------------------

	


	// deallocate temp arrays --------------------------------------------------
	fftwf_free(in_complex_ptr_0);
	fftwf_free(out_complex_ptr_0);
	fftwf_destroy_plan(p);
	fftwf_cleanup_threads();

	// END_deallocate temp arrays --------------------------------------------------


} // end function

