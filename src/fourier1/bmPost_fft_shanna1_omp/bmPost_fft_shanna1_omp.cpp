// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include <cmath>
#include <omp.h>
#include "bmPost_fft_shanna1_omp.h"
// #include "mex.h"

void bmPost_fft_shanna1_omp(float* out_complex_ptr_0_shared, 
							float* out_real_ptr_0_shared, 
							float* out_imag_ptr_0_shared,
							int Nx_u_shared,
							int nCh_shared)
{
	if (nCh_shared < omp_get_max_threads())
	{
		omp_set_num_threads(nCh_shared);
	}
	else
	{
		omp_set_num_threads(omp_get_max_threads());
	}

#pragma omp parallel shared(out_complex_ptr_0_shared, out_real_ptr_0_shared, out_imag_ptr_0_shared, Nx_u_shared, nCh_shared)
	{

		// mexPrintf("This is thread number %d .\n", omp_get_thread_num()); 


		float* out_real_ptr_0 = out_real_ptr_0_shared;
		float* out_real_ptr_run = out_real_ptr_0;

		float* out_imag_ptr_0 = out_imag_ptr_0_shared;
		float* out_imag_ptr_run = out_imag_ptr_0;

		float* out_complex_ptr_0 = out_complex_ptr_0_shared;
		float* out_complex_ptr_run = out_complex_ptr_0;


		int Nx_u = Nx_u_shared; 
		double Nx_u_double = (double)Nx_u;
		int cx_max;
		int cx;
		int x_incr = (int)ceil(Nx_u_double / 2);
		int nCh = nCh_shared; 


#pragma omp for
		for (int ind_ch = 0; ind_ch < nCh; ind_ch++)
		{

			out_complex_ptr_run = out_complex_ptr_0 + ((long long)2) * ((long long)x_incr) + ((long long)2)*((long long)ind_ch)*((long long)Nx_u);
			out_real_ptr_run = out_real_ptr_0 + ((long long)ind_ch)*((long long)Nx_u); 
			out_imag_ptr_run = out_imag_ptr_0 + ((long long)ind_ch)*((long long)Nx_u);

			cx_max = Nx_u - x_incr;
			for (cx = 0; cx < cx_max; cx++)
			{
				*out_real_ptr_run++ = (*out_complex_ptr_run++);
				*out_imag_ptr_run++ = (*out_complex_ptr_run++);
			}

			out_complex_ptr_run -= 2 * Nx_u;
			cx_max = Nx_u - cx_max;
			for (cx = 0; cx < cx_max; cx++)
			{
				*out_real_ptr_run++ = (*out_complex_ptr_run++);
				*out_imag_ptr_run++ = (*out_complex_ptr_run++);
			}
		}
		out_real_ptr_run = (float*)out_real_ptr_0;
		out_imag_ptr_run = (float*)out_imag_ptr_0;
		out_complex_ptr_run = (float*)out_complex_ptr_0;

	} // END_parallel
} // END_function

