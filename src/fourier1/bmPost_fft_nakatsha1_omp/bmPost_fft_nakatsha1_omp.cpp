// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include <cmath>
#include <omp.h>
#include "bmPost_fft_nakatsha1_omp.h"

void bmPost_fft_nakatsha1_omp(	float* out_complex_ptr_0_shared,
								float* out_real_ptr_0_shared,
								float* out_imag_ptr_0_shared,
								float* inter_real_ptr_0_shared,
								float* inter_imag_ptr_0_shared,
								float* KFC_conj_real_ptr_0_shared,
								float* KFC_conj_imag_ptr_0_shared,
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

#pragma omp parallel shared(out_complex_ptr_0_shared, \
							out_real_ptr_0_shared, \
							out_imag_ptr_0_shared, \
							inter_real_ptr_0_shared, \
							inter_imag_ptr_0_shared, \
							KFC_conj_real_ptr_0_shared, \
							KFC_conj_imag_ptr_0_shared, \
							Nx_u_shared, \
							nCh_shared)
	{

		float* out_real_ptr_0 = out_real_ptr_0_shared;
		float* out_real_ptr_run = out_real_ptr_0;

		float* out_imag_ptr_0 = out_imag_ptr_0_shared;
		float* out_imag_ptr_run = out_imag_ptr_0;

		float* inter_real_ptr_0 = inter_real_ptr_0_shared;
		float* inter_real_ptr_run = inter_real_ptr_0;

		float* inter_imag_ptr_0 = inter_imag_ptr_0_shared;
		float* inter_imag_ptr_run = inter_imag_ptr_0;

		float* KFC_conj_real_ptr_0 = KFC_conj_real_ptr_0_shared;
		float* KFC_conj_real_ptr_run = KFC_conj_real_ptr_0;

		float* KFC_conj_imag_ptr_0 = KFC_conj_imag_ptr_0_shared;
		float* KFC_conj_imag_ptr_run = KFC_conj_imag_ptr_0;

		float* out_complexR_ptr_0 = (float*)out_complex_ptr_0_shared;
		float* out_complexR_ptr_run = (float*)out_complexR_ptr_0;

		float* out_complexI_ptr_0 = out_complexR_ptr_0 + 1;
		float* out_complexI_ptr_run = (float*)out_complexI_ptr_0;

		int nCh = nCh_shared;

		int Nx_u = Nx_u_shared;
		double Nx_u_double = (double)Nx_u;
		int cx_max;
		int cx;
		int x_half = (int)ceil(Nx_u_double / 2);

		long long start_incr = ((long long)2) * ((long long)x_half);

		float F = 1.0f / ((float)Nx_u);

		int ind_ch;

#pragma omp for
		for (ind_ch = 0; ind_ch < nCh; ind_ch++)
		{
			inter_real_ptr_run		= inter_real_ptr_0		+ ((long long)ind_ch)*((long long)Nx_u);
			inter_imag_ptr_run		= inter_imag_ptr_0		+ ((long long)ind_ch)*((long long)Nx_u);
			KFC_conj_real_ptr_run	= KFC_conj_real_ptr_0	+ ((long long)ind_ch)*((long long)Nx_u);
			KFC_conj_imag_ptr_run	= KFC_conj_imag_ptr_0	+ ((long long)ind_ch)*((long long)Nx_u);
			out_complexR_ptr_run	= out_complexR_ptr_0	+ start_incr + ((long long)2)*((long long)ind_ch)*((long long)Nx_u);
			out_complexI_ptr_run	= out_complexR_ptr_run	+ 1;

			cx_max = Nx_u - x_half;
			for (cx = 0; cx < cx_max; cx++)
			{
				*inter_real_ptr_run++ = ((*out_complexR_ptr_run)*(*KFC_conj_real_ptr_run) - (*out_complexI_ptr_run)*(*KFC_conj_imag_ptr_run))*F;
				*inter_imag_ptr_run++ = ((*out_complexR_ptr_run)*(*KFC_conj_imag_ptr_run++) + (*out_complexI_ptr_run)*(*KFC_conj_real_ptr_run++))*F;
				out_complexR_ptr_run += 2;
				out_complexI_ptr_run += 2;
			}
			out_complexR_ptr_run -= 2 * Nx_u;
			out_complexI_ptr_run -= 2 * Nx_u;
			cx_max = Nx_u - cx_max;
			for (cx = 0; cx < cx_max; cx++)
			{
				*inter_real_ptr_run++ = ((*out_complexR_ptr_run)*(*KFC_conj_real_ptr_run) - (*out_complexI_ptr_run)*(*KFC_conj_imag_ptr_run))*F;
				*inter_imag_ptr_run++ = ((*out_complexR_ptr_run)*(*KFC_conj_imag_ptr_run++) + (*out_complexI_ptr_run)*(*KFC_conj_real_ptr_run++))*F;
				out_complexR_ptr_run += 2;
				out_complexI_ptr_run += 2;
			}
		}
		out_complexR_ptr_run = out_complexR_ptr_0;
		out_complexI_ptr_run = out_complexI_ptr_0;
		KFC_conj_real_ptr_run = KFC_conj_real_ptr_0; 
		KFC_conj_imag_ptr_run = KFC_conj_imag_ptr_0; 
		inter_real_ptr_run = inter_real_ptr_0;
		inter_imag_ptr_run = inter_imag_ptr_0;
	} // END_parallel


	// zero_initial -------------------------------------------------------
	{
		int zero_ind;
		int zero_ind_max = Nx_u_shared;
		float* out_real_ptr_run = out_real_ptr_0_shared;
		float* out_imag_ptr_run = out_imag_ptr_0_shared;
		for (zero_ind = 0; zero_ind < zero_ind_max; zero_ind++)
		{
			*out_real_ptr_run++ = 0;
			*out_imag_ptr_run++ = 0;
		}
		out_real_ptr_run = out_real_ptr_0_shared;
		out_imag_ptr_run = out_imag_ptr_0_shared;
	}
	// END_zero_initial --------------------------------------------------





	omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel shared(out_real_ptr_0_shared, out_imag_ptr_0_shared, inter_real_ptr_0_shared, inter_imag_ptr_0_shared, Nx_u_shared, nCh_shared)
	{

		// mexPrintf("This is thread number %d .\n", omp_get_thread_num()); 

		int M = Nx_u_shared;
		int N = nCh_shared;
		int j;

		float* v_ptr_0		= out_real_ptr_0_shared;
		float* v_ptr_run	= v_ptr_0;

		float* w_ptr_0		= inter_real_ptr_0_shared;
		float* w_ptr_run	= w_ptr_0;

#pragma omp for
		for (int i = 0; i < M; i++)
		{
			w_ptr_run = w_ptr_0 + i;
			v_ptr_run = v_ptr_0 + i;
			for (j = 0; j < N; j++)
			{
				*v_ptr_run += *w_ptr_run;
				w_ptr_run += M;
			}
		}

#pragma omp barrier

		v_ptr_0 = out_imag_ptr_0_shared;
		v_ptr_run = v_ptr_0;
		w_ptr_0 = inter_imag_ptr_0_shared;
		w_ptr_run = w_ptr_0;

#pragma omp for
		for (int i = 0; i < M; i++)
		{
			w_ptr_run = w_ptr_0 + i;
			v_ptr_run = v_ptr_0 + i;
			for (j = 0; j < N; j++)
			{
				*v_ptr_run += *w_ptr_run;
				w_ptr_run += M;
			}
		}

	} // END_parallel
} // END_function
