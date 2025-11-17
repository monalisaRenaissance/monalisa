// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include <cmath>
#include <omp.h>
#include "bmPre_fft_nakatsha3_omp.h"
// #include "mex.h"

void bmPre_fft_nakatsha3_omp(	float* in_complex_ptr_0_shared, 
								float* in_real_ptr_0_shared,
								float* in_imag_ptr_0_shared,
								int Nx_u_shared,
								int Ny_u_shared,
								int Nz_u_shared,
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


#pragma omp parallel shared(in_complex_ptr_0_shared, in_real_ptr_0_shared, in_imag_ptr_0_shared, Nx_u_shared, Ny_u_shared, Nz_u_shared, nCh_shared)
	{

		// mexPrintf("This is thread number %d .\n", omp_get_thread_num()); 

		float* in_real_ptr_0 = in_real_ptr_0_shared;
		float* in_real_ptr_run = in_real_ptr_0;

		float* in_imag_ptr_0 = in_imag_ptr_0_shared;
		float* in_imag_ptr_run = in_imag_ptr_0;

		float* in_complex_ptr_0 = in_complex_ptr_0_shared;
		float* in_complex_ptr_run = in_complex_ptr_0;

		int Nx_u = Nx_u_shared; 
		double Nx_u_double = (double)Nx_u;
		int cx_max;
		int cx;
		int Nx_half = (int)floor(Nx_u_double / 2);

		int Ny_u = Ny_u_shared;
		double Ny_u_double = (double)Ny_u;
		int cy_max;
		int cy;
		int Ny_half = (int)floor(Ny_u_double / 2);

		int Nz_u = Nz_u_shared;
		double Nz_u_double = (double)Nz_u;
		int cz_max;
		int cz;
		int Nz_half = (int)floor(Nz_u_double / 2);

		int nCh = nCh_shared; 

		long long start_incr = ((long long)Nx_half) + ((long long)Nx_u)*((long long)Ny_half) + ((long long)Nx_u)*((long long)Ny_u)*((long long)Nz_half);

#pragma omp for
		for (int ind_ch = 0; ind_ch < nCh; ind_ch++)
		{
			in_real_ptr_run = in_real_ptr_0 + start_incr + ((long long)ind_ch)*((long long)Nx_u)*((long long)Ny_u)*((long long)Nz_u);
			in_imag_ptr_run = in_imag_ptr_0 + start_incr + ((long long)ind_ch)*((long long)Nx_u)*((long long)Ny_u)*((long long)Nz_u);
			in_complex_ptr_run = in_complex_ptr_0 + ((long long)2)*((long long)ind_ch)*((long long)Nx_u)*((long long)Ny_u)*((long long)Nz_u);

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
		}
		in_real_ptr_run = in_real_ptr_0;
		in_imag_ptr_run = in_imag_ptr_0;
		in_complex_ptr_run = in_complex_ptr_0;
	} // END_parallel
} // END_function
