// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include <cmath>
#include <omp.h>
#include "bmPre_fft_shanna2_omp.h"
// #include "mex.h"

void bmPre_fft_shanna2_omp(	float* in_complex_ptr_0_shared, 
							float* in_real_ptr_0_shared,
							float* in_imag_ptr_0_shared,
							int Nx_u_shared,
							int Ny_u_shared, 
							int nCh_shared,
							float* KFC_real_ptr_0_shared,
							float* KFC_imag_ptr_0_shared)
{

	if (nCh_shared < omp_get_max_threads())
	{
		omp_set_num_threads(nCh_shared);
	}
	else
	{
		omp_set_num_threads(omp_get_max_threads());
	}

#pragma omp parallel shared(in_complex_ptr_0_shared, in_real_ptr_0_shared, in_imag_ptr_0_shared, Nx_u_shared, Ny_u_shared, nCh_shared, KFC_real_ptr_0_shared, KFC_imag_ptr_0_shared)
	{

		// mexPrintf("This is thread number %d .\n", omp_get_thread_num()); 


		float* in_real_ptr_0	= in_real_ptr_0_shared; 
		float* in_real_ptr_run	= in_real_ptr_0;

		float* in_imag_ptr_0	= in_imag_ptr_0_shared;
		float* in_imag_ptr_run	= in_imag_ptr_0; 

		float* in_complex_ptr_0 = in_complex_ptr_0_shared;
		float* in_complex_ptr_run = in_complex_ptr_0;
		
		float* KFC_real_ptr_0 = KFC_real_ptr_0_shared;
		float* KFC_real_ptr_run = KFC_real_ptr_0;
		
		float* KFC_imag_ptr_0 = KFC_imag_ptr_0_shared;
		float* KFC_imag_ptr_run = KFC_imag_ptr_0;



		int Nx_u = Nx_u_shared; 
		int Ny_u = Ny_u_shared;

		double Nx_u_double = (double)Nx_u;
		double Ny_u_double = (double)Ny_u;

		int cx_max;
		int cx;
		int cy_max;
		int cy;

		int Nx_half = (int)floor(Nx_u_double / 2);
		int Ny_half = (int)floor(Ny_u_double / 2);
		
		int nCh = nCh_shared; 

		long long start_incr = ((long long)Nx_half) + ((long long)Nx_u)*((long long)Ny_half); 

		in_real_ptr_run += start_incr;
		in_imag_ptr_run += start_incr;


#pragma omp for
		for (int ind_ch = 0; ind_ch < nCh; ind_ch++)
		{

			KFC_real_ptr_run = KFC_real_ptr_0 + start_incr + ((long long)ind_ch)*((long long)Nx_u)*((long long)Ny_u);
			KFC_imag_ptr_run = KFC_imag_ptr_0 + start_incr + ((long long)ind_ch)*((long long)Nx_u)*((long long)Ny_u);
			in_complex_ptr_run = in_complex_ptr_0 + ((long long)2)*((long long)ind_ch)*((long long)Nx_u)*((long long)Ny_u);


			cy_max = Ny_u - Ny_half;
			for (cy = 0; cy < cy_max; cy++)
			{
				cx_max = Nx_u - Nx_half;
				for (cx = 0; cx < cx_max; cx++)
				{
					*in_complex_ptr_run++ = (*in_real_ptr_run)*(*KFC_real_ptr_run) - (*in_imag_ptr_run)*(*KFC_imag_ptr_run);
					*in_complex_ptr_run++ = (*in_real_ptr_run++)*(*KFC_imag_ptr_run++) + (*in_imag_ptr_run++)*(*KFC_real_ptr_run++);
				}

				in_real_ptr_run -= Nx_u;
				in_imag_ptr_run -= Nx_u;
				KFC_real_ptr_run -= Nx_u;
				KFC_imag_ptr_run -= Nx_u;
				cx_max = Nx_u - cx_max;
				for (cx = 0; cx < cx_max; cx++)
				{
					*in_complex_ptr_run++ = (*in_real_ptr_run)*(*KFC_real_ptr_run) - (*in_imag_ptr_run)*(*KFC_imag_ptr_run);
					*in_complex_ptr_run++ = (*in_real_ptr_run++)*(*KFC_imag_ptr_run++) + (*in_imag_ptr_run++)*(*KFC_real_ptr_run++);
				}
				in_real_ptr_run += Nx_u;
				in_imag_ptr_run += Nx_u;
				KFC_real_ptr_run += Nx_u;
				KFC_imag_ptr_run += Nx_u;
			}
			in_real_ptr_run  -= Nx_u*Ny_u;
			in_imag_ptr_run  -= Nx_u*Ny_u;
			KFC_real_ptr_run -= Nx_u*Ny_u;
			KFC_imag_ptr_run -= Nx_u*Ny_u;
			cy_max = Ny_u - cy_max;
			for (cy = 0; cy < cy_max; cy++)
			{
				cx_max = Nx_u - Nx_half;
				for (cx = 0; cx < cx_max; cx++)
				{
					*in_complex_ptr_run++ = (*in_real_ptr_run)*(*KFC_real_ptr_run) - (*in_imag_ptr_run)*(*KFC_imag_ptr_run);
					*in_complex_ptr_run++ = (*in_real_ptr_run++)*(*KFC_imag_ptr_run++) + (*in_imag_ptr_run++)*(*KFC_real_ptr_run++);
				}

				in_real_ptr_run -= Nx_u;
				in_imag_ptr_run -= Nx_u;
				KFC_real_ptr_run -= Nx_u;
				KFC_imag_ptr_run -= Nx_u;
				cx_max = Nx_u - cx_max;
				for (cx = 0; cx < cx_max; cx++)
				{
					*in_complex_ptr_run++ = (*in_real_ptr_run)*(*KFC_real_ptr_run) - (*in_imag_ptr_run)*(*KFC_imag_ptr_run);
					*in_complex_ptr_run++ = (*in_real_ptr_run++)*(*KFC_imag_ptr_run++) + (*in_imag_ptr_run++)*(*KFC_real_ptr_run++);
				}
				in_real_ptr_run += Nx_u;
				in_imag_ptr_run += Nx_u;
				KFC_real_ptr_run += Nx_u;
				KFC_imag_ptr_run += Nx_u;
			}
		} // END_par_for
		in_real_ptr_run = in_real_ptr_0;
		in_imag_ptr_run = in_imag_ptr_0;
		KFC_real_ptr_run = KFC_real_ptr_0;
		KFC_imag_ptr_run = KFC_imag_ptr_0;
		in_complex_ptr_run = in_complex_ptr_0;

	} // END_parallel
}// END_function

