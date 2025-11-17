// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#ifndef __bmPre_fft_shanna1_omp_H__
#define __bmPre_fft_shanna1_omp_H__
 
extern void bmPre_fft_shanna1_omp(	float* in_complex_ptr_0,
									float* in_real_ptr_0,
									float* in_imag_ptr_0,
									int Nx_u, 
									int nCh,
									float* KFC_real_ptr_0,
									float* KFC_imag_ptr_0
									);

#endif // __bmPre_fft_shanna1_omp_H__