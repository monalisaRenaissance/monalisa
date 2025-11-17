// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#ifndef __bmPost_fft_nakatsha2_omp_H__
#define __bmPost_fft_nakatsha2_omp_H__
 
extern void bmPost_fft_nakatsha2_omp(	float* out_complex_ptr_0,
										float* out_real_ptr_0,
										float* out_imag_ptr_0,
										float* inter_real_ptr_0,
										float* inter_imag_ptr_0,
										float* KFC_conj_real_ptr_0,
										float* KFC_conj_imag_ptr_0,
										int Nx_u,
										int Ny_u,
										int nCh);

#endif // __bmPost_fft_nakatsha2_omp_H__