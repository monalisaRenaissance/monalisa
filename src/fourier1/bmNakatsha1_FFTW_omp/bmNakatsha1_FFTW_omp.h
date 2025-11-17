// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#ifndef __bmNakatsha1_FFTW_omp_H__
#define __bmNakatsha1_FFTW_omp_H__
 
extern void bmNakatsha1_FFTW_omp(float* x_real_ptr_0,
								 float* x_imag_ptr_0,
								 float* y_real_ptr_0,
								 float* y_imag_ptr_0,
								 float* KFC_conj_real_ptr_0,
								 float* KFC_conj_imag_ptr_0,
								 int Nx_u,
								 int nCh,
								 int r_size,
								 int* r_jump_ptr_0,
								 int* r_nJump_ptr_0,
								 float* m_val_ptr_0,
								 int l_size,
								 int* l_jump_ptr_0,
								 int l_nJump
								 );

#endif // __bmNakatsha1_FFTW_omp_H__