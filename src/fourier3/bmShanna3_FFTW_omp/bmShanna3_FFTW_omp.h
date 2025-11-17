// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#ifndef __bmShanna3_FFTW_omp_H__
#define __bmShanna3_FFTW_omp_H__
 
extern void bmShanna3_FFTW_omp(	float* y_real_ptr_0,
								float* y_imag_ptr_0,
								float* x_real_ptr_0,
								float* x_imag_ptr_0,
								float* KFC_real_ptr_0,
								float* KFC_imag_ptr_0,
								int Nx_u,
								int Ny_u,
								int Nz_u,
								int nCh,
								int r_size,
								int* r_jump_ptr_0,
								int* r_nJump_ptr_0,
								float* m_val_ptr_0,
								int l_size,
								int* l_jump_ptr_0,
								int l_nJump
								);

#endif // __bmShanna3_FFTW_omp_H__