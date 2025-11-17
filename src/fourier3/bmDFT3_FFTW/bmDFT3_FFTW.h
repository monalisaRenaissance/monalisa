// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#ifndef __bmDFT3_FFTW_H__
#define __bmDFT3_FFTW_H__
 
extern void bmDFT3_FFTW(float* out_real_ptr_0, float* out_imag_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u, int Nz_u, float dKx_u, float dKy_u, float dKz_u, int nCh);

#endif // __bmDFT3_FFTW_H__