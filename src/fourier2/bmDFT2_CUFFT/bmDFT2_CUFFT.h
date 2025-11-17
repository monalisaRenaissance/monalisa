// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#ifndef __bmDFT2_CUFFT_H__
#define __bmDFT2_CUFFT_H__
 
extern void bmDFT2_CUFFT(float* out_real_ptr_0, float* out_imag_ptr_0, float* in_real_ptr_0, float* in_imag_ptr_0, int Nx_u, int Ny_u,  float dKx_u, float dKy_u, int nCh);

#endif // __bmDFT2_CUFFT_H__