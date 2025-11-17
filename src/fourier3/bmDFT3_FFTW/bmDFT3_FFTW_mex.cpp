// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmDFT3_FFTW.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// initial ----------------------------------------------------------------

	if (nrhs != 4)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 2)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("Input vectors must be single. ");

	int m			= (int)mxGetM(prhs[0]);
	int n			= (int)mxGetN(prhs[0]);

	float* in_real_ptr = (float*)mxGetData(prhs[0]);
	float* in_imag_ptr = (float*)mxGetData(prhs[1]);

	int* N_u = (int*)mxGetPr(prhs[2]);
	int Nx_u = N_u[0];
	int Ny_u = N_u[1];
	int Nz_u = N_u[2];
	int nCh = n / Ny_u / Nz_u; 

	float* dK_u = (float*)mxGetPr(prhs[3]);
	float dKx_u = dK_u[0];
	float dKy_u = dK_u[1];
	float dKz_u = dK_u[2];


	mwSize* out_size = new mwSize[4];
	out_size[0] = (mwSize)Nx_u;
	out_size[1] = (mwSize)Ny_u;
	out_size[2] = (mwSize)Nz_u;
	out_size[3] = (mwSize)nCh;
	mwSize out_ndims = (mwSize)4;
	plhs[0] = mxCreateNumericArray(out_ndims, out_size, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(out_ndims, out_size, mxSINGLE_CLASS, mxREAL);
	float* out_real_ptr = (float*)mxGetData(plhs[0]);
	float* out_imag_ptr = (float*)mxGetData(plhs[1]);


	// END_initial ------------------------------------------------------------
	
	// FFT call
	bmDFT3_FFTW(out_real_ptr, out_imag_ptr, in_real_ptr, in_imag_ptr, Nx_u, Ny_u, Nz_u, dKx_u, dKy_u, dKz_u, nCh);

	delete[] out_size;
}


