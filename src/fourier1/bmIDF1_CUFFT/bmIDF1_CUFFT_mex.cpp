// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmIDF1_CUFFT.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// initial ----------------------------------------------------------------

	if (nrhs != 4)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 2)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("Input vectors must be single. ");

	int m	= (int)mxGetM(prhs[0]);
	int nCh	= (int)mxGetN(prhs[0]);


	float* in_real_ptr = (float*)mxGetData(prhs[0]);
	float* in_imag_ptr = (float*)mxGetData(prhs[1]);

	int Nx_u	= (int)mxGetScalar(prhs[2]); 
	float dKx_u = (float)mxGetScalar(prhs[3]);

	if (Nx_u != m)
		mexErrMsgTxt("Invalid input size. ");

	plhs[0] = mxCreateNumericMatrix(Nx_u, nCh, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(Nx_u, nCh, mxSINGLE_CLASS, mxREAL);
	float* out_real_ptr = (float*)mxGetData(plhs[0]);
	float* out_imag_ptr = (float*)mxGetData(plhs[1]);

	// END_initial ------------------------------------------------------------
	
	// FFT call
	bmIDF1_CUFFT(out_real_ptr, out_imag_ptr, in_real_ptr, in_imag_ptr, Nx_u, dKx_u, nCh);

}


