// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmPre_fft_shanna2_omp.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

	// input arguments definition
	float* x_real_ptr_0		= (float*)mxGetPr(prhs[0]);
	float* x_imag_ptr_0		= (float*)mxGetPi(prhs[0]);

	float* KFC_real_ptr_0	= (float*)mxGetPr(prhs[1]);
	float* KFC_imag_ptr_0	= (float*)mxGetPi(prhs[1]);

	int* N_u				= (int*)mxGetData(prhs[2]);
	int Nx_u				= (int)N_u[0]; 
	int Ny_u				= (int)N_u[1];
	int nCh					= (int)mxGetN(prhs[1]);

	int M					= (int)mxGetM(prhs[0]); 
	int N					= (int)mxGetN(prhs[1]);

	

	if (nrhs != 3)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]))
		mexErrMsgTxt("Input vectors must be single. ");

	if (!mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]))
		mexErrMsgTxt("Input vectors must be complex. ");

	if (Nx_u * Ny_u != M || nCh != N)
		mexErrMsgTxt("Wrong input size. ");




	// output arguments definition
	mwSize* y_size = new mwSize[2];
	y_size[0] = (mwSize)Nx_u*Ny_u;
	y_size[1] = (mwSize)nCh;
	mwSize y_ndims = (mwSize)2; 
	plhs[0] = mxCreateNumericArray(y_ndims, y_size, mxSINGLE_CLASS, mxCOMPLEX);
	float* y_real_ptr_0  = (float*)mxGetPr(plhs[0]);
	float* y_imag_ptr_0  = (float*)mxGetPi(plhs[0]);

	float* y_real_ptr_run = y_real_ptr_0; 
	float* y_imag_ptr_run = y_imag_ptr_0;

	float* y_complex_ptr_0 = new float[2 * nCh * Nx_u * Ny_u];
	float* y_complex_ptr_run = y_complex_ptr_0; 
	
	// function call
	bmPre_fft_shanna2_omp(	y_complex_ptr_0,
							x_real_ptr_0,
							x_imag_ptr_0,
							Nx_u,
							Ny_u,
							nCh,
							KFC_real_ptr_0,
							KFC_imag_ptr_0
							);


	for (int i = 0; i < nCh * Nx_u * Ny_u; i++)
	{
		*y_real_ptr_run++ = *y_complex_ptr_run++;
		*y_imag_ptr_run++ = *y_complex_ptr_run++;
	}

	// delete[]
	delete[] y_size; 
	delete[] y_complex_ptr_0; 
}
