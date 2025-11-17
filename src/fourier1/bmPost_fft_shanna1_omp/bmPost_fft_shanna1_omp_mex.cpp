// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmPost_fft_shanna1_omp.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

	// input arguments definition
	float* x_real_ptr_0		= (float*)mxGetPr(prhs[0]);
	float* x_imag_ptr_0		= (float*)mxGetPi(prhs[0]);

	int Nx_u				= (int)mxGetScalar(prhs[1]);
	int nCh					= (int)mxGetN(prhs[0]);

	int M					= (int)mxGetM(prhs[0]); 
	int N					= (int)mxGetN(prhs[0]);

	if (nrhs != 2)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("Input vectors must be single. ");

	if (Nx_u != M || nCh != N)
		mexErrMsgTxt("Wrong input size. ");




	// output arguments definition
	mwSize* y_size = new mwSize[2];
	y_size[0] = (mwSize)Nx_u;
	y_size[1] = (mwSize)nCh;
	mwSize y_ndims = (mwSize)2; 
	plhs[0] = mxCreateNumericArray(y_ndims, y_size, mxSINGLE_CLASS, mxCOMPLEX);
	float* y_real_ptr_0  = (float*)mxGetPr(plhs[0]);
	float* y_imag_ptr_0  = (float*)mxGetPi(plhs[0]);


	float* x_complex_ptr_0 = new float[2 * Nx_u*nCh];
	float* x_complex_ptr_run = x_complex_ptr_0;
	
	float* x_real_ptr_run = x_real_ptr_0; 
	float* x_imag_ptr_run = x_imag_ptr_0;

	for (int i = 0; i < Nx_u * nCh; i++)
	{
		*x_complex_ptr_run++ = *x_real_ptr_run++;
		*x_complex_ptr_run++ = *x_imag_ptr_run++;
	}
	x_complex_ptr_run = x_complex_ptr_0; 
	x_real_ptr_run = x_real_ptr_0; 
	x_imag_ptr_run = x_imag_ptr_0; 

	// function call
	bmPost_fft_shanna1_omp(	x_complex_ptr_0,
							y_real_ptr_0,
							y_imag_ptr_0,
							Nx_u,
							nCh); 


	// delete[]
	delete[] y_size; 
	delete[] x_complex_ptr_0;

}
