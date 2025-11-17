// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmShanna1_FFTW_omp.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

	// input arguments definition
	float* x_real_ptr_0		= (float*)mxGetData(prhs[0]);
	float* x_imag_ptr_0		= (float*)mxGetData(prhs[1]);

	float* KFC_real_ptr_0	= (float*)mxGetData(prhs[2]);
	float* KFC_imag_ptr_0	= (float*)mxGetData(prhs[3]);

	int Nx_u				= (int)mxGetScalar(prhs[4]);
	int nCh					= (int)mxGetN(prhs[2]);

	int M					= (int)mxGetM(prhs[2]); 
	int N					= (int)mxGetN(prhs[2]);

	int r_size				= (int)mxGetScalar(prhs[5]);
	int* r_jump_ptr_0		= (int*)mxGetData(prhs[6]);
	int* r_nJump_ptr_0		= (int*)mxGetData(prhs[7]);

	float* m_val_ptr_0		= (float*)mxGetData(prhs[8]);

	int l_size				= (int)mxGetScalar(prhs[9]);
	int* l_jump_ptr_0		= (int*)mxGetData(prhs[10]);
	int l_nJump				= (int)mxGetScalar(prhs[11]);





	if (nrhs != 12)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 2)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]))
		mexErrMsgTxt("Input vectors must be single. ");

	if (Nx_u != M || nCh != N)
		mexErrMsgTxt("Wrong input size. ");




	// output arguments definition
	mwSize* y_size = new mwSize[2];
	y_size[0] = (mwSize)l_size;
	y_size[1] = (mwSize)nCh;
	mwSize y_ndims = (mwSize)2; 
	plhs[0] = mxCreateNumericArray(y_ndims, y_size, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(y_ndims, y_size, mxSINGLE_CLASS, mxREAL);
	float* y_real_ptr_0  = (float*)mxGetData(plhs[0]);
	float* y_imag_ptr_0  = (float*)mxGetData(plhs[1]);

	
	// function call
	bmShanna1_FFTW_omp( y_real_ptr_0, 
						y_imag_ptr_0, 
						x_real_ptr_0, 
						x_imag_ptr_0,
						KFC_real_ptr_0, 
						KFC_imag_ptr_0, 
						Nx_u, 
						nCh, 
						r_size, 
						r_jump_ptr_0, 
						r_nJump_ptr_0,
						m_val_ptr_0, 
						l_size, 
						l_jump_ptr_0, 
						l_nJump
		);

	// delete[]
	delete[] y_size; 
}
