// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmPre_fft_nakatsha3_omp.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	// input arguments definition
	float* z_real_ptr_0			= (float*)mxGetPr(prhs[0]);
	float* z_imag_ptr_0			= (float*)mxGetPi(prhs[0]);

	int* N_u					= (int*)mxGetData(prhs[1]); 
	int Nx_u					= (int)N_u[0];
	int Ny_u					= (int)N_u[1];
	int Nz_u					= (int)N_u[2];
	int nCh						= (int)mxGetN(prhs[0]);

	int M						= (int)mxGetM(prhs[0]); 
	int N						= (int)mxGetN(prhs[0]);


	if (nrhs != 2)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("Input vectors must be single. ");

	if (Nx_u*Ny_u*Nz_u != M || nCh != N)
		mexErrMsgTxt("Wrong input size. ");

	// output arguments definition
	mwSize* out_size = new mwSize[2];
	out_size[0] = (mwSize)Nx_u*Ny_u*Nz_u;
	out_size[1] = (mwSize)nCh;
	mwSize out_ndims = (mwSize)2; 
	plhs[0] = mxCreateNumericArray(out_ndims, out_size, mxSINGLE_CLASS, mxCOMPLEX);
	float* out_real_ptr_run  = (float*)mxGetPr(plhs[0]);
	float* out_imag_ptr_run  = (float*)mxGetPi(plhs[0]);

	float* z_complex_ptr_0		= new float[2*nCh*Nx_u*Ny_u*Nz_u];
	float* z_complex_ptr_run	= z_complex_ptr_0; 
	
	// function call
	bmPre_fft_nakatsha3_omp(z_complex_ptr_0,
							z_real_ptr_0,
							z_imag_ptr_0,
							Nx_u,
							Ny_u,
							Nz_u,
							nCh);


	for (int i = 0; i < nCh*Nx_u*Ny_u*Nz_u; i++)
	{
		*out_real_ptr_run++ = *z_complex_ptr_run++;
		*out_imag_ptr_run++ = *z_complex_ptr_run++;
	}

	// delete[]
	delete[] out_size; 
	delete[] z_complex_ptr_0; 
}
