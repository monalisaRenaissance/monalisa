// Bastien Milani
// CHUV and UNIL
// Lausanne - Switzerland
// May 2023

#include "mex.h"
#include "bmPost_fft_nakatsha2_omp.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

	// input arguments definition
	float* z_real_ptr_0			= (float*)mxGetPr(prhs[0]);
	float* z_imag_ptr_0			= (float*)mxGetPi(prhs[0]);

	float* KFC_conj_real_ptr_0	= (float*)mxGetPr(prhs[1]);
	float* KFC_conj_imag_ptr_0	= (float*)mxGetPi(prhs[1]);

	int* N_u					= (int*)mxGetData(prhs[2]);
	int Nx_u					= (int)N_u[0];
	int Ny_u					= (int)N_u[1];
	int nCh						= (int)mxGetN(prhs[1]);

	int M						= (int)mxGetM(prhs[1]); 
	int N						= (int)mxGetN(prhs[1]);


	if (nrhs != 3)
		mexErrMsgTxt("Invaid number of inputs. ");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of outputs. ");

	if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]))
		mexErrMsgTxt("Input vectors must be single. ");

	if (Nx_u*Ny_u != M || nCh != N)
		mexErrMsgTxt("Wrong input size. ");


	// output arguments definition
	mwSize* out_size = new mwSize[2];
	out_size[0] = (mwSize)(Nx_u * Ny_u);
	out_size[1] = (mwSize)1;
	mwSize out_ndims = (mwSize)2; 
	plhs[0] = mxCreateNumericArray(out_ndims, out_size, mxSINGLE_CLASS, mxCOMPLEX);
	float* out_real_ptr_0  = (float*)mxGetPr(plhs[0]);
	float* out_imag_ptr_0  = (float*)mxGetPi(plhs[0]);



	float* z_real_ptr_run = z_real_ptr_0; 
	float* z_imag_ptr_run = z_imag_ptr_0;
	float* z_complex_ptr_0   = new float[2 * nCh * Nx_u *Ny_u];
	float* z_complex_ptr_run = z_complex_ptr_0; 
	float* inter_complex_ptr_0 = new float[2 * nCh * Nx_u * Ny_u];
	float* inter_real_ptr_0 = inter_complex_ptr_0;
	float* inter_imag_ptr_0 = inter_complex_ptr_0 + nCh * Nx_u * Ny_u;
	

	for (int i = 0; i < nCh * Nx_u * Ny_u; i++)
	{
		*z_complex_ptr_run++ = *z_real_ptr_run++; 
		*z_complex_ptr_run++ = *z_imag_ptr_run++;
	}
	z_complex_ptr_run	= z_complex_ptr_0; 
	z_real_ptr_run		= z_real_ptr_0;
	z_imag_ptr_run		= z_imag_ptr_0;
	

	// function call
	bmPost_fft_nakatsha2_omp(	z_complex_ptr_0,
								out_real_ptr_0,
								out_imag_ptr_0,
								inter_real_ptr_0,
								inter_imag_ptr_0,
								KFC_conj_real_ptr_0,
								KFC_conj_imag_ptr_0,
								Nx_u,
								Ny_u, 
								nCh); 

	// delete[]
	delete[] out_size; 
}
