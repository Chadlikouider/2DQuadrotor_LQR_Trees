/*=========================================================
 * getC2G.c - rapid evaluation of cost-to-go for all elements
 * in Tree given a sample state
 *
 * ji = getC2G(xsamp,x0,S,nElem)
 *
 **INPUT***************************************************
 * xsamp :      sample state         (nx x 1) 
 * x0 :         nominal states       (nx x nn) 
 * S :          cost-to-go matrices  (nx x nx*nn)
 * nElem :      number elem. in arr. (1)
 *
 **OUTPUT**************************************************
 * ji:     cost-to-go for each node i:
 *         j(i) = (xsamp-x0(:,i))'*S(:,((i-1)*nx + 1):i*nx)*(xsamp-x0(:,i));
 *
 * compile with
 * mex('-v', '-largeArrayDims', blaslib)
 * where blaslib contains string pointing to libmwblas.lib
 * of the compiler you set in mex -setup
 * for example on a 64-bit Win7 machine:
 * blaslib = [matlabroot, '\extern\lib\win64\microsoft\libmwblas.lib'];
 *=======================================================*/

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* pointers to input matrices and output scalars*/
    double *xsamp, *x0, *S, *ji, *vec, *dx;
    mxArray *vecA, *dxA;
    mwSize nx,nx2,oneI = 1, i = 0, j = 0, nElem;      /* matrix dimensions, nx = nstates */
    /* form of op(A) & op(B) to use in matrix multiplication */
    char *chn = "N"; // do A*b
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0;

    xsamp = mxGetPr(prhs[0]); /* sample state (nx x 1) */
    x0 = mxGetPr(prhs[1]); /* nominal states (nx x nn) */
    S  = mxGetPr(prhs[2]); /* weight matrices (nx x nx*nn) */
    nElem = (mwSize) mxGetScalar(prhs[3]); /* number of elements to analyze (1) */
   
    /* dimensions of input matrices */
    nx = (mwSize)mxGetM(prhs[0]); // dimension of the state
    nx2 = nx*nx; // precompute number of elements in weight matrices

    /* create output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, nElem, mxREAL);

    // create dummy vector
    vecA = mxCreateDoubleMatrix(nx, 1, mxREAL);
    dxA = mxCreateDoubleMatrix(nx, 1, mxREAL);
    vec = mxGetPr(vecA);
    dx = mxGetPr(dxA);
    ji = mxGetPr(plhs[0]); /* cost to go */
        
    // for all nodes, calculate cost to go
    for(i = 0; i < nElem; ++i) {
        for(j = 0; j < nx; ++j){
            *(dx + j) = *(xsamp + j) - *(x0 + i*nx + j);
        }
        // calculate weighted vector norm
        dgemv(chn, &nx, &nx, &one, (S + nx2*i), &nx, dx, &oneI, &zero, vec, &oneI);
        // get dot product of resulting vector with difference vector
        *(ji + i) = ddot(&nx,dx,&oneI,vec,&oneI);
    }
}
