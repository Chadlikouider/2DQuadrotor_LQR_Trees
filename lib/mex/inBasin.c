/*=========================================================
 * inBasin.c - rapid evaluation of whether sample is in tree
 * returns closest node in tree in any case
 *
 * [inBasin, closest] = inBasin(xsamp,x0,S,rho,nElem)
 *
 **INPUT***************************************************
 * xsamp :      sample state         (nx x 1) 
 * x0 :         nominal states       (nx x nn) 
 * S :          cost-to-go matrices  (nx x nx*nn) 
 * rho :        funnel parameters    (1 x nn) 
 * nElem :      number elem. in arr. (1)
 *
 **OUTPUT**************************************************
 * inBasin:     bool that indicates if sample is in basin of tree
 * closest :    index of closest node that has minimal
 *              (xsamp - xn)^T * Sn *(xsamp-xn)
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
    double *xsamp, *x0, *S, *rho, *inBasin, *closest, *vec, *dx;
    mxArray *vecA, *dxA;
    mwSize nx,nx2,oneI = 1, i = 0, j = 0, nElem = 0;      /* matrix dimensions, nx = nstates */
    /* form of op(A) & op(B) to use in matrix multiplication */
    char *chn = "N"; // do A*b
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0, value = 0.0, minVal = 0.0;

    xsamp = mxGetPr(prhs[0]); /* sample state (nx x 1) */
    x0 = mxGetPr(prhs[1]); /* nominal states (nx x nn) */
    S  = mxGetPr(prhs[2]); /* weight matrices (nx x nx*nn) */
    rho = mxGetPr(prhs[3]); /* ellipsis "radii" (1 x nn) */
    nElem = (mwSize) mxGetScalar(prhs[4]); /* number of elements to analyze (1 x nn) */
   
    /* dimensions of input matrices */
    nx = (mwSignedIndex)mxGetM(prhs[0]); // dimension of the state
    nx2 = nx*nx; // precompute number of elements in weight matrices

    /* create output scalars ind, indE */
    plhs[0] = mxCreateDoubleScalar(0.0);
    plhs[1] = mxCreateDoubleScalar(-1.0);

    // create dummy vector
    vecA = mxCreateDoubleMatrix(nx, 1, mxREAL);
    dxA = mxCreateDoubleMatrix(nx, 1, mxREAL);
    vec = mxGetPr(vecA);
    dx = mxGetPr(dxA);
    inBasin = mxGetPr(plhs[0]); /* bool stating if sample is in tree basin */
    closest = mxGetPr(plhs[1]); /* index of closest node (scalar) */
 
    // for all nodes, calculate cost to go, check if rho - cost is larger
    // than zero and update best node if smaller than minimum
    for(i = 1; i < nElem; ++i) {
        for(j = 0; j < nx; ++j){
            *(dx + j) = *(xsamp + j) - *(x0 + i*nx + j);
        }
        // calculate weighted vector norm
        dgemv(chn, &nx, &nx, &one, (S + nx2*i), &nx, dx, &oneI, &zero, vec, &oneI);
        // get dot product of resulting vector with difference vector
        value = ddot(&nx,dx,&oneI,vec,&oneI);
        // check if sample is in the node's funnel:
        if( *(rho + i) > value) {
            *inBasin = 1;
            return;
        }
        // keep track of the overall minimum:
        if(*closest < 0 || minVal > value) {
            minVal = value;
            *closest = (double) (i + 1);
        }
    }
}
