/*=========================================================
 * inBasin.c - rapid evaluation of whether sample is in tree
 * returns closest node in tree in any case
 *
 * [newRho] = inBasin(xsim, x0, S, rho, nElem)
 *
 **INPUT***************************************************
 * xsim :       failed traj. from sim.                  (nx x ntraj) 
 * x0 :         nominal states of trajectory followed   (nx x ntraj) 
 * S :          cost-to-go matrices of traj. followed   (nx x nx*ntraj) 
 * rho :        funnel parameters                       (1 x ntraj)
 *
 **OUTPUT**************************************************
 * newRho:     updated funnel parameters    (1 x ntraj)
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
    double *xtape, *x0, *S, *rho, *newRho, *vec, *dx;
    double value = 0;
    mxArray *vecA, *dxA;
    mwSize nx,nx2,oneI = 1, i = 0, j = 0, nElem; /* matrix dimensions, nx = nstates */
    /* form of op(A) & op(B) to use in matrix multiplication */
    char *chn = "N"; // do A*b
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0;

    xtape = mxGetPr(prhs[0]);   /* sample state       (nx x 1) */
    x0 = mxGetPr(prhs[1]);      /* nominal states     (nx x nn) */
    S  = mxGetPr(prhs[2]);      /* weight matrices    (nx x nx*nn) */
    rho = mxGetPr(prhs[3]);     /* funnel parameters  (1 x nn) */
    
    nElem = (mwSize)mxGetN(prhs[3]); /* number of elements to analyze (1) */
   
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
    
    // init pointer to output array
    newRho = mxGetPr(plhs[0]); /* cost to go */
        
    // for all nodes of the nominal trajectory, update funnel parameter if
    // necessary
    for(i = 0; i < nElem; ++i) {
        for(j = 0; j < nx; ++j){
            *(dx + j) = *(xtape + i*nx + j) - *(x0 + i*nx + j);
        }
        // calculate weighted vector norm
        dgemv(chn, &nx, &nx, &one, (S + nx2*i), &nx, dx, &oneI, &zero, vec, &oneI);
        // get dot product of resulting vector with difference vector
        value = ddot(&nx,dx,&oneI,vec,&oneI);
        // only update the funnel parameter if the state along the simulated,
        // failed trajectory is outside the funnel:
        if(*(rho + i) > value){
            *(newRho + i) = value;
        }else{
            *(newRho + i) = -1;
        }
    }
}
