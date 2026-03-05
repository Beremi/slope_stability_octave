/*
 * constitutive_problem_2D_SDS_mex.c
 *
 * Mex wrapper for stress + consistent-tangent evaluation of the 2D
 * Mohr-Coulomb constitutive model. OpenMP-parallel over integration points.
 *
 * Usage:
 *   [S, DS] = constitutive_problem_2D_SDS_mex(E, c_bar, sin_phi,
 *                                              shear, bulk, lame)
 *
 *   E        - 3 x n_int
 *   c_bar    - 1 x n_int
 *   sin_phi  - 1 x n_int
 *   shear    - 1 x n_int
 *   bulk     - 1 x n_int
 *   lame     - 1 x n_int
 *
 *   S        - 4 x n_int
 *   DS       - 9 x n_int   (column-major 3x3 blocks)
 */

#include "mex.h"
#include "constitutive_2D_kernel.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
        mexErrMsgTxt("constitutive_problem_2D_SDS_mex: "
                     "6 inputs required (E, c_bar, sin_phi, shear, bulk, lame)");
    if (nlhs > 2)
        mexErrMsgTxt("constitutive_problem_2D_SDS_mex: up to 2 outputs (S, DS)");

    const double *E = mxGetPr(prhs[0]);
    const double *c_bar = mxGetPr(prhs[1]);
    const double *sin_phi = mxGetPr(prhs[2]);
    const double *shear = mxGetPr(prhs[3]);
    const double *bulk = mxGetPr(prhs[4]);
    const double *lam = mxGetPr(prhs[5]);

    if (mxGetM(prhs[0]) != 3)
        mexErrMsgTxt("constitutive_problem_2D_SDS_mex: E must have 3 rows.");

    mwSize n_int = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(4, n_int, mxREAL);
    double *S = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(9, n_int, mxREAL);
    double *DS = mxGetPr(plhs[1]);

#pragma omp parallel for schedule(static)
    for (mwSize p = 0; p < n_int; p++)
    {
        constitutive_2D_point(E + 3 * p,
                              c_bar[p], sin_phi[p],
                              shear[p], bulk[p], lam[p],
                              S + 4 * p, DS + 9 * p);
    }
}
