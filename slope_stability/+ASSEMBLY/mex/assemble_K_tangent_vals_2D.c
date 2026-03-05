/*
 * assemble_K_tangent_vals_2D.c
 *
 * Element-level tangent stiffness assembly for 2D FEM.
 *
 * Computes V_tang = values of K_tangent(Q,Q) at the precomputed sparsity
 * pattern positions, using element-local quadrature:
 *
 *   K_tangent = sum_e sum_q  B_eq' * (w_q * D_eq) * B_eq
 *
 * where B_eq is reconstructed from DPhi1/2 at each integration point,
 * D_eq = reshape(DS(:,g), 3, 3), and w_q = WEIGHT(g).
 *
 * Usage:
 *   V_tang = assemble_K_tangent_vals_2D(DPhi1, DPhi2, DS, WEIGHT, ...
 *                                        scatter_map, n_q, nnz_out);
 */

#include "mex.h"
#include <string.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Fixed dimensions for 2D: n_strain=3, dim=2 */
#define NS 3
#define DIM 2

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 7)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:nrhs",
                          "Seven inputs required: DPhi1, DPhi2, DS, WEIGHT, scatter_map, n_q, nnz_out");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:nlhs",
                          "At most one output.");

    const double *DPhi1 = mxGetPr(prhs[0]);  /* n_p x n_int */
    const double *DPhi2 = mxGetPr(prhs[1]);  /* n_p x n_int */
    const double *DS = mxGetPr(prhs[2]);     /* 9 x n_int */
    const double *WEIGHT = mxGetPr(prhs[3]); /* 1 x n_int */

    if (!mxIsInt64(prhs[4]))
        mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:type",
                          "scatter_map must be int64.");
    const int64_t *scatter_map = (const int64_t *)mxGetData(prhs[4]);

    const int n_q = (int)mxGetScalar(prhs[5]);
    const int nnz_out = (int)mxGetScalar(prhs[6]);

    const mwSize n_p = mxGetM(prhs[0]);
    const mwSize n_int = mxGetN(prhs[0]);
    const mwSize n_e = mxGetN(prhs[4]);
    const int n_local_dof = (int)(DIM * n_p);
    const int n_ld2 = n_local_dof * n_local_dof;

    if (mxGetM(prhs[1]) != n_p || mxGetN(prhs[1]) != n_int)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:size",
                          "DPhi1 and DPhi2 must have the same size.");
    if (mxGetM(prhs[2]) != NS * NS || mxGetN(prhs[2]) != n_int)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:size",
                          "DS must have size (9, n_int) for 2D.");
    if (mxGetM(prhs[4]) != (mwSize)n_ld2)
        mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:size",
                          "scatter_map row count must be n_local_dof^2.");

    plhs[0] = mxCreateDoubleMatrix((mwSize)nnz_out, 1, mxREAL);
    double *V_tang = mxGetPr(plhs[0]);

#pragma omp parallel
    {
        double *Ke = (double *)calloc((size_t)n_ld2, sizeof(double));
        double *B_eq = (double *)calloc((size_t)NS * (size_t)n_local_dof, sizeof(double));
        double *D_eq = (double *)calloc((size_t)NS * (size_t)NS, sizeof(double));
        double *tmp = (double *)calloc((size_t)NS * (size_t)n_local_dof, sizeof(double));

        if (!Ke || !B_eq || !D_eq || !tmp) {
#pragma omp critical
            {
                mexErrMsgIdAndTxt("assemble_K_tangent_vals_2D:alloc",
                                  "Failed to allocate thread-local buffers.");
            }
        }

#pragma omp for schedule(dynamic, 256)
        for (mwSize e = 0; e < n_e; e++)
        {
            memset(Ke, 0, sizeof(double) * (size_t)n_ld2);
            const mwSize g_base = e * (mwSize)n_q;

            for (int q = 0; q < n_q; q++)
            {
                const mwSize g = g_base + (mwSize)q;
                const double w = WEIGHT[g];

                const double *ds_col = DS + g * NS * NS;
                for (int k = 0; k < NS * NS; k++) {
                    D_eq[k] = w * ds_col[k];
                }

                memset(B_eq, 0, sizeof(double) * (size_t)NS * (size_t)n_local_dof);
                const double *dp1 = DPhi1 + g * n_p;
                const double *dp2 = DPhi2 + g * n_p;

                for (mwSize i = 0; i < n_p; i++)
                {
                    const double dN1 = dp1[i];
                    const double dN2 = dp2[i];
                    const int c = (int)i * DIM;

                    /* eps_11 */
                    B_eq[0 + (c + 0) * NS] = dN1;
                    /* eps_22 */
                    B_eq[1 + (c + 1) * NS] = dN2;
                    /* gamma_12 */
                    B_eq[2 + (c + 0) * NS] = dN2;
                    B_eq[2 + (c + 1) * NS] = dN1;
                }

                for (int j = 0; j < n_local_dof; j++)
                {
                    for (int ii = 0; ii < NS; ii++)
                    {
                        double s = 0.0;
                        for (int kk = 0; kk < NS; kk++) {
                            s += D_eq[ii + kk * NS] * B_eq[kk + j * NS];
                        }
                        tmp[ii + j * NS] = s;
                    }
                }

                for (int j = 0; j < n_local_dof; j++)
                {
                    for (int ii = 0; ii < n_local_dof; ii++)
                    {
                        double s = 0.0;
                        for (int kk = 0; kk < NS; kk++) {
                            s += B_eq[kk + ii * NS] * tmp[kk + j * NS];
                        }
                        Ke[ii + j * n_local_dof] += s;
                    }
                }
            }

            const int64_t *smap = scatter_map + e * (mwSize)n_ld2;
            for (int k = 0; k < n_ld2; k++)
            {
                const int64_t idx = smap[k];
                if (idx > 0)
                {
#pragma omp atomic update
                    V_tang[idx - 1] += Ke[k];
                }
            }
        }

        free(Ke);
        free(B_eq);
        free(D_eq);
        free(tmp);
    }
}
