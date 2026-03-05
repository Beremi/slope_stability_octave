/*
 * constitutive_2D_kernel.h
 *
 * Per-integration-point kernel for 2D Mohr-Coulomb elastic-perfectly-plastic
 * constitutive model (plane strain-like 2D formulation used in
 * constitutive_problem_2D.m).
 *
 * Provides:
 *   constitutive_2D_point(E_in, c_bar, sin_phi, shear, bulk, lame, S, DS)
 *
 * Inputs:
 *   E_in  - 3 components (e11, e22, gamma12)
 * Outputs:
 *   S     - 4 components
 *   DS    - 9 components (3x3, column-major) or NULL
 */

#ifndef CONSTITUTIVE_2D_KERNEL_H
#define CONSTITUTIVE_2D_KERNEL_H

#include <math.h>
#include <string.h>

static inline void outer3_scaled(const double *a, double scale, double *out9)
{
    int i, j;
    for (j = 0; j < 3; ++j)
    {
        for (i = 0; i < 3; ++i)
        {
            out9[i + 3 * j] = scale * a[i] * a[j];
        }
    }
}

static inline void constitutive_2D_point(
    const double *E_in,
    double c_bar, double sin_phi, double shear, double bulk, double lame,
    double *S, double *DS)
{
    /* Trial strain in 4-component representation */
    double E_tr[4] = {E_in[0], E_in[1], E_in[2], 0.0};

    /* Inordered eigenvalues of trial strain */
    const double I1 = E_tr[0] + E_tr[1];
    const double I2 = sqrt((E_tr[0] - E_tr[1]) * (E_tr[0] - E_tr[1]) + E_tr[2] * E_tr[2]);
    const double eig0_1 = 0.5 * (I1 + I2);
    const double eig0_2 = 0.5 * (I1 - I2);
    const double eig0_3 = E_tr[3];

    /* Inordered eigenprojections (first derivatives) */
    double Eig0_1[4] = {0.0, 0.0, 0.0, 0.0};
    if (I2 != 0.0)
    {
        Eig0_1[0] = (E_tr[0] - eig0_2) / I2;
        Eig0_1[1] = (E_tr[1] - eig0_2) / I2;
        Eig0_1[2] = 0.5 * E_tr[2] / I2;
    }
    else
    {
        Eig0_1[0] = 1.0;
        Eig0_1[1] = 1.0;
    }

    const double Eig0_2[4] = {1.0 - Eig0_1[0], 1.0 - Eig0_1[1], -Eig0_1[2], 0.0};
    const double Eig0_3[4] = {0.0, 0.0, 0.0, 1.0};

    /* Inordered second derivatives of eigenvalues */
    double EIG0_1[9] = {0.0};
    double EIG0_2[9] = {0.0};
    double EIG0_3[9] = {0.0};

    if (I2 != 0.0)
    {
        const double a1 = Eig0_1[0], a2 = Eig0_1[1], a3 = Eig0_1[2];
        const double b1 = Eig0_2[0], b2 = Eig0_2[1], b3 = Eig0_2[2];
        const double invI2 = 1.0 / I2;

        EIG0_1[0] = (1.0 - a1 * a1 - b1 * b1) * invI2;
        EIG0_1[1] = (-a2 * a1 - b2 * b1) * invI2;
        EIG0_1[2] = (-a3 * a1 - b3 * b1) * invI2;
        EIG0_1[3] = (-a1 * a2 - b1 * b2) * invI2;
        EIG0_1[4] = (1.0 - a2 * a2 - b2 * b2) * invI2;
        EIG0_1[5] = (-a3 * a2 - b3 * b2) * invI2;
        EIG0_1[6] = (-a1 * a3 - b1 * b3) * invI2;
        EIG0_1[7] = (-a2 * a3 - b2 * b3) * invI2;
        EIG0_1[8] = (0.5 - a3 * a3 - b3 * b3) * invI2;

        {
            int k;
            for (k = 0; k < 9; ++k)
            {
                EIG0_2[k] = -EIG0_1[k];
            }
        }
    }

    /* Reordering of eigenvalues and derivatives */
    double eig_1 = eig0_1, eig_2 = eig0_2, eig_3 = eig0_3;
    double Eig_1[4], Eig_2[4], Eig_3[4];
    double EIG_1[9], EIG_2[9], EIG_3[9];
    memcpy(Eig_1, Eig0_1, 4 * sizeof(double));
    memcpy(Eig_2, Eig0_2, 4 * sizeof(double));
    memcpy(Eig_3, Eig0_3, 4 * sizeof(double));
    memcpy(EIG_1, EIG0_1, 9 * sizeof(double));
    memcpy(EIG_2, EIG0_2, 9 * sizeof(double));
    memcpy(EIG_3, EIG0_3, 9 * sizeof(double));

    if ((eig0_1 >= eig0_3) && (eig0_3 > eig0_2))
    {
        eig_2 = eig0_3;
        eig_3 = eig0_2;
        memcpy(Eig_2, Eig0_3, 4 * sizeof(double));
        memcpy(Eig_3, Eig0_2, 4 * sizeof(double));
        memcpy(EIG_2, EIG0_3, 9 * sizeof(double));
        memcpy(EIG_3, EIG0_2, 9 * sizeof(double));
    }

    if (eig0_3 > eig0_1)
    {
        eig_1 = eig0_3;
        eig_2 = eig0_1;
        eig_3 = eig0_2;
        memcpy(Eig_1, Eig0_3, 4 * sizeof(double));
        memcpy(Eig_2, Eig0_1, 4 * sizeof(double));
        memcpy(Eig_3, Eig0_2, 4 * sizeof(double));
        memcpy(EIG_1, EIG0_3, 9 * sizeof(double));
        memcpy(EIG_2, EIG0_1, 9 * sizeof(double));
        memcpy(EIG_3, EIG0_2, 9 * sizeof(double));
    }

    /* Criteria and multipliers */
    const double trace_E = eig_1 + eig_2 + eig_3;
    const double f_tr = 2.0 * shear * ((1.0 + sin_phi) * eig_1 - (1.0 - sin_phi) * eig_3) +
                        2.0 * lame * sin_phi * trace_E - c_bar;

    const double gamma_sl = (eig_1 - eig_2) / (1.0 + sin_phi);
    const double gamma_sr = (eig_2 - eig_3) / (1.0 - sin_phi);
    const double gamma_la = (eig_1 + eig_2 - 2.0 * eig_3) / (3.0 - sin_phi);
    const double gamma_ra = (2.0 * eig_1 - eig_2 - eig_3) / (3.0 + sin_phi);

    const double denom_s = 4.0 * lame * sin_phi * sin_phi +
                           2.0 * shear * (1.0 + sin_phi) * (1.0 + sin_phi) +
                           2.0 * shear * (1.0 - sin_phi) * (1.0 - sin_phi);
    const double denom_l = 4.0 * lame * sin_phi * sin_phi +
                           shear * (1.0 + sin_phi) * (1.0 + sin_phi) +
                           2.0 * shear * (1.0 - sin_phi) * (1.0 - sin_phi);
    const double denom_r = 4.0 * lame * sin_phi * sin_phi +
                           2.0 * shear * (1.0 + sin_phi) * (1.0 + sin_phi) +
                           shear * (1.0 - sin_phi) * (1.0 - sin_phi);
    const double denom_a = 4.0 * bulk * sin_phi * sin_phi;

    const double lambda_s = f_tr / denom_s;
    const double lambda_l = (shear * ((1.0 + sin_phi) * (eig_1 + eig_2) - 2.0 * (1.0 - sin_phi) * eig_3) +
                             2.0 * lame * sin_phi * trace_E - c_bar) / denom_l;
    const double lambda_r = (shear * (2.0 * (1.0 + sin_phi) * eig_1 - (1.0 - sin_phi) * (eig_2 + eig_3)) +
                             2.0 * lame * sin_phi * trace_E - c_bar) / denom_r;
    const double lambda_a = (2.0 * bulk * sin_phi * trace_E - c_bar) / denom_a;
    (void)lambda_a;

    /* Branch selection (same order as MATLAB) */
    const int test_el = (f_tr <= 0.0);
    const int test_s = (!test_el) && (lambda_s <= ((gamma_sl < gamma_sr) ? gamma_sl : gamma_sr));
    const int test_l = (!test_el && !test_s) && (gamma_sl < gamma_sr) &&
                       (lambda_l >= gamma_sl) && (lambda_l <= gamma_la);
    const int test_r = (!test_el && !test_s) && (gamma_sl > gamma_sr) &&
                       (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra);
    const int test_a = (!test_el && !test_s && !test_l && !test_r);

    /* Stress */
    if (test_el)
    {
        const double vol = E_tr[0] + E_tr[1] + E_tr[3];
        const double s2 = 2.0 * shear;
        S[0] = lame * vol + s2 * E_tr[0];
        S[1] = lame * vol + s2 * E_tr[1];
        S[2] = shear * E_tr[2];
        S[3] = lame * vol + s2 * E_tr[3];
    }
    else if (test_s)
    {
        const double sigma_1 = lame * trace_E + 2.0 * shear * eig_1 -
                               lambda_s * (2.0 * lame * sin_phi + 2.0 * shear * (1.0 + sin_phi));
        const double sigma_2 = lame * trace_E + 2.0 * shear * eig_2 -
                               lambda_s * (2.0 * lame * sin_phi);
        const double sigma_3 = lame * trace_E + 2.0 * shear * eig_3 -
                               lambda_s * (2.0 * lame * sin_phi - 2.0 * shear * (1.0 - sin_phi));
        {
            int i;
            for (i = 0; i < 4; ++i)
            {
                S[i] = sigma_1 * Eig_1[i] + sigma_2 * Eig_2[i] + sigma_3 * Eig_3[i];
            }
        }
    }
    else if (test_l)
    {
        const double sigma_1 = lame * trace_E + shear * (eig_1 + eig_2) -
                               lambda_l * (2.0 * lame * sin_phi + shear * (1.0 + sin_phi));
        const double sigma_3 = lame * trace_E + 2.0 * shear * eig_3 -
                               lambda_l * (2.0 * lame * sin_phi - 2.0 * shear * (1.0 - sin_phi));
        {
            int i;
            for (i = 0; i < 4; ++i)
            {
                S[i] = sigma_1 * (Eig_1[i] + Eig_2[i]) + sigma_3 * Eig_3[i];
            }
        }
    }
    else if (test_r)
    {
        const double sigma_1 = lame * trace_E + 2.0 * shear * eig_1 -
                               lambda_r * (2.0 * lame * sin_phi + 2.0 * shear * (1.0 + sin_phi));
        const double sigma_3 = lame * trace_E + shear * (eig_2 + eig_3) -
                               lambda_r * (2.0 * lame * sin_phi - shear * (1.0 - sin_phi));
        {
            int i;
            for (i = 0; i < 4; ++i)
            {
                S[i] = sigma_1 * Eig_1[i] + sigma_3 * (Eig_2[i] + Eig_3[i]);
            }
        }
    }
    else /* apex */
    {
        const double sigma_1 = c_bar / (2.0 * sin_phi);
        S[0] = sigma_1;
        S[1] = sigma_1;
        S[2] = 0.0;
        S[3] = sigma_1;
    }

    if (!DS)
    {
        return;
    }

    /* Tangent DS (3x3, column-major) */
    if (test_el)
    {
        DS[0] = lame + 2.0 * shear;
        DS[1] = lame;
        DS[2] = 0.0;
        DS[3] = lame;
        DS[4] = lame + 2.0 * shear;
        DS[5] = 0.0;
        DS[6] = 0.0;
        DS[7] = 0.0;
        DS[8] = shear;
        return;
    }

    if (test_s)
    {
        const double sigma_1 = lame * trace_E + 2.0 * shear * eig_1 -
                               lambda_s * (2.0 * lame * sin_phi + 2.0 * shear * (1.0 + sin_phi));
        const double sigma_2 = lame * trace_E + 2.0 * shear * eig_2 -
                               lambda_s * (2.0 * lame * sin_phi);
        const double sigma_3 = lame * trace_E + 2.0 * shear * eig_3 -
                               lambda_s * (2.0 * lame * sin_phi - 2.0 * shear * (1.0 - sin_phi));

        double v1[3] = {Eig_1[0], Eig_1[1], Eig_1[2]};
        double v2[3] = {Eig_2[0], Eig_2[1], Eig_2[2]};
        double v3[3] = {Eig_3[0], Eig_3[1], Eig_3[2]};

        double mat3[9], mat4[9], mat5[9], mat6[9];
        double Eig6[3] = {
            2.0 * shear * (1.0 + sin_phi) * v1[0] - 2.0 * shear * (1.0 - sin_phi) * v3[0] + 2.0 * lame * sin_phi,
            2.0 * shear * (1.0 + sin_phi) * v1[1] - 2.0 * shear * (1.0 - sin_phi) * v3[1] + 2.0 * lame * sin_phi,
            2.0 * shear * (1.0 + sin_phi) * v1[2] - 2.0 * shear * (1.0 - sin_phi) * v3[2]};

        outer3_scaled(v1, 2.0 * shear, mat3);
        outer3_scaled(v2, 2.0 * shear, mat4);
        outer3_scaled(v3, 2.0 * shear, mat5);
        outer3_scaled(Eig6, 1.0 / denom_s, mat6);

        {
            int k;
            for (k = 0; k < 9; ++k)
            {
                const double vol_k = (k == 0 || k == 1 || k == 3 || k == 4) ? 1.0 : 0.0;
                DS[k] = sigma_1 * EIG_1[k] + sigma_2 * EIG_2[k] + sigma_3 * EIG_3[k] +
                        lame * vol_k + mat3[k] + mat4[k] + mat5[k] - mat6[k];
            }
        }
        return;
    }

    if (test_l)
    {
        const double sigma_1 = lame * trace_E + shear * (eig_1 + eig_2) -
                               lambda_l * (2.0 * lame * sin_phi + shear * (1.0 + sin_phi));
        const double sigma_3 = lame * trace_E + 2.0 * shear * eig_3 -
                               lambda_l * (2.0 * lame * sin_phi - 2.0 * shear * (1.0 - sin_phi));

        double Eig12[3] = {Eig_1[0] + Eig_2[0], Eig_1[1] + Eig_2[1], Eig_1[2] + Eig_2[2]};
        double v3[3] = {Eig_3[0], Eig_3[1], Eig_3[2]};
        double EIG12[9];
        double mat3[9], mat5[9], mat6[9];
        double Eig6[3] = {
            shear * (1.0 + sin_phi) * Eig12[0] - 2.0 * shear * (1.0 - sin_phi) * v3[0] + 2.0 * lame * sin_phi,
            shear * (1.0 + sin_phi) * Eig12[1] - 2.0 * shear * (1.0 - sin_phi) * v3[1] + 2.0 * lame * sin_phi,
            shear * (1.0 + sin_phi) * Eig12[2] - 2.0 * shear * (1.0 - sin_phi) * v3[2]};

        {
            int k;
            for (k = 0; k < 9; ++k)
            {
                EIG12[k] = EIG_1[k] + EIG_2[k];
            }
        }
        outer3_scaled(Eig12, shear, mat3);
        outer3_scaled(v3, 2.0 * shear, mat5);
        outer3_scaled(Eig6, 1.0 / denom_l, mat6);

        {
            int k;
            for (k = 0; k < 9; ++k)
            {
                const double vol_k = (k == 0 || k == 1 || k == 3 || k == 4) ? 1.0 : 0.0;
                DS[k] = sigma_1 * EIG12[k] + sigma_3 * EIG_3[k] +
                        lame * vol_k + mat3[k] + mat5[k] - mat6[k];
            }
        }
        return;
    }

    if (test_r)
    {
        const double sigma_1 = lame * trace_E + 2.0 * shear * eig_1 -
                               lambda_r * (2.0 * lame * sin_phi + 2.0 * shear * (1.0 + sin_phi));
        const double sigma_3 = lame * trace_E + shear * (eig_2 + eig_3) -
                               lambda_r * (2.0 * lame * sin_phi - shear * (1.0 - sin_phi));

        double v1[3] = {Eig_1[0], Eig_1[1], Eig_1[2]};
        double Eig23[3] = {Eig_2[0] + Eig_3[0], Eig_2[1] + Eig_3[1], Eig_2[2] + Eig_3[2]};
        double EIG23[9];
        double mat3[9], mat5[9], mat6[9];
        double Eig6[3] = {
            2.0 * shear * (1.0 + sin_phi) * v1[0] - shear * (1.0 - sin_phi) * Eig23[0] + 2.0 * lame * sin_phi,
            2.0 * shear * (1.0 + sin_phi) * v1[1] - shear * (1.0 - sin_phi) * Eig23[1] + 2.0 * lame * sin_phi,
            2.0 * shear * (1.0 + sin_phi) * v1[2] - shear * (1.0 - sin_phi) * Eig23[2]};

        {
            int k;
            for (k = 0; k < 9; ++k)
            {
                EIG23[k] = EIG_2[k] + EIG_3[k];
            }
        }
        outer3_scaled(v1, 2.0 * shear, mat3);
        outer3_scaled(Eig23, shear, mat5);
        outer3_scaled(Eig6, 1.0 / denom_r, mat6);

        {
            int k;
            for (k = 0; k < 9; ++k)
            {
                const double vol_k = (k == 0 || k == 1 || k == 3 || k == 4) ? 1.0 : 0.0;
                DS[k] = sigma_1 * EIG_1[k] + sigma_3 * EIG23[k] +
                        lame * vol_k + mat3[k] + mat5[k] - mat6[k];
            }
        }
        return;
    }

    /* Apex */
    memset(DS, 0, 9 * sizeof(double));
}

#endif /* CONSTITUTIVE_2D_KERNEL_H */
