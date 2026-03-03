function [U_it, lambda_it, flag_N, it, history] = newton_ind_SSR(U_ini, omega, lambda_ini, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, constitutive_matrix_builder, linear_solver)
%--------------------------------------------------------------------------
% newton_ind_SSR solves the system
%
%       F_lambda(U) = f   and   f' * U = omega,
%
% using a semismooth Newton method with regularization.
%
% OUTPUTS:
%   U_it      - Updated displacement field after convergence or termination.
%   lambda_it - Updated value of lambda.
%   flag_N    - Newton success/failure flag (0: success, 1: failure).
%   it        - Number of Newton iterations performed.
%   history   - Struct with fields:
%                 .residual  - History of relative residuals.
%                 .r         - History of regularization parameters.
%                 .alpha     - History of damping factors.
%--------------------------------------------------------------------------

% Initialization.
n_n = size(U_ini, 2);
dim = size(U_ini, 1);
V = zeros(dim, n_n);
W = zeros(dim, n_n);
U_it = U_ini;
lambda_it = lambda_ini;
eps = tol / 1000;
norm_f = norm(f(Q));
flag_N = 0;
compute_diffs = 1;
r = r_min;
it = 0;

% Preallocate history arrays
residual_history = nan(1, it_newt_max);
r_history = nan(1, it_newt_max);
alpha_history = nan(1, it_newt_max);

% One-time: restrict B to free DOFs and build maximal sparsity pattern.
if ~constitutive_matrix_builder.pattern_initialized
    constitutive_matrix_builder.init_K_r_pattern(Q);
end
n_Q = constitutive_matrix_builder.n_Q;

% Precompute upper-triangle mask for sparsersb "unique","sym" constructor.
% ref_I/ref_J are the same every call, so we mask once.
upper_mask = (constitutive_matrix_builder.ref_I <= constitutive_matrix_builder.ref_J);
K_I_sym = constitutive_matrix_builder.ref_I(upper_mask);
K_J_sym = constitutive_matrix_builder.ref_J(upper_mask);

% Flag for HYPRE IJV pattern reuse: first call does full setup,
% subsequent calls only update values (same sparsity pattern).
hypre_pattern_initialized = false;

% Semismooth Newton loop.
while true
    it = it + 1;
    t_step = tic;

    %% Compute constitutive response (F and DS, skip sparse K_tangent).
    if compute_diffs
        F = constitutive_matrix_builder.build_F_and_DS_all(lambda_it, U_it);
        criterion = norm(F(Q) - f(Q));
        rel_resid = criterion / norm_f;

        residual_history(it) = rel_resid;

        if (rel_resid < tol) && (it > 1)
            fprintf('Newton method converges: iteration = %d, rel. resid = %e\n', it, rel_resid);
            break;
        end
    end

    r_history(it) = r;

    %% Build K_r(Q,Q) as [I, J, V] triplets, then sparsersb matrix.
    [K_I, K_J, K_V] = constitutive_matrix_builder.build_K_r_QQ_vals(r);

    K_rQQ = sparsersb(K_I_sym, K_J_sym, K_V(upper_mask), n_Q, n_Q, "unique", "sym");

    %% Estimate Newton increment.
    lambda_eps = lambda_it + eps;
    F_eps = constitutive_matrix_builder.build_F_all(lambda_eps, U_it);
    G = (F_eps - F) / eps;

    %% Setup/update preconditioner (HYPRE reuses sparsity pattern after first call).
    if ~hypre_pattern_initialized
        linear_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Q);
        hypre_pattern_initialized = true;
    else
        linear_solver.update_preconditioner_values(K_V);
    end

    linear_solver.A_orthogonalize(K_rQQ);

    W(Q) = linear_solver.solve(K_rQQ, -G(Q));
    V(Q) = linear_solver.solve(K_rQQ, f(Q) - F(Q));

    d_l = - (f(Q)' * V(Q)) / (f(Q)' * W(Q));
    d_U = V + d_l * W;

    %% Line search damping.
    alpha = NEWTON.damping_ALG5(it_damp_max, U_it, lambda_it, d_U, d_l, f, criterion, Q, constitutive_matrix_builder);
    alpha_history(it) = alpha;

    %% Regularization adjustments.
    compute_diffs = 1;
    if alpha < 1e-1
        if alpha == 0
            compute_diffs = 0;
            r = r * 2;
        else
            r = r * 2^(1/4);
        end
    else
        linear_solver.expand_deflation_basis(W(Q));
        linear_solver.expand_deflation_basis(V(Q));
        if alpha > 0.5
            r = max(r / sqrt(2), r_min);
        end
    end

    if (alpha == 0) && (r > 1)
        fprintf('\nNewton solver does not converge: iteration = %d, stopping criterion=%e\n', it, rel_resid);
        if rel_resid > 10 * tol
            flag_N = 1;
        end
        break;
    end

    %% Update variables.
    U_it = U_it + alpha * d_U;
    U_it = omega * U_it / (f(Q)' * U_it(Q));  % Enforce constraint
    lambda_it = lambda_it + alpha * d_l;

    fprintf('  newton_ind_SSR it=%d  resid=%.2e  alpha=%.2f  r=%.3g  step_time=%.2f s\n', ...
        it, rel_resid, alpha, r, toc(t_step));

    %% Check maximum iterations.
    if it == it_newt_max
        fprintf('Newton solver does not converge: iteration = %d, stopping criterion=%e\n', it, rel_resid);
        if rel_resid > 10 * tol
            flag_N = 1;
        end
        break;
    end
end % while

% Trim unused history entries
history.residual = residual_history(1:it);
history.r = r_history(1:it);
history.alpha = alpha_history(1:it);

end % function
