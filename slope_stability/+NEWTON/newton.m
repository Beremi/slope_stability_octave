function [U_it, flag_N, it, history] = newton(U_ini, tol, it_newt_max, it_damp_max, r_min, ...
    K_elast, Q, f, constitutive_matrix_builder, linear_system_solver)
%--------------------------------------------------------------------------
% newton applies a semismooth Newton method with damping and regularization
% to solve the nonlinear system:
%
%       F(U) = f,
%
% where F(U) is the internal force vector computed from the displacement field U.
%
% The method updates U iteratively until the relative residual
% norm(F(Q)-f(Q))/norm(f(Q)) is below the prescribed tolerance tol.
%
% INPUT ARGUMENTS:
%   U_ini               - Initial displacement field, size (dim, n_n).
%   tol                 - Relative tolerance for Newton's convergence.
%   it_newt_max         - Maximum number of Newton iterations.
%   it_damp_max         - Maximum number of damping iterations.
%   r_min               - Base regularization parameter.
%   K_elast             - (unused, kept for API compatibility; elastic K
%                          is precomputed inside constitutive_matrix_builder).
%   Q                   - Logical free-DOF mask (dim x n_n).
%   f                   - External load vector.
%   constitutive_matrix_builder - CONSTITUTIVE object (pattern must be
%                         initialised via init_K_r_pattern before calling).
%   linear_system_solver- DFGMRES solver object with IJV support.
%
% OUTPUTS:
%   U_it    - Computed displacement field approximating the solution.
%   flag_N  - Newton solver flag: 0 if converged, 1 if not converged.
%   it      - Number of Newton iterations performed.
%   history - Struct (currently empty; timing is handled by PROFILING.Profiler).
%--------------------------------------------------------------------------
%
% Initialization.
n_n = size(U_ini, 2);         % Number of nodes.
dim = size(U_ini, 1);         % Spatial dimension.
dU = zeros(dim, n_n);         % Newton increment.
U_it = U_ini;                 % Current displacement.
it = 0;                       % Iteration counter.
flag_N = 0;                   % Convergence flag.
r = r_min;                    % Regularization parameter.
compute_diffs = 1;            % Flag: must recompute stress/tangent?
norm_f = norm(f(Q));

% Ensure pattern is initialised.
if ~constitutive_matrix_builder.pattern_initialized
    constitutive_matrix_builder.init_K_r_pattern(Q);
end
n_Q = constitutive_matrix_builder.n_Q;

% Precompute upper-triangle mask for sparsersb "unique","sym" constructor.
upper_mask = (constitutive_matrix_builder.ref_I <= constitutive_matrix_builder.ref_J);
K_I_sym = constitutive_matrix_builder.ref_I(upper_mask);
K_J_sym = constitutive_matrix_builder.ref_J(upper_mask);

% IJV preconditioner path (available for HYPRE BoomerAMG only).
can_reuse_ijv = ~isempty(linear_system_solver.preconditioner_initializator) && ...
                ~isempty(linear_system_solver.preconditioner_updater);
ijv_pattern_initialized = false;

% Preallocate history arrays.
residual_history = nan(1, it_newt_max);
r_history = nan(1, it_newt_max);
alpha_history = nan(1, it_newt_max);
rel_resid = NaN;
alpha = NaN;
status_msg = 'unknown';
last_progress_chars = 0;
linear_iters_total = 0;
t_newton_start = tic;

% Newton's iteration loop.
while true
    it = it + 1;
    t_step = tic;

    % Compute internal force and tangent moduli.
    if compute_diffs
        F = constitutive_matrix_builder.build_F_and_DS_reduced(U_it);
        criterion = norm(F(Q) - f(Q));
        rel_resid = criterion / norm_f;
        residual_history(it) = rel_resid;
        if (rel_resid < tol)
            status_msg = 'converged';
            break;
        end
    end

    r_history(it) = r;

    % Build K_r(Q,Q) as [I, J, V] triplets  →  sparsersb.
    [K_I, K_J, K_V] = constitutive_matrix_builder.build_K_r_QQ_vals(r);

    K_rQQ = sparsersb(K_I_sym, K_J_sym, K_V(upper_mask), n_Q, n_Q, "unique", "sym");

    % Setup preconditioner from IJV when available, otherwise fallback to A-based setup.
    if can_reuse_ijv
        if ~ijv_pattern_initialized
            linear_system_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Q);
            ijv_pattern_initialized = true;
        else
            linear_system_solver.update_preconditioner_values(K_V);
        end
    else
        linear_system_solver.setup_preconditioner(K_rQQ);
    end

    linear_system_solver.A_orthogonalize(K_rQQ);

    % Solve for the Newton increment.
    [dU_Q, lin_it] = linear_system_solver.solve(K_rQQ, f(Q) - F(Q));
    dU(Q) = dU_Q;
    linear_iters_total = linear_iters_total + lin_it;

    % Determine damping factor via line search.
    alpha = NEWTON.damping(it_damp_max, U_it, dU, F, f, constitutive_matrix_builder);
    alpha_history(it) = alpha;

    % Regularization adjustments based on alpha.
    compute_diffs = 1;
    if alpha < 1e-1
        if alpha == 0
            compute_diffs = 0;  % Skip recomputation if step is rejected.
            r = r * 2;
        else
            r = r * 2^(1/4);
        end
    else
        % Expand deflation basis with the computed increment.
        linear_system_solver.expand_deflation_basis(dU(Q));
        if alpha > 0.5
            r = max(r / sqrt(2), r_min);
        end
    end

    if (alpha == 0) && (r > 1)
        status_msg = 'stalled_regularization';
        if rel_resid > 10 * tol
            flag_N = 1;
        end
        break;
    end

    % Update displacement.
    U_it = U_it + alpha * dU;

    progress_line = sprintf('  newton it=%d  rel_resid=%.2e  alpha=%.2f  r=%.3g  lin_it=%d  step_time=%.2f s', ...
        it, rel_resid, alpha, r, lin_it, toc(t_step));
    last_progress_chars = local_print_progress(progress_line, last_progress_chars);

    % Check for divergence or maximum iterations.
    if isnan(rel_resid) || (it == it_newt_max)
        if rel_resid > 10 * tol
            flag_N = 1;
        end
        if isnan(rel_resid)
            status_msg = 'nan_residual';
        else
            status_msg = 'max_iterations';
        end
        break;
    end

end % while

newton_wall_time = toc(t_newton_start);
local_finish_progress(last_progress_chars);
fprintf('newton summary: status=%s, it=%d, rel_resid=%e, lin_it_total=%d, wall_time=%.2f s\n', ...
    status_msg, it, rel_resid, linear_iters_total, newton_wall_time);

% Build output history.
history.residual = residual_history(1:it);
history.r = r_history(1:it);
history.alpha = alpha_history(1:it);
history.linear_iters_total = linear_iters_total;
history.wall_time = newton_wall_time;

end % function

function last_chars = local_print_progress(msg, last_chars)
pad = max(0, last_chars - numel(msg));
fprintf('\r%s%s', msg, repmat(' ', 1, pad));
fflush(stdout);
last_chars = numel(msg);
end

function local_finish_progress(last_chars)
if last_chars > 0
    fprintf('\r%s\r', repmat(' ', 1, last_chars));
end
end
