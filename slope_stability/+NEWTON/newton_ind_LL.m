function [U_it, t_it, flag_N, it] = newton_ind_LL( ...
    U_ini, t_ini, omega, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f, ...
    constitutive_matrix_builder, linear_solver)
%--------------------------------------------------------------------------
% newton_ind_LL applies a semismooth Newton method to solve:
%
%       F(U) = t * f,    and    f' * U = omega,
%
% where F(U) is the internal force vector depending on the displacement U.
%
% The method adaptively updates both the displacement field U and the load
% factor t, employing a line search damping strategy and regularization
% based on a mixture of the elastic stiffness and the consistent tangent
% operator.
%
% INPUT ARGUMENTS:
%   U_ini               - Initial displacement field.
%   t_ini               - Initial load factor.
%   omega               - Prescribed control parameter.
%   it_newt_max         - Maximum number of Newton iterations.
%   it_damp_max         - Maximum number of damping (line search) steps.
%   tol                 - Convergence tolerance (relative).
%   r_min               - Base regularization factor for the stiffness matrix.
%   K_elast             - Elastic stiffness matrix.
%   Q                   - Logical array (mask) for restricted Dirichlet dofs.
%   f                   - External load vector.
%   constitutive_matrix_builder - Object to build constitutive matrices.
%   linear_solver       - Linear solver object with deflation capabilities.
%
% OUTPUT:
%   U_it  - Updated displacement field upon convergence or termination.
%   t_it  - Updated load factor.
%   flag_N- Newton success/failure flag (0: success, 1: failure).
%   it    - Number of Newton iterations performed.
%
%--------------------------------------------------------------------------
%
% Basic size and geometry definitions.
n_n = size(U_ini, 2);    % Number of nodes.
dim = size(U_ini, 1);    % Dimension (2D or 3D).

% Initialize auxiliary vectors for iterative updates.
V = zeros(dim, n_n);
W = zeros(dim, n_n);

% Initialize iteration variables.
U_it  = U_ini;         % Current displacement.
t_it  = t_ini;         % Current load factor.
flag_N = 0;            % Success flag: 0 if converged, 1 if failure.
it    = 0;             % Newton iteration counter.
norm_f = norm(f(Q));   % Norm of the restricted external load.
r = r_min;             % Regularization parameter.

compute_diffs = 1;     % Flag to decide whether to recompute differences.
rel_resid = NaN;
alpha = NaN;
status_msg = 'unknown';
last_progress_chars = 0;
lin_it_w_total = 0;
lin_it_v_total = 0;
t_newton_start = tic;

% One-time: initialize K_r(Q,Q) sparsity pattern on free DOFs.
if ~constitutive_matrix_builder.pattern_initialized
    constitutive_matrix_builder.init_K_r_pattern(Q);
end
n_Q = constitutive_matrix_builder.n_Q;

% Precompute upper-triangle mask for sparsersb "unique","sym" constructor.
upper_mask = (constitutive_matrix_builder.ref_I <= constitutive_matrix_builder.ref_J);
K_I_sym = constitutive_matrix_builder.ref_I(upper_mask);
K_J_sym = constitutive_matrix_builder.ref_J(upper_mask);

% IJV preconditioner path (available for HYPRE BoomerAMG only).
can_reuse_ijv = ~isempty(linear_solver.preconditioner_initializator) && ...
                ~isempty(linear_solver.preconditioner_updater);
ijv_pattern_initialized = false;

% Start semismooth Newton iteration.
while true
    it = it + 1;
    t_step = tic;
    
    % Compute internal force vector and constitutive tangent moduli.
    if compute_diffs
        F_int = constitutive_matrix_builder.build_F_and_DS_reduced(U_it);
    end

    % Compute convergence criterion.
    criterion = norm(F_int(Q) - t_it * f(Q));
    rel_resid = criterion / norm_f;
    if (criterion < tol * norm_f) && (it > 1)
        status_msg = 'converged';
        break;
    end

    % Build K_r(Q,Q) as [I, J, V] triplets and construct symmetric matrix.
    [K_I, K_J, K_V] = constitutive_matrix_builder.build_K_r_QQ_vals(r);
    K_rQQ = sparsersb(K_I_sym, K_J_sym, K_V(upper_mask), n_Q, n_Q, "unique", "sym");
    
    % Setup/update preconditioner from IJV when available, otherwise fallback.
    if can_reuse_ijv
        if ~ijv_pattern_initialized
            linear_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Q);
            ijv_pattern_initialized = true;
        else
            linear_solver.update_preconditioner_values(K_V);
        end
    else
        linear_solver.setup_preconditioner(K_rQQ);
    end

    linear_solver.A_orthogonalize(K_rQQ);
    
    % Solve for W(Q): used for load increment computation.
    [W_Q, lin_it_w] = linear_solver.solve(K_rQQ, f(Q));
    W(Q) = W_Q;
    % Solve for V(Q): represents the displacement correction.
    [V_Q, lin_it_v] = linear_solver.solve(K_rQQ, t_it * f(Q) - F_int(Q));
    V(Q) = V_Q;
    lin_it_w_total = lin_it_w_total + lin_it_w;
    lin_it_v_total = lin_it_v_total + lin_it_v;
    
    % Compute the load increment d_t.
    d_t = - (f(Q)' * V(Q)) / (f(Q)' * W(Q));
    % Compute the displacement increment dU.
    dU  = V + d_t * W;
    
    % Perform a line search damping to obtain a suitable step length alpha.
    alpha = NEWTON.damping(it_damp_max, U_it, dU, F_int, F_int*0, constitutive_matrix_builder);
    
    % Adjust regularization based on the damping result.
    compute_diffs = 1;
    if alpha < 1e-1
        if alpha == 0
            compute_diffs = 0;  % Skip recomputing K_tangent if step rejected.
            r = r * 2;
        else
            r = r * 2^(1/4);
        end
    else
        % Expand deflation bases with current correction vectors.
        linear_solver.expand_deflation_basis(W(Q));
        linear_solver.expand_deflation_basis(V(Q));
        
        if alpha > 0.5
            r = max(r / sqrt(2), r_min);
        end
    end
    if alpha == 0 && r > 1
        status_msg = 'stalled_regularization';
        flag_N = 1;
        break;
    end

    % Update displacement and load factor.
    U_it = U_it + alpha * dU;
    % Enforce the constraint f' * U = omega.
    U_it = omega * U_it / (f(Q)' * U_it(Q));
    t_it = t_it + d_t;

    progress_line = sprintf(['  newton_ind_LL it=%d  rel_resid=%.2e  alpha=%.2f  r=%.3g', ...
        '  lin_it=[W:%d,V:%d]  step_time=%.2f s'], ...
        it, rel_resid, alpha, r, lin_it_w, lin_it_v, toc(t_step));
    last_progress_chars = local_print_progress(progress_line, last_progress_chars);
    
    % Check if maximum iterations reached.
    if isnan(rel_resid) || (it == it_newt_max)
        flag_N = 1;
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
fprintf(['newton_ind_LL summary: status=%s, it=%d, rel_resid=%e, ', ...
    'lin_it_total=[W:%d,V:%d], wall_time=%.2f s\n'], ...
    status_msg, it, rel_resid, lin_it_w_total, lin_it_v_total, newton_wall_time);

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
