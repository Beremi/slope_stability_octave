function pw = newton_flow(pw_init, conduct0, Q_w, weight, Bw, C, K_D, wc, ...
    elem, coord, HatP, WF, eps_int, grho, it_max, tol, linear_system_solver)
%--------------------------------------------------------------------------
% newton_flow applies a damped Newton method to the nonlinear unconfined
% seepage problem.
%
% When a DFGMRES / HYPRE solver object is supplied the tangent system is
% solved iteratively (with IJV pattern-reuse, deflation, and damping
% analogous to NEWTON.newton).  Otherwise the classical backslash (\) path
% is used as a fallback.
%
% INPUT ARGUMENTS:
%   pw_init    - Initial pressure field (1 x n_n).
%   conduct0   - Hydraulic conductivity at integration points (1 x n_int).
%   Q_w        - Logical free-DOF mask for pressure (1 x n_n).
%   weight     - Quadrature weights at integration points.
%   Bw         - Gradient operator (dim*n_int x n_n).
%   C          - Function-value operator (n_int x n_n).
%   K_D        - Darcy stiffness matrix (n_n x n_n, sparse).
%   wc         - weight .* conduct0 product (1 x n_int).
%   elem       - Element connectivity.
%   coord      - Nodal coordinates (dim x n_n).
%   HatP       - Basis function values at quadrature points.
%   WF         - Quadrature weights for volume integration.
%   eps_int    - Penalty parameters at integration points.
%   grho       - Specific weight of water.
%   it_max     - Maximum number of Newton iterations.
%   tol        - Relative tolerance for convergence.
%   linear_system_solver - (optional) DFGMRES handle or [].
%                          If empty / missing, the direct solver (\) is used.
%
% OUTPUT:
%   pw         - Converged pressure field (1 x n_n).
%--------------------------------------------------------------------------

use_iterative = (nargin >= 17) && ~isempty(linear_system_solver);
can_reuse_ijv = use_iterative && ...
    ~isempty(linear_system_solver.preconditioner_initializator) && ...
    ~isempty(linear_system_solver.preconditioner_updater);

% Dimensions.
n_n   = size(coord, 2);
dim   = size(coord, 1);
n_e   = size(elem, 2);
n_q   = length(WF);
n_int = n_e * n_q;
n_p   = size(elem, 1);
n_Qw  = sum(Q_w);

% Initialisation.
pw = pw_init;
it = 0;
status_msg = 'unknown';
last_progress_chars = 0;
linear_iters_total = 0;
t_newton_start = tic;
norm_ref = max(1, norm(pw_init));
rel_update = NaN;

% Damping parameters (same style as NEWTON.newton).
it_damp_max = 10;

% Pre-compute gravity gradient (constant across iterations).
q2 = Bw * coord(2,:)';
wc2 = repmat(wc, dim, 1);

while true
    it = it + 1;
    t_step = tic;

    % ------------------------------------------------------------------
    %  Evaluate relative permeability and its derivative
    % ------------------------------------------------------------------
    pw_e   = reshape(pw(elem(:)), n_p, n_e);
    pw_int = sum(repmat(HatP, 1, n_e) .* kron(pw_e, ones(1, n_q)));

    perm_r     = ones(1, n_int);
    perm_r_der = zeros(1, n_int);
    part1 = (pw_int < eps_int) & (pw_int > 0);
    part2 = (pw_int <= 0);
    perm_r(part1)     = pw_int(part1) ./ eps_int(part1);
    perm_r(part2)     = 0;
    perm_r_der(part1) = 1 ./ eps_int(part1);

    % ------------------------------------------------------------------
    %  Assemble constitutive matrix E  (dim*n_int x n_int)
    % ------------------------------------------------------------------
    vE = repmat(conduct0 .* perm_r_der .* weight, dim, 1) .* ...
        reshape(Bw * (grho * coord(2,:)'), dim, n_int);
    iE = reshape(1:dim*n_int, dim, n_int);
    jE = repmat(1:n_int, dim, 1);
    E  = sparse(iE, jE, vE);

    % ------------------------------------------------------------------
    %  Tangent matrix K and right-hand side f
    % ------------------------------------------------------------------
    K = K_D + Bw' * E * C;

    kappa = repmat(perm_r, dim, 1);
    q1 = Bw * pw';
    q3 = q1 + grho * kappa(:) .* q2;
    f  = -Bw' * (wc2(:) .* q3);

    % ------------------------------------------------------------------
    %  Solve for Newton increment
    % ------------------------------------------------------------------
    dp = zeros(1, n_n);

    if use_iterative
        K_QQ = K(Q_w, Q_w);
        rhs  = f(Q_w);
        if can_reuse_ijv
            % IJV for HYPRE.  The sparsity pattern can change (E depends on
            % the active zone), so we re-do setup_preconditioner_ijv each time.
            [K_I, K_J, K_V] = find(K_QQ);
            linear_system_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Qw);
        else
            linear_system_solver.setup_preconditioner(K_QQ);
        end
        linear_system_solver.A_orthogonalize(K_QQ);
        [dp_q, lin_it] = linear_system_solver.solve(K_QQ, rhs);
        dp(Q_w) = dp_q;
        linear_iters_total = linear_iters_total + lin_it;
    else
        dp(Q_w) = K(Q_w, Q_w) \ f(Q_w);
        lin_it = 0;
    end

    % ------------------------------------------------------------------
    %  Damping (bisection on residual norm, analogous to damping_ALG5)
    % ------------------------------------------------------------------
    criterion = norm(f(Q_w));
    alpha     = 1;
    it_damp   = 0;

    while it_damp < it_damp_max
        it_damp  = it_damp + 1;
        pw_alpha = pw + alpha * dp;

        % Residual at trial point.
        pw_e_a   = reshape(pw_alpha(elem(:)), n_p, n_e);
        pw_int_a = sum(repmat(HatP, 1, n_e) .* kron(pw_e_a, ones(1, n_q)));
        perm_r_a = ones(1, n_int);
        pa1 = (pw_int_a < eps_int) & (pw_int_a > 0);
        perm_r_a(pa1) = pw_int_a(pa1) ./ eps_int(pa1);
        perm_r_a(pw_int_a <= 0) = 0;

        kappa_a = repmat(perm_r_a, dim, 1);
        q3_a    = Bw * pw_alpha' + grho * kappa_a(:) .* q2;
        f_a     = -Bw' * (wc2(:) .* q3_a);

        if norm(f_a(Q_w)) < criterion
            break;
        end
        alpha = alpha / 2;
    end
    if it_damp >= it_damp_max
        alpha = 0;
    end

    % ------------------------------------------------------------------
    %  Update
    % ------------------------------------------------------------------
    pw = pw + alpha * dp;

    % Expand deflation basis on accepted steps.
    if use_iterative && (alpha > 0)
        linear_system_solver.expand_deflation_basis(dp(Q_w)');
    end

    % ------------------------------------------------------------------
    %  Stopping criteria
    % ------------------------------------------------------------------
    rel_update = norm(dp) / norm_ref;
    progress_line = sprintf('  seepage_newton it=%d  rel_update=%.2e  alpha=%.2f  lin_it=%d  step_time=%.2f s', ...
        it, rel_update, alpha, lin_it, toc(t_step));
    last_progress_chars = local_print_progress(progress_line, last_progress_chars);

    if rel_update < tol
        status_msg = 'converged';
        break;
    end
    if isnan(rel_update)
        status_msg = 'nan_update';
        warning('Seepage Newton produced NaN update norm.');
        break;
    end
    if it > it_max
        status_msg = 'max_iterations';
        warning('Seepage Newton does not converge.');
        break;
    end
end

newton_wall_time = toc(t_newton_start);
local_finish_progress(last_progress_chars);
fprintf(['seepage_newton summary: status=%s, it=%d, rel_update=%e, ', ...
    'lin_it_total=%d, wall_time=%.2f s\n'], ...
    status_msg, it, rel_update, linear_iters_total, newton_wall_time);

end

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
