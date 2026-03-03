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
first_ijv_setup = true;

% Damping parameters (same style as NEWTON.newton).
it_damp_max = 10;

% Pre-compute gravity gradient (constant across iterations).
q2 = Bw * coord(2,:)';
wc2 = repmat(wc, dim, 1);

% Pre-restrict operators to free DOFs (done once, avoids full n_n assembly).
Bw_Q   = Bw(:, Q_w);          % (dim*n_int x n_Qw)
C_Q    = C(:, Q_w);           % (n_int     x n_Qw)
K_D_QQ = K_D(Q_w, Q_w);       % (n_Qw x n_Qw)

while true
    it = it + 1;

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
        reshape(grho * q2, dim, n_int);
    iE = reshape(1:dim*n_int, dim, n_int);
    jE = repmat(1:n_int, dim, 1);
    E  = sparse(iE, jE, vE);

    % ------------------------------------------------------------------
    %  Tangent matrix K_QQ and right-hand side f_Q (restricted to free DOFs)
    % ------------------------------------------------------------------
    K_QQ = K_D_QQ + Bw_Q' * E * C_Q;

    kappa = repmat(perm_r, dim, 1);
    q1 = Bw * pw';
    q3 = q1 + grho * kappa(:) .* q2;
    f_Q = -Bw_Q' * (wc2(:) .* q3);

    % ------------------------------------------------------------------
    %  Solve for Newton increment
    % ------------------------------------------------------------------
    dp = zeros(1, n_n);

    if use_iterative
        % IJV for HYPRE.  The sparsity pattern can change (E depends on
        % the active zone), so we re-do setup_preconditioner_ijv each time.
        [K_I, K_J, K_V] = find(K_QQ);
        if first_ijv_setup
            linear_system_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Qw);
            first_ijv_setup = false;
        else
            linear_system_solver.setup_preconditioner_ijv(K_I, K_J, K_V, n_Qw);
        end
        K_QQ_rsb = sparsersb(K_QQ);   % fast matvec for GMRES
        linear_system_solver.A_orthogonalize(K_QQ_rsb);
        dp(Q_w) = linear_system_solver.solve(K_QQ_rsb, f_Q);
    else
        dp(Q_w) = K_QQ \ f_Q;
    end

    % ------------------------------------------------------------------
    %  Damping (bisection on residual norm, analogous to damping_ALG5)
    % ------------------------------------------------------------------
    criterion = norm(f_Q);
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
        f_a_Q   = -Bw_Q' * (wc2(:) .* q3_a);

        if norm(f_a_Q) < criterion
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
    crit = norm(dp) / norm(pw_init);
    fprintf('  seepage newton it=%d  resid=%.2e  alpha=%.2f\n', it, crit, alpha);

    if crit < tol
        fprintf('Seepage Newton converges: iteration = %d, stopping criterion = %e\n', it, crit);
        break;
    end
    if it > it_max
        warning('Seepage Newton does not converge.');
        fprintf('  stopping criterion = %e\n', crit);
        break;
    end
end

end