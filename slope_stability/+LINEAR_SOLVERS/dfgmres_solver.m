function [x, iters, res_hist] = dfgmres_solver(A, b, M, W, maxits, tol, x0, profiler)
%DFGMRES  Deflated Flexible GMRES iterative solver.
%
%   [x, iters, res_hist] = DFGMRES(A, b, M, W, maxits, tol, x0, profiler)
%   solves the linear system A*x = b using a deflated version of the
%   flexible GMRES method.
%
%   The algorithm incorporates deflation through a coarse space provided via
%   the deflation basis matrix W. When W is non-empty, an initial coarse
%   solve is applied along with a projection operator to accelerate
%   convergence. If W is empty, the method reverts to standard FGMRES.
%
%   Inputs:
%       A      - Matrix (or object with * overload) for the matrix-vector product.
%       b      - Right-hand side vector
%       M      - Function handle for the preconditioner, i.e., M(x)
%       W      - Deflation basis matrix. If empty, no deflation is performed.
%       maxits - Maximum number of iterations (default: 1000)
%       tol    - Relative residual tolerance (default: 1e-6)
%       x0     - Initial guess (default: zero vector)
%       profiler - Optional PROFILING.Profiler handle ([] = disabled)
%
%   Outputs:
%       x        - Computed solution vector
%       iters    - Number of iterations performed
%       res_hist - History of the relative residual norms
%
%   See also GMRES, PCG.

if nargin < 8, profiler = []; end

% --- Timing accumulators ---
t_precond    = 0;
t_project    = 0;
t_matvec     = 0;
t_ortho      = 0;
t_leastsq    = 0;
n_matvecs    = 0;
n_prec_applies = 0;

t_init_start = tic;

if isempty(maxits)
    maxits = 1000;
end
if isempty(tol)
    tol = 1e-6;
end
if isempty(x0)
    x = zeros(size(b));
else
    x = x0;
end

proj_has_matvec = false;
if isempty(W)
    Proj_fct = @(x) x;
else
    % Define coarse solve and projection operators based on deflation basis W.
    Q_coarse_solve = @(x) W * (W' * x);
    Proj_fct = @(x) x - W * (W' * (A * x));
    proj_has_matvec = true;
    x = Q_coarse_solve(b);
end

n = length(b);
r0 = b - A*x;
n_matvecs = n_matvecs + 1;
b_norm = norm(b);
if b_norm == 0, b_norm = 1; end
res_norm = norm(r0);

% Preallocate the residual history vector
res_hist = zeros(maxits+1, 1);
res_hist(1) = res_norm / b_norm;

if res_hist(1) < tol
    iters = 0;
    res_hist = res_hist(1);
    t_init = toc(t_init_start);
    timing = struct('t_init', t_init, 't_precond', 0, 't_project', 0, ...
        't_matvec', 0, 't_ortho', 0, 't_leastsq', 0, ...
        't_reconstruct', 0, 'n_iters', 0, ...
        'n_matvecs', n_matvecs, 'n_prec_applies', 0);  %#ok<NASGU>
    if ~isempty(profiler)
        profiler.add_time('DFGMRES.solve/init', t_init);
        profiler.add_count('DFGMRES.solve.n_matvecs', n_matvecs);
    end
    return;
end

% Preallocate arrays for Arnoldi vectors, Hessenberg matrix, and preconditioned vectors
V = zeros(n, maxits+1);
H = zeros(maxits+1, maxits);
W = zeros(n, maxits);

V(:,1) = r0 / res_norm;
g = zeros(maxits+1, 1);
g(1) = res_norm;

t_init = toc(t_init_start);

for j = 1:maxits
    % Apply preconditioner to the j-th Krylov basis vector
    t_tmp = tic;
    w = M(V(:,j));
    t_precond = t_precond + toc(t_tmp);
    n_prec_applies = n_prec_applies + 1;

    % Apply deflation projection
    t_tmp = tic;
    w = Proj_fct(w);
    t_project = t_project + toc(t_tmp);
    if proj_has_matvec
        n_matvecs = n_matvecs + 1;
    end

    % Compute A*w
    t_tmp = tic;
    u = A*w;
    t_matvec = t_matvec + toc(t_tmp);
    n_matvecs = n_matvecs + 1;

    % Arnoldi process (orthogonalisation)
    t_tmp = tic;
    for i = 1:j
        H(i,j) = V(:,i)' * u;
        u = u - H(i,j) * V(:,i);
    end
    H(j+1,j) = norm(u);
    t_ortho = t_ortho + toc(t_tmp);

    W(:,j) = w;  % store the preconditioned vector

    % Check for happy breakdown
    if H(j+1,j) < 1e-14
        iters = j;
        break;
    end

    V(:,j+1) = u / H(j+1,j);

    % Solve the least-squares problem and check convergence
    t_tmp = tic;
    y = H(1:j+1, 1:j) \ g(1:j+1);

    % Compute the true residual norm explicitly
    res_norm = norm(g(1:j+1) - H(1:j+1, 1:j) * y);
    res_hist(j+1) = res_norm / b_norm;
    t_leastsq = t_leastsq + toc(t_tmp);

    if res_hist(j+1) < tol
        iters = j;
        break;
    end
    iters = j;  % update iteration count
end

% Trim the residual history to the iterations performed
res_hist = res_hist(1:iters+1);

% Reconstruct the final solution: x = x0 + sum_{j=1}^{iters} y(j)*W(:,j)
t_tmp = tic;
x = x + W(:,1:iters) * y;
t_reconstruct = toc(t_tmp);

% Feed profiler if provided
if ~isempty(profiler)
    profiler.add_time('DFGMRES.solve/init', t_init);
    profiler.add_time('DFGMRES.solve/precond_apply', t_precond);
    profiler.add_time('DFGMRES.solve/deflation_project', t_project);
    profiler.add_time('DFGMRES.solve/matvec', t_matvec);
    profiler.add_time('DFGMRES.solve/arnoldi_ortho', t_ortho);
    profiler.add_time('DFGMRES.solve/least_squares', t_leastsq);
    profiler.add_time('DFGMRES.solve/reconstruct', t_reconstruct);
    profiler.add_count('DFGMRES.solve.n_gmres_iters', iters);
    profiler.add_count('DFGMRES.solve.n_matvecs', n_matvecs);
    profiler.add_count('DFGMRES.solve.n_prec_applies', n_prec_applies);
end

end
