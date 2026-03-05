% Auto-generated from slope_stability_3D_hetero_LL.ipynb.

% Update the notebook first, then re-run scripts/sync_notebooks_to_m.py.

%% 3D Heterogeneous LL (HYPRE + Profiler)

%% Notes Before Running
% - Activate local environment: `source setup/activate_optimized_octave.sh`
% - In Jupyter, select `Octave (local-rsb)` kernel.
% - This notebook follows the COMSOL notebook structure with HYPRE + profiler instrumentation.

%plot -f png -r 600
% Load sparsersb for multithreaded sparse matrix-vector products
pkg load sparsersb;

% Configure toolkit for both GUI and headless Jupyter environments.
try
    graphics_toolkit('qt');
catch
    graphics_toolkit('gnuplot');
end
set(0, 'defaultfigurevisible', 'off');

fprintf('Working directory: %s\n', pwd);
fprintf('Graphics toolkit: %s\n', graphics_toolkit());

%% 1) Main Input Data

elem_type = 'P2';
Davis_type = 'B';
mat_props = [15, 30,  0, 10000, 0.33, 19, 19;  % Cover layer
    15, 38,  0, 50000, 0.30, 22, 22;  % General foundation
    10, 35,  0, 50000, 0.30, 21, 21;  % Relatively weak foundation
    18, 32,  0, 20000, 0.33, 20, 20]; % General slope mass

file_path = 'meshes/LL_hetero_ada_L1.h5';

%% 2) Reference Element Data and Mesh Load

[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

[coord, elem, surf, Q, material_identifier] = MESH.load_mesh_P2(file_path);
[coord, elem, surf, Q] = MESH.reorder_mesh(coord, elem, surf, Q);

n_n = size(coord, 2);
n_unknown = length(coord(Q));
n_e = size(elem, 2);
n_q = length(WF);
n_int = n_e * n_q;

fprintf('\nMesh data:');
fprintf('  number of nodes =%d ', n_n);
fprintf('  number of unknowns =%d ', n_unknown);
fprintf('  number of elements =%d ', n_e);
fprintf('  number of integration points =%d ', n_int);
fprintf('\n');

%% Mesh Preview

plotter = VIZ.SolutionPlotter(coord, elem, surf, 'comsol');
plotter.plot_mesh();
drawnow;

%% 3) Mechanical Material Fields and Assembly

fields = {'c0', 'phi', 'psi', 'young', 'poisson', 'gamma_sat', 'gamma_unsat'};
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

saturation = true(1, n_int);
[c0, phi, psi, shear, bulk, lame, gamma] = ...
    ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);

[K_elast, B, WEIGHT, DPhi1_out, DPhi2_out, DPhi3_out] = ASSEMBLY.elastic_stiffness_matrix_3D(...
    elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
f = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

plotter.set_strain_data(B, Xi);

%% 4) Continuation, Newton, and Linear Solver Parameters

d_t_min = 1e-3;
step_max = 100;
LL_omega_max = 8e7;

it_newt_max = 200;
it_damp_max = 10;
tol = 1e-4;
r_min = 1e-4;

% agmg folder is baked into LINEAR_SOLVERS.set_linear_solver
solver_type = 'DFGMRES_HYPRE_BOOMERAMG';

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
    'use_as_preconditioner', true);

linear_system_solver = LINEAR_SOLVERS.set_linear_solver(solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, coord, boomeramg_opts);

n_strain = 6;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(...
    B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, 3);

constitutive_matrix_builder.set_element_data(elem, DPhi1_out, DPhi2_out, DPhi3_out);

profiler = PROFILING.Profiler();
constitutive_matrix_builder.profiler = profiler;
linear_system_solver.profiler = profiler;

%% 5) Run LL Indirect Continuation

fprintf('\n Indirect continuation method for the LL method\n');
tic;

constitutive_matrix_builder.reduction(1.0);
U_elast = zeros(3, n_n);
linear_system_solver.setup_preconditioner(K_elast(Q, Q));
U_elast(Q) = linear_system_solver.solve(K_elast(Q, Q), f(Q));
linear_system_solver.expand_deflation_basis(U_elast(Q));
omega_el = f(Q)' * U_elast(Q);

d_omega_ini = omega_el / 2;
U_elast = U_elast / 2;

[U, t_hist, omega_hist, U_max_hist] = CONTINUATION.LL_indirect_continuation(...
    d_omega_ini, d_t_min, step_max, LL_omega_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, U_elast, Q, f, ...
    constitutive_matrix_builder, linear_system_solver.copy());

time_run = toc;
fprintf('Running_time = %f \n', time_run);

if ~isempty(strfind(upper(char(solver_type)), 'BOOMERAMG'))
    LINEAR_SOLVERS.hypre_boomeramg_clear();
end

%% Profiler Summary

profiler.print_summary();

%% 6) Mechanical Results and Convergence

plotter.add_solution('indirect', U, t_hist, omega_hist, U_max_hist);
plotter.plot_displacements();
drawnow; pause(0.2);
plotter.plot_deviatoric_strain(0.75);
drawnow; pause(0.2);

figure; hold on; box on; grid on;
plot(omega_hist, t_hist, '-o');
title('Indirect continuation method for the LL method', 'Interpreter', 'latex');
xlabel('Control variable - $\omega$', 'Interpreter', 'latex');
ylabel('Load factor - $t$', 'Interpreter', 'latex');

fprintf('Notebook workflow completed.\n');
