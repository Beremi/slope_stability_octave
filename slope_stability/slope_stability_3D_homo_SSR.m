% Auto-generated from slope_stability_3D_homo_SSR.ipynb.

% Update the notebook first, then re-run scripts/sync_notebooks_to_m.py.

%% 3D Homogeneous SSR (HYPRE + Profiler)

%% Notes Before Running
% - Activate local environment: `source setup/activate_optimized_octave.sh`
% - In Jupyter, select `Octave (local-rsb)` kernel.
% - This notebook follows the same workflow style as the COMSOL notebook, adapted to the homogeneous SSR benchmark.

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

% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';

% Mechanical parameters per subdomain:
% [c0, phi, psi, young, poisson, gamma_sat, gamma_unsat]
mat_props = [6, 45, 0, 40000, 0.3, 20, 20];

% Mesh file (benchmark family)
%   SSR_homo_uni.h5:    25415 / 15770
%   SSR_homo_ada_L1.h5: 25591 / 16615
%   SSR_homo_ada_L2.h5: 48449 / 32985
%   SSR_homo_ada_L3.h5: 91659 / 63968
%   SSR_homo_ada_L4.h5: 174902 / 124300
%   SSR_homo_ada_L5.h5: 336774 / 242716
file_path = 'meshes/SSR_homo_ada_L1.h5';

%% 2) Reference Element Data and Mesh Load

% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);

% local basis functions and their derivatives
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

switch(elem_type)
    case 'P1'
        error('Prepared meshes are only for P2 elements.');
    case 'P2'
        [coord, elem, surf, Q, ~] = MESH.load_mesh_P2(file_path);
        fprintf('P2 elements: homogeneous slope\n');
    otherwise
        error('Bad choice of element type.');
end

% Same default as the COMSOL-style notebook: apply reordering for solver performance.
[coord, elem, surf, Q] = MESH.reorder_mesh(coord, elem, surf, Q);

% Mesh statistics
n_n = size(coord,2);
n_unknown = length(coord(Q));
n_e = size(elem,2);
n_q = length(WF);
n_int = n_e * n_q;

fprintf('\nMesh data:');
fprintf('  number of nodes =%d ', n_n);
fprintf('  number of unknowns =%d ', n_unknown);
fprintf('  number of elements =%d ', n_e);
fprintf('  number of integration points =%d ', n_int);
fprintf('\n');

% Homogeneous identifier
material_identifier = zeros(1, n_e);

%% Mesh Preview

% Use explicit COMSOL camera preset for visual parity.
plotter = VIZ.SolutionPlotter(coord, elem, surf, 'comsol');
plotter.plot_mesh();
drawnow;

%% 3) Mechanical Material Fields and Assembly

fields = {'c0', 'phi', 'psi', 'young', 'poisson', 'gamma_sat', 'gamma_unsat'};
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% Homogeneous benchmark uses fully saturated body.
saturation = true(1, n_int);

[c0, phi, psi, shear, bulk, lame, gamma] = ...
    ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);

% Elastic stiffness matrix + element-level derivative data (DPhi1/2/3)
[K_elast, B, WEIGHT, DPhi1_out, DPhi2_out, DPhi3_out] = ASSEMBLY.elastic_stiffness_matrix_3D(
    elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Gravity-only volume force
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

% Attach strain data to the plotter for later visualisation.
plotter.set_strain_data(B, Xi);

%% 4) Continuation, Newton, and Linear Solver Parameters

% Continuation parameters
lambda_init = 0.9;
d_lambda_init = 0.1;
d_lambda_min = 1e-3;
d_lambda_diff_scaled_min = 0.001;
omega_max_stop = 3000;
step_max = 100;

% Newton parameters
it_newt_max = 100;
it_damp_max = 10;
tol = 1e-4;
r_min = 1e-4;

% Linear solver settings (HYPRE BoomerAMG + DFGMRES)
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

% Constitutive model object
n_strain = 6;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(
    B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, 3);

% Provide element geometry for fast element-level tangent assembly.
constitutive_matrix_builder.set_element_data(elem, DPhi1_out, DPhi2_out, DPhi3_out);

% Shared profiler
profiler = PROFILING.Profiler();
constitutive_matrix_builder.profiler = profiler;
linear_system_solver.profiler = profiler;

%% 5) Run SSR Continuation (Indirect Only)

direct_on = 0;
indirect_on = 1;

if direct_on
    fprintf('\n Direct continuation method\n');
    tic;
    [U2, lambda_hist2, omega_hist2, Umax_hist2] = CONTINUATION.SSR_direct_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf('Running_time = %f \n', time_run);
end

if indirect_on
    fprintf('\n Indirect continuation method\n');
    tic;
    [U3, lambda_hist3, omega_hist3, Umax_hist3, stats] = CONTINUATION.SSR_indirect_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        omega_max_stop, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf('Running_time = %f \n', time_run);
end

if ~isempty(strfind(upper(char(solver_type)), 'BOOMERAMG'))
    LINEAR_SOLVERS.hypre_boomeramg_clear();
end

%% Profiler Summary

profiler.print_summary();

%% 6) Mechanical Results and Convergence

% Register computed solutions.
if direct_on
    plotter.add_solution('direct', U2, lambda_hist2, omega_hist2, Umax_hist2);
end
if indirect_on
    plotter.add_solution('indirect', U3, lambda_hist3, omega_hist3, Umax_hist3);
end

% Mechanical field plots.
plotter.plot_displacements();
drawnow; pause(0.2);

plotter.plot_deviatoric_strain(0.25);
drawnow; pause(0.2);

plotter.plot_convergence();

fprintf('Notebook workflow completed.\n');
