% Auto-generated from slope_stability_3D_hetero_seepage_SSR.ipynb.

% Update the notebook first, then re-run scripts/sync_notebooks_to_m.py.

%% 3D Heterogeneous Seepage SSR (Concave Waterlevels Mesh)

%% Notes Before Running
% - Activate the local environment first: `source setup/activate_optimized_octave.sh`
% - In Jupyter, pick the `Octave (local-rsb)` kernel.
% - This notebook is an Octave-first workflow aligned with `slope_stability_3D_hetero_seepage_SSR.m` and uses the same module APIs as the COMSOL notebook style.

%plot -f png -r 600
% Load sparsersb for multithreaded sparse matrix-vector products
pkg load sparsersb;

% Configure toolkit for both GUI and headless Jupyter environments.
graphics_toolkit('qt');
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
mat_props = [15, 38,  0, 50000, 0.30, 22, 22;  % General foundation
    10, 35,  0, 50000, 0.30, 21, 21;  % Relatively weak foundation
    18, 32,  0, 20000, 0.33, 20, 20;  % General slope mass
    15, 30,  0, 10000, 0.33, 19, 19;  % Cover layer
    ];

% Hydraulic conductivity for each subdomain [m/s]
k = [1; 1; 1; 1];

%% 2) Reference Element Data and Mesh Load

% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);

% local basis functions and their derivatives
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

% Available slope_with_waterlevels_concave*.h5 meshes (nodes / elements):
%   slope_with_waterlevels_concave.h5:    106520 / 72673
%   slope_with_waterlevels_concave_L2.h5:  69733 / 48205
%   slope_with_waterlevels_concave_L3.h5: 244609 / 174745
%   slope_with_waterlevels_concave_L4.h5: 473420 / 339843
file_path = 'meshes/slope_with_waterlevels_concave_L2.h5';
[coord, elem, surf, Q, material_identifier, triangle_labels] = ...
    MESH.load_mesh_gmsh_waterlevels(file_path);

% Uncomment to reduce matrix bandwidth via Reverse Cuthill-McKee node reordering:
[coord, elem, surf, Q] = MESH.reorder_mesh(coord, elem, surf, Q);

% Mesh statistics
n_n = size(coord,2);
n_unknown = length(coord(Q));
n_e = size(elem,2);
n_q = length(WF);
n_int = n_e * n_q;

fprintf('\n');
fprintf('Mesh data:');
fprintf('  number of nodes =%d ', n_n);
fprintf('  number of unknowns =%d ', n_unknown);
fprintf('  number of elements =%d ', n_e);
fprintf('  number of integration points =%d ', n_int);
fprintf('\n');

%% Mesh Preview

% Create the plotter early so we can visualise intermediate results.
plotter = VIZ.SolutionPlotter(coord, elem, surf);
plotter.plot_mesh();
drawnow;

%% 3) Seepage Problem

% Hydraulic conductivity at each integration point
conduct0 = SEEPAGE.heter_conduct(material_identifier, n_q, k);

% Specific weight of water [kPa]
grho = 9.81;

% Dirichlet boundary conditions for seepage (problem dependent)
[Q_w, pw_D] = MESH.seepage_boundary_3D_hetero(coord, surf, triangle_labels, grho);

% Set up iterative solver for seepage (scalar BoomerAMG, no nullspace)
seepage_boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
    'use_as_preconditioner', true);
seepage_solver_type = 'DFGMRES_HYPRE_BOOMERAMG';
seepage_linear_solver_tolerance = 1e-6;
seepage_linear_solver_maxit = 200;
seepage_deflation_basis_tolerance = 1e-3;
seepage_linear_solver_printing = 0;
seepage_solver = LINEAR_SOLVERS.set_linear_solver(seepage_solver_type, ...
    seepage_linear_solver_tolerance, seepage_linear_solver_maxit, ...
    seepage_deflation_basis_tolerance, seepage_linear_solver_printing, ...
    Q_w, [], seepage_boomeramg_opts);

% Solve seepage (DFGMRES + HYPRE inside Newton iterations)
[pw, grad_p, mater_sat] = SEEPAGE.seepage_problem_3D( ...
    coord, elem, Q_w, pw_D, grho, conduct0, HatP, DHatP1, DHatP2, DHatP3, WF, seepage_solver);

% Clean up HYPRE instance used for seepage
LINEAR_SOLVERS.hypre_boomeramg_clear();

% Integration-point saturation mask
mater_sat_ext = repmat(mater_sat, n_q, 1);
saturation = mater_sat_ext(:);

%% Pore Pressure Preview

plotter.set_pore_pressure(pw);
plotter.plot_pore_pressure();
drawnow;

%% 4) Mechanical Material Fields and Assembly

fields = {'c0', 'phi', 'psi', 'young', 'poisson', 'gamma_sat', 'gamma_unsat'};
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

[c0, phi, psi, shear, bulk, lame, gamma] = ...
    ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);

% Elastic stiffness matrix + element-level derivative data (DPhi1/2/3)
[K_elast, B, WEIGHT, DPhi1_out, DPhi2_out, DPhi3_out] = ASSEMBLY.elastic_stiffness_matrix_3D(
    elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

% Volume forces at integration points: seepage gradient + gravity
f_V_int = [-grad_p(1,:); -grad_p(2,:)-gamma; -grad_p(3,:)];
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

% Attach strain data to the plotter for later slice visualisation.
plotter.set_strain_data(B, Xi);

%% 5) Continuation, Newton, and Linear Solver Parameters

% Continuation parameters
lambda_init = 0.7;
d_lambda_init = 0.1;
d_lambda_min = 1e-5;
d_lambda_diff_scaled_min = 0.01
omega_max_stop = 2e7;           % Maximum omega, then stop (notebook runtime)
step_max = 100;                 % Maximum number of continuation steps (notebook runtime)

% Newton parameters
it_newt_max = 50;
it_damp_max = 10;
tol = 1e-4;
r_min = 1e-4;

% Linear solver settings
% agmg folder is baked into LINEAR_SOLVERS.set_linear_solver
solver_type = 'DFGMRES_HYPRE_BOOMERAMG';

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0

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

% Attach shared profiler to both objects (class-level instrumentation)
profiler = PROFILING.Profiler();
constitutive_matrix_builder.profiler = profiler;
linear_system_solver.profiler = profiler;

%% 6) Run SSR Continuation

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

%% 7) Mechanical Results and Convergence

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

plotter.plot_deviatoric_strain(0.75);
drawnow; pause(0.2);

plane_vals = {[], [35, 40, 55], [21.6506, 43.3013]};
plotter.plot_deviatoric_slices(plane_vals, 1);
drawnow; pause(0.2);

plotter.plot_convergence();

fprintf('Notebook workflow completed.\n');
