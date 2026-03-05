% Auto-generated from slope_stability_2D_homo_LL.ipynb.

% Update the notebook first, then re-run scripts/sync_notebooks_to_m.py.

%% Homogeneous slope and its stability (via LL method)

%% Notes Before Running
% - Activate local environment: `source setup/activate_optimized_octave.sh`
% - Activate Jupyter venv: `source .venv/bin/activate`
% - In Jupyter, select kernel **Octave (local-rsb)**

close all; clearvars; clc;
%plot -f png -r 600
pkg load sparsersb;

try
    graphics_toolkit('qt');
catch
    graphics_toolkit('gnuplot');
end
set(0, 'defaultfigurevisible', 'off');

fprintf('Working directory: %s\n', pwd);
fprintf('Graphics toolkit: %s\n', graphics_toolkit());

%% Homogeneous slope and its stability (via LL method)

% =========================================================================
%
%  This program solves a 2D slope stability problem by the limit
%  load (LL) method described in (Sysala et al., CAS 2025). The Mohr-
%  Coulomb yield criterion, Davis approach, standard finite elements
%  (either P1 or P2 elements) and meshes with different densities are
%  considered. For P2 elements, the 7-point Gauss quadrature
%  is used. To find the safety factor of the LL method, the indirect
%  continuation technique is used. A benchmark with a homogeneous slope
%  is considered. It is possible to change slope inclination and other
%  geometrical parameters.
%
% =========================================================================

%% Main input data

% elem_type - type of finite elements; available choices: 'P1', 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A', 'B', 'C'
Davis_type = 'B';

% Material parameters for each subdomain. In the following table, we
% specify in each column the following material parameters, respectively:
% [c0, phi, psi, young, poisson, gamma_sat, gamma_unsat], where
%    c0 ... Cohesion (c)
%    phi ... Friction angle (phi in degrees)
%    psi ... Dilatancy angle (psi in degrees)
%    young ... Young's modulus (E)
%    poisson ...  Poisson's ratio (nu)
%    gamma_sat ...   Specific weight - saturated (gamma_sat in kN/m^3)
%    gamma_unsat ... Specific weight - unsaturated (gamma_unsat in kN/m^3)
% If gamma_sat and gamma_unsat are not distinguished, use the same values
% for these parameters. Each row of the table represents one subdomain. If
% a homogeneous body is considered, only one row is prescribed.
mat_props = [6, 45, 0, 40000, 0.3, 20, 20];

% Geometrical parameters
x1 = 15;         % Length of the body in front of the slope
x3 = 15;         % Length of the body behind the slope
y1 = 10;         % Height of the body below the slope
y2 = 10;         % Height of the slope
beta = pi/4;     % Slope angle
x2 = y2/tan(beta); % Length of the slope in the x-direction

% Mesh data
h = 1/2;         % Discretization parameter

%% Data from the reference element

% Quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_2D(elem_type);
% Local basis functions and their derivatives
[HatP, DHatP1, DHatP2] = ASSEMBLY.local_basis_volume_2D(elem_type, Xi);

%% Creation of the uniform finite element mesh

switch(elem_type)
    case 'P1'
        [coord, elem, ELEM_ED, EDGE_EL, Q] = MESH.mesh_P1_2D(h, x1, x2, x3, y1, y2);
        fprintf('P1 elements: \n')
    case 'P2'
        [coord, elem, ELEM_ED, EDGE_EL, Q] = MESH.mesh_P2_2D(h, x1, x2, x3, y1, y2);
        fprintf('P2 elements: \n')
    otherwise
        error('Bad choice of element type');
end

% Uncomment to reduce matrix bandwidth via Reverse Cuthill-McKee node reordering:
% [coord, elem, ~, Q] = MESH.reorder_mesh(coord, elem, [], Q);

% Number of nodes, elements, and integration points + print
n_n = size(coord,2);          % Number of nodes
n_unknown = length(coord(Q)); % Number of unknowns
n_e = size(elem,2);           % Number of elements
n_ed = size(EDGE_EL,2);       % Number of edges
n_q = length(WF);             % Number of quadrature points
n_int = n_e * n_q;            % Total number of integration points

fprintf('\n The mesh data:');
fprintf('  Number of nodes = %d ', n_n);
fprintf('  Number of unknowns = %d ', n_unknown);
fprintf('  Number of elements = %d ', n_e);
fprintf('  Number of edges = %d ', n_ed);
fprintf('  Number of integration points = %d \n', n_int);

% The array material_identifier for a homogeneous body
material_identifier = zeros(1,n_e);

%% Material parameters at integration points

% Fields with prescribed material properties
fields = {'c0',      ... % Cohesion (c)
    'phi',     ... % Friction angle (phi in degrees)
    'psi',     ... % Dilatancy angle (psi in degrees)
    'young',   ... % Young's modulus (E)
    'poisson', ... % Poisson's ratio (nu)
    'gamma_sat', ... % Specific weight - saturated (gamma_sat in kN/m^3)
    'gamma_unsat'};  % Specific weight - unsaturated (gamma_unsat in kN/m^3)

% Convert properties to structured format.
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% saturation - a prescribed logical array indicating integration points
%              where the body is saturated. If gamma_sat and gamma_unsat
%              are the same, set saturation=true(1,n_int). Otherwise,
%              this logical array is derived from a given phreatic surface.
saturation = true(1,n_int);

% Material parameters at integration points.
[c0, phi, psi, shear, bulk, lame, gamma] = ...
    ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);

%% Assembling of the elastic stiffness matrix

[K_elast, B, WEIGHT, DPhi1_out, DPhi2_out] = ASSEMBLY.elastic_stiffness_matrix_2D(elem, coord, DHatP1, DHatP2, WF, shear, lame);

%% Assembling of the vector of volume forces

% Volume forces at integration points, size(f_V_int) = (2, n_int)
f_V_int = [zeros(1, n_int); -gamma];
% Vector of volume forces
f_V = ASSEMBLY.vector_volume_2D(elem, coord, f_V_int, HatP, WEIGHT);

%% Input parameters for the indirect continuation method

d_t_min = 1e-3;                % Minimal increment of load factor t.
step_max = 100;                % Maximum number of continuation steps.
LL_omega_max = 2000;            % Maximum value of the control parameter omega.

%% Input parameters for Newton's solvers

it_newt_max = 100;               % Number of Newton's iterations
it_damp_max = 10;               % Number of iterations within line search
tol = 1e-4;                     % Relative tolerance for Newton's solvers
r_min = 1e-4;                   % Basic minimal regularization of the stiffness matrix

%% Defining linear solver

% agmg folder is baked into LINEAR_SOLVERS.set_linear_solver
solver_type = 'DFGMRES_HYPRE_BOOMERAMG'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG", "DFGMRES_HYPRE_BOOMERAMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

% Optional BoomerAMG options (used when solver_type contains BOOMERAMG).
boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
    'use_as_preconditioner', true);

[linear_system_solver] = LINEAR_SOLVERS.set_linear_solver(solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing, Q, coord, boomeramg_opts);

%% Constitutive problem and matrix builder

dim = 2;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);
% Enable element-level tangent assembly path for 2D B*D*B values.
constitutive_matrix_builder.set_element_data(elem, DPhi1_out, DPhi2_out);
disp(['2D element-level tangent mex enabled = ' , num2str(constitutive_matrix_builder.elem_use_mex)]);
disp(['2D constitutive mex enabled = ' , num2str(constitutive_matrix_builder.use_2D_mex)]);


%--------------------------------------------------------------------------

%% Profiler Setup

profiler = PROFILING.Profiler();
constitutive_matrix_builder.profiler = profiler;
linear_system_solver.profiler = profiler;

%% Computation of the limit load factor by the indirect continuation

fprintf('\n Indirect continuation method for the LL method\n');
tic;

% Compute the elastic displacement field as a starting point.
constitutive_matrix_builder.reduction(1.0);
U_elast = zeros(2, n_n);              % Elastic displacement vector.
linear_system_solver.setup_preconditioner(K_elast(Q, Q));
U_elast(Q) = linear_system_solver.solve(K_elast(Q, Q), f_V(Q));
linear_system_solver.expand_deflation_basis(U_elast(Q));
omega_el = f_V(Q)' * U_elast(Q);        % Work of external forces.

% Set the initial increment of omega.
d_omega_ini = omega_el / 5;
% Scale the elastic displacement field.
U_elast = U_elast / 5;

% Run the indirect continuation method for the LL method.
[U, t_hist, omega_hist, U_max_hist] = CONTINUATION.LL_indirect_continuation(...
    d_omega_ini, d_t_min, step_max, LL_omega_max, ...
    it_newt_max, it_damp_max, tol, r_min, K_elast, U_elast, Q, f_V, ...
    constitutive_matrix_builder, linear_system_solver.copy());

time_run = toc;
fprintf("Running_time = %f \n", time_run);


if ~isempty(strfind(upper(char(solver_type)), 'BOOMERAMG'))
    LINEAR_SOLVERS.hypre_boomeramg_clear();
end

%% Postprocessing - visualization of selected results

if false
VIZ.plot_deviatoric_strain_2D(U,coord,elem,B);
VIZ.plot_displacements_2D(U,coord,elem);
% Visualization of the continuation curve: omega -> lambda.
figure; hold on; box on; grid on;
plot(omega_hist, t_hist, '-o');
title('Indirect continuation method for the LL method', 'Interpreter', 'latex');
xlabel('Control variable - $\omega$', 'Interpreter', 'latex');
ylabel('Load factor - $t$', 'Interpreter', 'latex');

end

%% Profiler Summary

if exist('profiler', 'var')
    profiler.print_summary();
end

%% Unified Postprocessing via SolutionPlotter

plotter = VIZ.SolutionPlotter(coord, elem, [], B, [], 'comsol');

if exist('material_identifier', 'var')
    plotter.set_material_identifier(material_identifier);
end
if exist('mater_sat', 'var')
    plotter.set_element_saturation(mater_sat);
end
if exist('pw', 'var')
    plotter.set_pore_pressure(pw);
end

if exist('U2', 'var') && exist('lambda_hist2', 'var') && exist('omega_hist2', 'var')
    if exist('Umax_hist2', 'var')
        Umax2 = Umax_hist2;
    else
        Umax2 = [];
    end
    plotter.add_solution('direct', U2, lambda_hist2, omega_hist2, Umax2, struct( ...
        'title', 'Direct continuation method', ...
        'xlabel', 'Control variable - $\omega$', ...
        'ylabel', 'strength reduction factor - $\lambda$', ...
        'marker', '-o'));
end

if exist('U3', 'var') && exist('lambda_hist3', 'var') && exist('omega_hist3', 'var')
    if exist('Umax_hist3', 'var')
        Umax3 = Umax_hist3;
    else
        Umax3 = [];
    end
    plotter.add_solution('indirect', U3, lambda_hist3, omega_hist3, Umax3, struct( ...
        'title', 'Indirect continuation method', ...
        'xlabel', 'Control variable - $\omega$', ...
        'ylabel', 'strength reduction factor - $\lambda$', ...
        'marker', '-o'));
end

if exist('U', 'var') && exist('t_hist', 'var') && exist('omega_hist', 'var')
    if exist('U_max_hist', 'var')
        UmaxLL = U_max_hist;
    else
        UmaxLL = [];
    end
    plotter.add_solution('indirect_LL', U, t_hist, omega_hist, UmaxLL, struct( ...
        'title', 'Indirect continuation method for the LL method', ...
        'xlabel', 'Control variable - $\omega$', ...
        'ylabel', 'Load factor - $t$', ...
        'marker', '-o'));
end

if plotter.n_solutions > 0
    plotter.plot_deviatoric_strain();
    plotter.plot_displacements();
    plotter.plot_convergence();
end

if exist('pw', 'var')
    plotter.plot_pore_pressure();
end
if exist('material_identifier', 'var')
    plotter.plot_material_map();
end
if exist('mater_sat', 'var')
    plotter.plot_saturation();
end
