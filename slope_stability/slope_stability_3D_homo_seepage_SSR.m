% Auto-generated from slope_stability_3D_homo_seepage_SSR.ipynb.

% Update the notebook first, then re-run scripts/sync_notebooks_to_m.py.

%% 3D Homogeneous Seepage SSR (HYPRE + Profiler)

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
mat_props = [6, 45, 0, 40000, 0.3, 20, 20];
k = 1.0;

x1 = 15; x2 = 10; x3 = 15;
y1 = 10; y2 = 10; z = 5;
N_h = 1;

%% 2) Reference Element Data and Mesh Load

[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);

switch(elem_type)
    case 'P1'
        [coord, elem, surf, Q] = MESH.mesh_P1_3D(N_h, x1, x2, x3, y1, y2, z);
    case 'P2'
        [coord, elem, surf, Q] = MESH.mesh_P2_3D(N_h, x1, x2, x3, y1, y2, z);
    otherwise
        error('Bad choice of element type');
end

[coord, elem, surf, Q] = MESH.reorder_mesh(coord, elem, surf, Q);

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

material_identifier = zeros(1, n_e);

%% Mesh Preview

plotter = VIZ.SolutionPlotter(coord, elem, surf);
plotter.plot_mesh();
drawnow;

%% 3) Seepage Problem

conduct0 = SEEPAGE.heter_conduct(material_identifier, n_q, k);
grho = 9.81;

Q_w = true(1, n_n);
Q_w(coord(1,:)<=0.001)=0;
Q_w(coord(1,:)>=x1+x2+x3-0.001)=0;
Q_w(coord(2,:)>=y1+y2-0.001)=0;
Q_w((coord(2,:)>=y1-0.001)&(coord(1,:)>=x1+x2-0.001))=0;
Q_w((coord(2,:)>=y1-0.001)&(coord(2,:)>=-(y2/x2)*coord(1,:)+y1+y2*(1+x1/x2)-0.001))=0;

y21=2;
y22=6;
pw_D=zeros(1,n_n);
x_bar=x1+(1-y21/y2)*x2;
part1=(coord(1,:)<x_bar)&(coord(2,:)<=-(y22/x_bar)*coord(1,:)+y1+y21+y22);
part2=coord(1,:)>=x_bar;
pw_D(part1)=grho*((y22/x_bar)*(x_bar-coord(1,part1))+y1+y21-coord(2,part1));
pw_D(part2)=grho*(y1+y21-coord(2,part2));

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

[pw, grad_p, mater_sat] = SEEPAGE.seepage_problem_3D(...
    coord, elem, Q_w, pw_D, grho, conduct0, HatP, DHatP1, DHatP2, DHatP3, WF, seepage_solver);

LINEAR_SOLVERS.hypre_boomeramg_clear();

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

[K_elast, B, WEIGHT, DPhi1_out, DPhi2_out, DPhi3_out] = ASSEMBLY.elastic_stiffness_matrix_3D(...
    elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

f_V_int = [-grad_p(1,:); -grad_p(2,:)-gamma; -grad_p(3,:)];
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

plotter.set_strain_data(B, Xi);

%% 5) Continuation, Newton, and Linear Solver Parameters

lambda_init = 0.7;
d_lambda_init = 0.1;
d_lambda_min = 1e-5;
d_lambda_diff_scaled_min = 0.001;
omega_max_stop = 7e7;
step_max = 100;

it_newt_max = 50;
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

%% 6) Run SSR Continuation (Indirect Only)

direct_on = 0;
indirect_on = 1;

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

if indirect_on
    plotter.add_solution('indirect', U3, lambda_hist3, omega_hist3, Umax_hist3);
end

plotter.plot_displacements();
drawnow; pause(0.2);
plotter.plot_deviatoric_strain(0.25);
drawnow; pause(0.2);
plotter.plot_convergence();

fprintf('Notebook workflow completed.\n');
