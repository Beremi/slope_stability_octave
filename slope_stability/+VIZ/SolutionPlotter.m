classdef SolutionPlotter < handle
    %SOLUTIONPLOTTER  High-level 3-D FE result visualizer with multi-solution comparison.
    %
    %   plotter = VIZ.SolutionPlotter(coord, elem, surf, B, Xi)
    %
    %   The object stores the mesh once and accepts multiple registered solutions.
    %   Each plot method creates side-by-side subfigures when more than one solution
    %   is present, providing easy visual comparison of e.g. direct vs indirect SSR.
    %
    %   Typical workflow
    %   ----------------
    %     plotter = VIZ.SolutionPlotter(coord, elem, surf, B, Xi);
    %     plotter.set_pore_pressure(pw);
    %     plotter.add_solution('indirect', U3, lambda_hist3, omega_hist3, Umax_hist3);
    %     plotter.plot_mesh();
    %     plotter.plot_pore_pressure();
    %     plotter.plot_displacements();
    %     plotter.plot_deviatoric_strain();
    %     plotter.plot_deviatoric_slices({[], [35], [1e-16, 21.6506]}, 1);
    %     plotter.plot_convergence();
    %
    %   See individual methods for details.

    properties (SetAccess = private)
        coord           % (3 x n_n) nodal coordinates
        elem            % (n_p x n_e) element connectivity
        surf            % (n_v x n_s) boundary surface connectivity (auto-filtered)
        B               % strain-displacement matrix  (6*n_int x 3*n_n)
        Xi              % (3 x n_q)  quadrature local coordinates

        edges_merged    % precomputed merged bounding edges  (M x 2)
        faces_struct    % precomputed face-detection struct  (from VIZ.detect_faces)

        pore_pressure   % (1 x n_n) or []
        solutions       % struct array: .name .U .lambda_hist .omega_hist .Umax_hist
        n_solutions     % number of registered solutions
    end

    % =====================================================================
    %  PUBLIC INTERFACE
    % =====================================================================
    methods

        function obj = SolutionPlotter(coord, elem, surf, B, Xi)
            %SOLUTIONPLOTTER  Constructor.
            %
            %   plotter = VIZ.SolutionPlotter(coord, elem, surf)
            %   plotter = VIZ.SolutionPlotter(coord, elem, surf, B, Xi)
            %
            %   coord  (3 x n_n) - nodal coordinates
            %   elem   (n_p x n_e) - element connectivity
            %   surf   (n_v x n_s) - surface connectivity (filtered to boundary automatically)
            %   B      (optional)  - strain-displacement matrix (needed for deviatoric plots)
            %   Xi     (optional)  - quadrature local coordinates (needed for slice plots)
            obj.coord = coord;
            obj.elem  = elem;
            obj.surf  = VIZ.SolutionPlotter.filter_boundary_surface(coord, elem, surf);
            if nargin >= 4; obj.B  = B;  else; obj.B  = []; end
            if nargin >= 5; obj.Xi = Xi; else; obj.Xi = []; end
            obj.pore_pressure = [];
            obj.solutions   = struct('name',{}, 'U',{}, 'lambda_hist',{}, ...
                'omega_hist',{}, 'Umax_hist',{});
            obj.n_solutions = 0;

            % Precompute bounding edges and face grouping (used by many plots).
            obj.edges_merged = VIZ.compute_bounding_edges(obj.surf, obj.coord, 0);
            obj.faces_struct = VIZ.detect_faces(obj.edges_merged, obj.coord);
        end

        % -----------------------------------------------------------------
        function set_strain_data(obj, B, Xi)
            %SET_STRAIN_DATA  Set strain-displacement matrix and quadrature points.
            %   plotter.set_strain_data(B, Xi)
            %   Call after construction if B/Xi were not provided to the constructor.
            obj.B  = B;
            obj.Xi = Xi;
        end

        % -----------------------------------------------------------------
        function set_pore_pressure(obj, pw)
            %SET_PORE_PRESSURE  Register nodal pore-pressure field.
            %   plotter.set_pore_pressure(pw)   pw is (1 x n_n).
            obj.pore_pressure = pw;
        end

        % -----------------------------------------------------------------
        function add_solution(obj, name, U, lambda_hist, omega_hist, Umax_hist)
            %ADD_SOLUTION  Register a displacement solution and its convergence history.
            %
            %   plotter.add_solution(name, U, lambda_hist, omega_hist)
            %   plotter.add_solution(name, U, lambda_hist, omega_hist, Umax_hist)
            %
            %   name         - string label (e.g. 'direct', 'indirect')
            %   U            - (3 x n_n) displacement field
            %   lambda_hist  - vector of strength-reduction-factor history
            %   omega_hist   - vector of control-variable history
            %   Umax_hist    - (optional) vector of max-displacement history
            if nargin < 6; Umax_hist = []; end
            idx = obj.n_solutions + 1;
            obj.solutions(idx).name        = name;
            obj.solutions(idx).U           = U;
            obj.solutions(idx).lambda_hist = lambda_hist;
            obj.solutions(idx).omega_hist  = omega_hist;
            obj.solutions(idx).Umax_hist   = Umax_hist;
            obj.n_solutions = idx;
        end

        % -----------------------------------------------------------------
        function fig = plot_mesh(obj)
            %PLOT_MESH  Draw the 3-D surface mesh with bounding edges.
            fig = figure('Name', 'Mesh');
            ax  = gca;
            hold(ax, 'on');
            verts = obj.coord';
            patch(ax, 'Faces', obj.surf(1:3,:)', 'Vertices', verts, ...
                'FaceColor', [0.93 0.93 0.93], 'FaceAlpha', 1.0, ...
                'EdgeColor', [0.30 0.45 0.80], 'LineWidth', 0.35);
            axis(ax, 'equal'); axis(ax, 'tight');
            VIZ.apply_standard_3d_camera(ax, verts);
            obj.draw_bounding_edges(ax);
            set(ax, 'box', 'off');
            hold(ax, 'off');
            title(ax, 'Mesh');
            drawnow;
        end

        % -----------------------------------------------------------------
        function fig = plot_pore_pressure(obj)
            %PLOT_PORE_PRESSURE  Visualize registered pore-pressure on the boundary surface.
            assert(~isempty(obj.pore_pressure), ...
                'VIZ.SolutionPlotter: no pore pressure registered (call set_pore_pressure first).');
            fig = figure('Name', 'Pore Pressure');
            ax  = gca;
            obj.draw_quantity_on_axes(ax, obj.coord, obj.surf, ...
                zeros(size(obj.coord)), obj.pore_pressure);
            obj.draw_bounding_edges(ax);
            title(ax, 'Pore pressures [kPa]');
            drawnow;
        end

        % -----------------------------------------------------------------
        function fig = plot_displacements(obj, scale_factor)
            %PLOT_DISPLACEMENTS  Displacement magnitude on deformed shape.
            %
            %   fig = plotter.plot_displacements()
            %   fig = plotter.plot_displacements(scale_factor)     default 0.05
            %
            %   Multiple solutions produce side-by-side subplots.
            if nargin < 2; scale_factor = 0.05; end
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');
            N   = obj.n_solutions;
            fig = figure('Name', 'Displacements');
            for k = 1:N
                ax = subplot(1, N, k);
                U  = obj.solutions(k).U;
                obj.draw_displacement_on_axes(ax, U, scale_factor);
                title(ax, sprintf('Displacements: %s', obj.solutions(k).name));
            end
            drawnow;
        end

        % -----------------------------------------------------------------
        function fig = plot_deviatoric_strain(obj, clim_scale)
            %PLOT_DEVIATORIC_STRAIN  Deviatoric strain norm on the boundary surface.
            %
            %   fig = plotter.plot_deviatoric_strain()
            %   fig = plotter.plot_deviatoric_strain(clim_scale)   default 0.25
            %
            %   clim_scale shrinks the upper color limit to clim_scale * max.
            %   Multiple solutions produce side-by-side subplots.
            if nargin < 2; clim_scale = []; end
            assert(~isempty(obj.B), 'VIZ.SolutionPlotter: B matrix not provided.');
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');
            N   = obj.n_solutions;
            fig = figure('Name', 'Deviatoric Strain');
            for k = 1:N
                ax = subplot(1, N, k);
                U  = obj.solutions(k).U;
                obj.draw_deviatoric_on_axes(ax, U);
                if ~isempty(clim_scale)
                    cl = caxis(ax);
                    caxis(ax, [cl(1), clim_scale * cl(2)]);
                end
                title(ax, sprintf('Deviatoric strain: %s', obj.solutions(k).name));
            end
            drawnow;
        end

        % -----------------------------------------------------------------
        function [figs, slice_info] = plot_deviatoric_slices(obj, plane_vals, h, clim)
            %PLOT_DEVIATORIC_SLICES  Deviatoric strain norm on axis-aligned cutting planes.
            %
            %   [figs, info] = plotter.plot_deviatoric_slices(plane_vals, h)
            %   [figs, info] = plotter.plot_deviatoric_slices(plane_vals, h, clim)
            %
            %   plane_vals  1x3 cell {x_vals, y_vals, z_vals}  ([] to skip an axis)
            %   h           target edge length for slice triangulation (default 1)
            %   clim        [cmin cmax] or [] for auto
            %
            %   For each plane, one figure is produced.  If multiple solutions are
            %   registered, each figure contains side-by-side subplots.
            if nargin < 3 || isempty(h);    h    = 1;  end
            if nargin < 4;                  clim = []; end
            assert(~isempty(obj.B),  'VIZ.SolutionPlotter: B matrix not provided.');
            assert(~isempty(obj.Xi), 'VIZ.SolutionPlotter: Xi not provided.');
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');

            warnState = warning('off', 'all');

            N = obj.n_solutions;

            % Deviatoric projector (constant).
            iota = [1; 1; 1; 0; 0; 0];
            VOL  = iota * iota.';
            DEV  = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;

            % Integration-point coordinates (solution-independent).
            transformed_points = ASSEMBLY.integration_points(obj.elem, obj.coord, obj.Xi);

            % Per-solution deviatoric strain norm at IPs.
            norm_E_all = cell(N, 1);
            for k = 1:N
                E     = obj.B * obj.solutions(k).U(:);
                E     = reshape(E, 6, []);
                dev_E = DEV * E;
                norm_E_all{k} = sqrt(max(0, sum(E .* dev_E, 1))).';
            end

            % Collect requested (plane_id, plane_val) pairs.
            req = [];
            for pid = 1:3
                vals = plane_vals{pid};
                if ~isempty(vals)
                    req = [req; [pid * ones(numel(vals), 1), vals(:)]]; %#ok<AGROW>
                end
            end

            cmap_name = local_colormap_name();
            figs = [];
            slice_info = struct('plane_id',{}, 'plane_val',{}, 'TR2',{}, 'P3',{}, 'poly3',{});

            for q = 1:size(req, 1)
                plane_id  = req(q, 1);
                plane_val = req(q, 2);

                % Geometry: slice polygon & triangulation (shared across solutions).
                [poly3, ~] = VIZ.slice_by_plane(obj.faces_struct, obj.edges_merged, ...
                    obj.coord, plane_id, plane_val);
                if isempty(poly3); continue; end
                [TR2, P3] = VIZ.triangulate_polygon_slice(poly3, h, struct('plot2d', false));
                [ttl_base, xl, yl, free_axes] = local_plane_labels(plane_id, plane_val);

                fig = figure('Name', sprintf('Slice %s', ttl_base));
                for k = 1:N
                    ax = subplot(1, N, k);
                    norm_E = norm_E_all{k};
                    obj.draw_slice_on_axes(ax, TR2, poly3, transformed_points, ...
                        norm_E, plane_id, plane_val, free_axes, h, clim, cmap_name);
                    xlabel(ax, xl, 'Interpreter', 'latex', 'FontSize', 14);
                    ylabel(ax, yl, 'Interpreter', 'latex', 'FontSize', 14);
                    title(ax, sprintf('%s: %s', ttl_base, obj.solutions(k).name), ...
                        'Interpreter', 'latex', 'FontSize', 14);
                end
                drawnow;

                figs(end + 1, 1) = fig; %#ok<AGROW>
                slice_info(end + 1) = struct('plane_id', plane_id, 'plane_val', plane_val, ...
                    'TR2', TR2, 'P3', P3, 'poly3', poly3); %#ok<AGROW>
            end
            warning(warnState);
        end

        % -----------------------------------------------------------------
        function fig = plot_convergence(obj)
            %PLOT_CONVERGENCE  Overlay omega-lambda curves for all registered solutions.
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');
            fig = figure('Name', 'SSR Convergence');
            ax  = gca;
            hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
            markers = {'-o', '-s', '-d', '-^', '-v', '->'};
            for k = 1:obj.n_solutions
                mk = markers{mod(k - 1, numel(markers)) + 1};
                plot(ax, obj.solutions(k).omega_hist, obj.solutions(k).lambda_hist, ...
                    mk, 'DisplayName', obj.solutions(k).name);
            end
            xlabel(ax, 'control variable - $\omega$', 'Interpreter', 'latex');
            ylabel(ax, 'strength reduction factor - $\lambda$', 'Interpreter', 'latex');
            if obj.n_solutions > 1
                legend(ax, 'Location', 'best', 'Interpreter', 'latex');
            end
            title(ax, 'SSR Convergence', 'Interpreter', 'latex');
            hold(ax, 'off');
            drawnow;
        end

        % -----------------------------------------------------------------
        function save_figure(obj, filename, width_cm, height_cm) %#ok<INUSL>
            %SAVE_FIGURE  Save the current figure to a PNG file at 600 dpi.
            %
            %   plotter.save_figure('output.png', 20, 12)
            %
            %   width_cm, height_cm are in centimetres.
            fig = gcf;
            set(fig, 'Units', 'centimeters');
            set(fig, 'Position', [2 2 width_cm height_cm]);
            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperSize', [width_cm height_cm]);
            set(fig, 'PaperPositionMode', 'manual');
            set(fig, 'PaperPosition', [0 0 width_cm height_cm]);
            print(fig, filename, '-dpng', '-r600');
        end

    end % public methods

    % =====================================================================
    %  PRIVATE HELPERS
    % =====================================================================
    methods (Access = private)

        function draw_bounding_edges(obj, ax)
            %DRAW_BOUNDING_EDGES  Overlay precomputed feature edges on a given axes.
            edges = obj.edges_merged;
            n_e   = size(edges, 1);
            X = zeros(3, n_e);
            Y = zeros(3, n_e);
            Z = zeros(3, n_e);
            for i = 1:n_e
                X(:, i) = [obj.coord(1, edges(i,1)); obj.coord(1, edges(i,2)); NaN];
                Y(:, i) = [obj.coord(2, edges(i,1)); obj.coord(2, edges(i,2)); NaN];
                Z(:, i) = [obj.coord(3, edges(i,1)); obj.coord(3, edges(i,2)); NaN];
            end
            hold(ax, 'on');
            plot3(ax, X(:), Y(:), Z(:), 'k-', 'LineWidth', 1);
        end

        function draw_quantity_on_axes(obj, ax, coord_in, surf_in, U_in, Q_node)
            %DRAW_QUANTITY_ON_AXES  Patch plot of a scalar field on (possibly deformed) surface.
            hold(ax, 'on');
            verts = (coord_in + U_in)';
            if size(surf_in, 1) > 3
                tri_pattern = [1 4 6; 4 2 5; 4 5 6; 6 5 3];
                surf_all = [surf_in(tri_pattern(1,:), :)';
                    surf_in(tri_pattern(2,:), :)';
                    surf_in(tri_pattern(3,:), :)';
                    surf_in(tri_pattern(4,:), :)'];
                patch(ax, 'Faces', surf_all, 'Vertices', verts, ...
                    'FaceVertexCData', Q_node', 'FaceColor', 'interp', 'EdgeColor', 'none');
            else
                patch(ax, 'Faces', surf_in(1:3,:)', 'Vertices', verts, ...
                    'FaceVertexCData', Q_node', 'FaceColor', 'interp', 'EdgeColor', 'none');
            end
            colorbar(ax);
            set(ax, 'box', 'off');
            axis(ax, 'equal'); axis(ax, 'tight');
            VIZ.apply_standard_3d_camera(ax, verts);
            hold(ax, 'off');
        end

        function draw_displacement_on_axes(obj, ax, U, scale_factor)
            %DRAW_DISPLACEMENT_ON_AXES  Displacement magnitude on deformed mesh.
            U_norm = sqrt(sum(U.^2));
            scale  = scale_factor * max(abs(obj.coord(:))) / max(abs(U(:)));

            elem_values = U_norm(obj.elem);
            [bnd_faces, bnd_coord, bnd_vals] = ...
                VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, elem_values);

            U_x = U(1,:);  U_x_elem = U_x(obj.elem);
            U_y = U(2,:);  U_y_elem = U_y(obj.elem);
            U_z = U(3,:);  U_z_elem = U_z(obj.elem);
            [~, ~, Ux_face] = VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, U_x_elem);
            [~, ~, Uy_face] = VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, U_y_elem);
            [~, ~, Uz_face] = VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, U_z_elem);

            obj.draw_quantity_on_axes(ax, bnd_coord, bnd_faces, ...
                [Ux_face; Uy_face; Uz_face] * scale, bnd_vals);
            obj.draw_bounding_edges(ax);
        end

        function draw_deviatoric_on_axes(obj, ax, U)
            %DRAW_DEVIATORIC_ON_AXES  Deviatoric strain norm on boundary faces.
            elem_values = VIZ.get_elem_stress_3D(U, obj.B);
            [bnd_faces, bnd_coord, bnd_vals] = ...
                VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, elem_values);
            obj.draw_quantity_on_axes(ax, bnd_coord, bnd_faces, ...
                zeros(size(bnd_coord)), bnd_vals);
            obj.draw_bounding_edges(ax);
        end

        function draw_slice_on_axes(obj, ax, TR2, poly3, transformed_points, ...
                norm_E, plane_id, plane_val, free_axes, h, clim, cmap_name) %#ok<INUSL>
            %DRAW_SLICE_ON_AXES  Interpolate IP values onto a 2-D slice mesh and plot.
            CL = TR2.ConnectivityList;
            V2 = TR2.Points;

            % Sample integration points close to the slice plane.
            d = abs(transformed_points(plane_id, :) - plane_val);
            slice_tol = max(1e-6, 0.75 * h);
            idx = d <= slice_tol;
            if ~any(idx)
                [~, ord] = sort(d, 'ascend');
                take = min(10000, numel(ord));
                idx = false(size(d));
                idx(ord(1:take)) = true;
            end
            pts2 = transformed_points(free_axes, idx).';
            vals = norm_E(idx);
            if size(pts2, 1) > 30000
                rng(0);
                take = randperm(size(pts2, 1), 30000);
                pts2 = pts2(take, :);
                vals = vals(take);
            end

            % Interpolate onto triangulation nodes.
            vals_on_nodes = local_interp2d(pts2, vals, V2);

            % Plot.
            hold(ax, 'on');
            patch(ax, 'Faces', CL, 'Vertices', V2, ...
                'FaceVertexCData', double(vals_on_nodes), ...
                'FaceColor', 'interp', 'EdgeColor', 'none');
            poly2 = poly3(free_axes, :).';
            plot(ax, poly2(:,1), poly2(:,2), 'k-', 'LineWidth', 1.2);
            hold(ax, 'off');

            axis(ax, 'equal'); axis(ax, 'tight');
            grid(ax, 'on');  box(ax, 'on');
            if ~isempty(clim) && isfloat(clim) && numel(clim) == 2
                caxis(ax, clim);
            end
            colormap(ax, cmap_name);
            colorbar(ax);
        end

    end % private methods

    % =====================================================================
    %  STATIC HELPERS
    % =====================================================================
    methods (Static, Access = private)

        function surf_out = filter_boundary_surface(coord, elem, surf)
            %FILTER_BOUNDARY_SURFACE  Keep only surface faces that lie on the mesh boundary.
            v  = elem(1:4, :);
            F  = [v([2 3 4], :), v([1 3 4], :), v([1 2 4], :), v([1 2 3], :)];
            Fs = sort(F, 1).';
            [~, ~, ic] = unique(Fs, 'rows');
            cnt   = accumarray(ic, 1);
            bndFs = Fs(cnt(ic) == 1, :);

            surfV  = surf(1:3, :);
            surfKs = sort(surfV, 1).';
            isBnd  = ismember(surfKs, bndFs, 'rows');
            surf_out = surf(:, isBnd);
        end

    end % static methods

end % classdef

% =========================================================================
%  FILE-LOCAL HELPER FUNCTIONS
% =========================================================================

function [ttl, xl, yl, free_axes] = local_plane_labels(plane_id, plane_val)
switch plane_id
    case 1;  ttl = sprintf('$x = %.6g$', plane_val); xl = '$y$'; yl = '$z$'; free_axes = [2 3];
    case 2;  ttl = sprintf('$y = %.6g$', plane_val); xl = '$x$'; yl = '$z$'; free_axes = [1 3];
    case 3;  ttl = sprintf('$z = %.6g$', plane_val); xl = '$x$'; yl = '$y$'; free_axes = [1 2];
    otherwise; error('Invalid plane_id: %d', plane_id);
end
end

function cmap_name = local_colormap_name()
if exist('parula', 'file') == 2 || exist('parula', 'builtin') == 5
    cmap_name = 'parula';
else
    cmap_name = 'jet';
end
end

function vals_q = local_interp2d(pts2, vals, q2)
%LOCAL_INTERP2D  Robust 2-D scattered interpolation with nearest fallback.
if size(pts2, 1) < 3
    vals_q = repmat(mean(vals), size(q2, 1), 1);
    return;
end
vals_q = griddata(pts2(:,1), pts2(:,2), vals, q2(:,1), q2(:,2), 'linear');
bad = isnan(vals_q);
if any(bad)
    vals_q(bad) = griddata(pts2(:,1), pts2(:,2), vals, q2(bad,1), q2(bad,2), 'nearest');
end
if any(isnan(vals_q))
    vals_q(isnan(vals_q)) = mean(vals);
end
end
