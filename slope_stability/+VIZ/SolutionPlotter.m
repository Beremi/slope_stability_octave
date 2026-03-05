classdef SolutionPlotter < handle
    %SOLUTIONPLOTTER High-level FE result visualizer for 2D/3D solutions.
    %
    %   Supported mesh types:
    %     - 3D: coord is (3 x n_n), elem is tetrahedral connectivity.
    %     - 2D: coord is (2 x n_n), elem is triangular connectivity.
    %
    %   Constructor forms:
    %     plotter = VIZ.SolutionPlotter(coord, elem, surf)
    %     plotter = VIZ.SolutionPlotter(coord, elem, surf, B, Xi)
    %     plotter = VIZ.SolutionPlotter(coord, elem, surf, 'comsol')
    %     plotter = VIZ.SolutionPlotter(coord, elem, surf, B, Xi, 'comsol')
    %
    %   For 2D meshes, pass surf as []:
    %     plotter = VIZ.SolutionPlotter(coord2d, elem2d, [], B)
    %
    %   Camera presets for 3D:
    %     'comsol'      - legacy camera orientation.
    %     'mirrored_x'  - mirrored orientation across x-center plane (default).

    properties (SetAccess = private)
        mesh_dim        % 2 or 3
        coord           % (dim x n_n) nodal coordinates
        elem            % (n_p x n_e) element connectivity
        surf            % (n_v x n_s) 3D boundary surface connectivity or [] for 2D
        B               % strain-displacement matrix
        Xi              % quadrature local coordinates

        edges_merged    % precomputed 3D merged boundary edges
        faces_struct    % precomputed 3D face-detection struct

        pore_pressure   % nodal pore-pressure field or []
        material_identifier % element-level material ids or []
        element_saturation  % element-level saturation flags or []
        solutions       % struct array with registered solutions
        n_solutions     % number of registered solutions
        camera_view_tag % camera preset tag for 3D views
    end

    methods

        function obj = SolutionPlotter(coord, elem, surf, varargin)
            if nargin < 3
                surf = [];
            end

            [B, Xi, camera_tag] = local_parse_solutionplotter_args(varargin{:});

            obj.coord = coord;
            obj.elem = elem;
            obj.mesh_dim = size(coord, 1);
            assert(obj.mesh_dim == 2 || obj.mesh_dim == 3, ...
                'VIZ.SolutionPlotter: coord must have 2 or 3 rows.');

            if obj.mesh_dim == 3
                if isempty(surf)
                    obj.surf = VIZ.SolutionPlotter.infer_boundary_surface_from_tetra(elem);
                else
                    obj.surf = VIZ.SolutionPlotter.filter_boundary_surface(coord, elem, surf);
                end
            else
                obj.surf = [];
            end

            obj.B = B;
            obj.Xi = Xi;
            obj.pore_pressure = [];
            obj.material_identifier = [];
            obj.element_saturation = [];
            obj.solutions = struct('name', {}, 'U', {}, 'lambda_hist', {}, ...
                'omega_hist', {}, 'Umax_hist', {}, 'curve_meta', {});
            obj.n_solutions = 0;
            obj.camera_view_tag = local_normalize_camera_tag(camera_tag);

            if obj.mesh_dim == 3 && ~isempty(obj.surf)
                obj.edges_merged = VIZ.compute_bounding_edges(obj.surf, obj.coord, 0);
                obj.faces_struct = VIZ.detect_faces(obj.edges_merged, obj.coord);
            else
                obj.edges_merged = [];
                obj.faces_struct = [];
            end
        end

        function set_strain_data(obj, B, Xi)
            obj.B = B;
            obj.Xi = Xi;
        end

        function set_pore_pressure(obj, pw)
            obj.pore_pressure = pw;
        end

        function set_material_identifier(obj, material_identifier)
            obj.material_identifier = material_identifier;
        end

        function set_element_saturation(obj, element_saturation)
            obj.element_saturation = element_saturation;
        end

        function add_solution(obj, name, U, lambda_hist, omega_hist, Umax_hist, curve_meta)
            if nargin < 6 || isempty(Umax_hist)
                Umax_hist = [];
            end
            if nargin < 7 || isempty(curve_meta)
                curve_meta = struct();
            end
            curve_meta = local_merge_curve_meta_defaults(curve_meta, name);

            idx = obj.n_solutions + 1;
            obj.solutions(idx).name = name;
            obj.solutions(idx).U = U;
            obj.solutions(idx).lambda_hist = lambda_hist;
            obj.solutions(idx).omega_hist = omega_hist;
            obj.solutions(idx).Umax_hist = Umax_hist;
            obj.solutions(idx).curve_meta = curve_meta;
            obj.n_solutions = idx;
        end

        function fig = plot_mesh(obj)
            fig = obj.new_figure('Mesh');
            ax = gca;

            if obj.mesh_dim == 3
                hold(ax, 'on');
                verts = obj.coord';
                patch(ax, 'Faces', obj.surf(1:3, :)', 'Vertices', verts, ...
                    'FaceColor', [0.93 0.93 0.93], 'FaceAlpha', 1.0, ...
                    'EdgeColor', [0.30 0.45 0.80], 'LineWidth', 0.35);
                axis(ax, 'equal'); axis(ax, 'tight');
                VIZ.apply_standard_3d_camera(ax, verts, obj.camera_opts());
                obj.draw_bounding_edges(ax);
                set(ax, 'box', 'off');
                hold(ax, 'off');
            else
                hold(ax, 'on');
                verts = obj.coord';
                faces = obj.triangulation_faces_2d(obj.elem);
                patch(ax, 'Faces', faces, 'Vertices', verts, ...
                    'FaceColor', [0.93 0.93 0.93], 'FaceAlpha', 1.0, ...
                    'EdgeColor', [0.30 0.45 0.80], 'LineWidth', 0.35);
                axis(ax, 'equal'); axis(ax, 'tight');
                grid(ax, 'on');
                set(ax, 'box', 'off');
                xlabel(ax, 'x', 'Interpreter', 'none');
                ylabel(ax, 'y', 'Interpreter', 'none');
                obj.draw_bounding_edges(ax, faces, verts);
                hold(ax, 'off');
            end

            title(ax, 'Mesh');
            drawnow;
        end

        function fig = plot_pore_pressure(obj)
            assert(~isempty(obj.pore_pressure), ...
                'VIZ.SolutionPlotter: no pore pressure registered (call set_pore_pressure first).');

            fig = obj.new_figure('Pore Pressure');
            ax = gca;

            if obj.mesh_dim == 3
                conn = obj.surf;
            else
                conn = obj.elem;
            end

            obj.draw_quantity_on_axes(ax, obj.coord, conn, ...
                zeros(size(obj.coord)), obj.pore_pressure);

            if obj.mesh_dim == 3
                obj.draw_bounding_edges(ax);
            end
            title(ax, 'Pore pressures [kPa]');
            drawnow;
        end

        function fig = plot_material_map(obj)
            assert(~isempty(obj.material_identifier), ...
                'VIZ.SolutionPlotter: no material_identifier registered (call set_material_identifier first).');
            fig = obj.new_figure('Material map');
            ax = gca;
            obj.draw_element_quantity_on_axes(ax, double(obj.material_identifier));
            title(ax, 'Material identifier');
            drawnow;
        end

        function fig = plot_saturation(obj)
            assert(~isempty(obj.element_saturation), ...
                'VIZ.SolutionPlotter: no element saturation registered (call set_element_saturation first).');
            fig = obj.new_figure('Saturation');
            ax = gca;
            obj.draw_element_quantity_on_axes(ax, double(obj.element_saturation));
            caxis(ax, [0 1]);
            title(ax, 'Saturation');
            drawnow;
        end

        function fig = plot_displacements(obj, scale_factor)
            if nargin < 2 || isempty(scale_factor)
                scale_factor = 0.05;
            end
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');

            N = obj.n_solutions;
            figs = zeros(N, 1);
            for k = 1:N
                sol_name = obj.solutions(k).name;
                figs(k) = obj.new_figure(sprintf('Displacements: %s', sol_name));
                ax = gca;
                obj.draw_displacement_on_axes(ax, obj.solutions(k).U, scale_factor);
                title(ax, sprintf('Displacements: %s', sol_name), 'Interpreter', 'none');
                drawnow;
            end

            if N == 1
                fig = figs(1);
            else
                fig = figs;
            end
        end

        function fig = plot_deviatoric_strain(obj, clim_scale)
            if nargin < 2
                clim_scale = [];
            end
            assert(~isempty(obj.B), 'VIZ.SolutionPlotter: B matrix not provided.');
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');

            N = obj.n_solutions;
            figs = zeros(N, 1);
            for k = 1:N
                sol_name = obj.solutions(k).name;
                figs(k) = obj.new_figure(sprintf('Deviatoric strain: %s', sol_name));
                ax = gca;
                obj.draw_deviatoric_on_axes(ax, obj.solutions(k).U);
                if ~isempty(clim_scale)
                    cl = caxis(ax);
                    caxis(ax, [cl(1), clim_scale * cl(2)]);
                end
                title(ax, sprintf('Deviatoric strain: %s', sol_name), 'Interpreter', 'none');
                drawnow;
            end

            if N == 1
                fig = figs(1);
            else
                fig = figs;
            end
        end

        function [figs, slice_info] = plot_deviatoric_slices(obj, plane_vals, h, clim)
            if nargin < 3 || isempty(h)
                h = 1;
            end
            if nargin < 4
                clim = [];
            end

            assert(obj.mesh_dim == 3, ...
                'VIZ.SolutionPlotter: plot_deviatoric_slices is available only for 3D meshes.');
            assert(~isempty(obj.B), 'VIZ.SolutionPlotter: B matrix not provided.');
            assert(~isempty(obj.Xi), 'VIZ.SolutionPlotter: Xi not provided.');
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');

            warnState = warning('off', 'all');

            N = obj.n_solutions;

            iota = [1; 1; 1; 0; 0; 0];
            VOL = iota * iota.';
            DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;

            transformed_points = ASSEMBLY.integration_points(obj.elem, obj.coord, obj.Xi);

            norm_E_all = cell(N, 1);
            for k = 1:N
                E = obj.B * obj.solutions(k).U(:);
                E = reshape(E, 6, []);
                dev_E = DEV * E;
                norm_E_all{k} = sqrt(max(0, sum(E .* dev_E, 1))).';
            end

            req = [];
            for pid = 1:3
                vals = plane_vals{pid};
                if ~isempty(vals)
                    req = [req; [pid * ones(numel(vals), 1), vals(:)]]; %#ok<AGROW>
                end
            end

            cmap_name = local_colormap_name();
            figs = [];
            slice_info = struct('plane_id', {}, 'plane_val', {}, 'TR2', {}, ...
                'P3', {}, 'poly3', {}, 'solution_index', {}, 'solution_name', {});

            for q = 1:size(req, 1)
                plane_id = req(q, 1);
                plane_val = req(q, 2);

                [poly3, ~] = VIZ.slice_by_plane(obj.faces_struct, obj.edges_merged, ...
                    obj.coord, plane_id, plane_val);
                if isempty(poly3)
                    continue;
                end

                [TR2, P3] = VIZ.triangulate_polygon_slice(poly3, h, struct('plot2d', false));
                [ttl_base, xl, yl, free_axes] = local_plane_labels(plane_id, plane_val);

                for k = 1:N
                    sol_name = obj.solutions(k).name;
                    fig_k = obj.new_figure(sprintf('Slice %s: %s', ttl_base, sol_name));
                    ax = gca;
                    norm_E = norm_E_all{k};
                    obj.draw_slice_on_axes(ax, TR2, poly3, transformed_points, ...
                        norm_E, plane_id, plane_val, free_axes, h, clim, cmap_name);
                    xlabel(ax, xl, 'Interpreter', 'latex', 'FontSize', 14);
                    ylabel(ax, yl, 'Interpreter', 'latex', 'FontSize', 14);
                    title(ax, sprintf('%s: %s', ttl_base, sol_name), ...
                        'Interpreter', 'none', 'FontSize', 14);
                    drawnow;

                    figs(end + 1, 1) = fig_k; %#ok<AGROW>
                    slice_info(end + 1) = struct('plane_id', plane_id, ...
                        'plane_val', plane_val, 'TR2', TR2, 'P3', P3, ...
                        'poly3', poly3, 'solution_index', k, ...
                        'solution_name', sol_name); %#ok<AGROW>
                end
            end

            warning(warnState);
        end

        function fig = plot_convergence(obj)
            assert(obj.n_solutions > 0, 'VIZ.SolutionPlotter: no solutions registered.');

            markers = {'-o', '-s', '-d', '-^', '-v', '->'};
            N = obj.n_solutions;
            figs = zeros(N, 1);

            for k = 1:N
                sol = obj.solutions(k);
                figs(k) = obj.new_figure(sprintf('Continuation: %s', sol.name));
                ax = gca;
                hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

                mk = sol.curve_meta.marker;
                if isempty(mk)
                    mk = markers{mod(k - 1, numel(markers)) + 1};
                end

                plot(ax, sol.omega_hist, sol.lambda_hist, mk, 'LineWidth', 1.2);
                xlabel(ax, sol.curve_meta.xlabel, 'Interpreter', 'latex');
                ylabel(ax, sol.curve_meta.ylabel, 'Interpreter', 'latex');
                title(ax, sol.curve_meta.title, 'Interpreter', 'latex');
                hold(ax, 'off');
                drawnow;
            end

            if N == 1
                fig = figs(1);
            else
                fig = figs;
            end
        end

        function save_figure(obj, filename, width_cm, height_cm) %#ok<INUSL>
            fig = gcf;
            set(fig, 'Units', 'centimeters');
            set(fig, 'Position', [2 2 width_cm height_cm]);
            set(fig, 'PaperUnits', 'centimeters');
            set(fig, 'PaperSize', [width_cm height_cm]);
            set(fig, 'PaperPositionMode', 'manual');
            set(fig, 'PaperPosition', [0 0 width_cm height_cm]);
            print(fig, filename, '-dpng', '-r600');
        end

    end

    methods (Access = private)

        function fig = new_figure(obj, name) %#ok<INUSD>
            fig = figure('Name', name, 'Color', 'w');
            set(fig, 'Units', 'pixels');
            pos = get(fig, 'Position');
            if numel(pos) == 4
                pos(3) = 900;
                pos(4) = 700;
                set(fig, 'Position', pos);
            end
        end

        function draw_bounding_edges(obj, ax, faces2d, verts2d)
            if obj.mesh_dim == 3
                edges = obj.edges_merged;
                if isempty(edges)
                    return;
                end

                n_e = size(edges, 1);
                X = zeros(3, n_e);
                Y = zeros(3, n_e);
                Z = zeros(3, n_e);
                for i = 1:n_e
                    X(:, i) = [obj.coord(1, edges(i, 1)); obj.coord(1, edges(i, 2)); NaN];
                    Y(:, i) = [obj.coord(2, edges(i, 1)); obj.coord(2, edges(i, 2)); NaN];
                    Z(:, i) = [obj.coord(3, edges(i, 1)); obj.coord(3, edges(i, 2)); NaN];
                end
                hold(ax, 'on');
                plot3(ax, X(:), Y(:), Z(:), 'k-', 'LineWidth', 1);
                return;
            end

            if nargin < 3 || isempty(faces2d)
                faces2d = obj.triangulation_faces_2d(obj.elem);
            end
            if nargin < 4 || isempty(verts2d)
                verts2d = obj.coord';
            end
            edges = obj.boundary_edges_from_faces_2d(faces2d);
            if isempty(edges)
                return;
            end

            n_e = size(edges, 1);
            X = zeros(3, n_e);
            Y = zeros(3, n_e);
            for i = 1:n_e
                X(:, i) = [verts2d(edges(i, 1), 1); verts2d(edges(i, 2), 1); NaN];
                Y(:, i) = [verts2d(edges(i, 1), 2); verts2d(edges(i, 2), 2); NaN];
            end
            hold(ax, 'on');
            plot(ax, X(:), Y(:), 'k-', 'LineWidth', 1.0);
        end

        function draw_quantity_on_axes(obj, ax, coord_in, surf_in, U_in, Q_node, draw_outline_2d)
            if nargin < 7 || isempty(draw_outline_2d)
                draw_outline_2d = true;
            end
            hold(ax, 'on');

            if obj.mesh_dim == 3
                verts = (coord_in + U_in)';
                if size(surf_in, 1) > 3
                    tri_pattern = [1 4 6; 4 2 5; 4 5 6; 6 5 3];
                    surf_all = [surf_in(tri_pattern(1, :), :)'; ...
                        surf_in(tri_pattern(2, :), :)'; ...
                        surf_in(tri_pattern(3, :), :)'; ...
                        surf_in(tri_pattern(4, :), :)'];
                    patch(ax, 'Faces', surf_all, 'Vertices', verts, ...
                        'FaceVertexCData', Q_node(:), ...
                        'FaceColor', 'interp', 'EdgeColor', 'none');
                else
                    patch(ax, 'Faces', surf_in(1:3, :)', 'Vertices', verts, ...
                        'FaceVertexCData', Q_node(:), ...
                        'FaceColor', 'interp', 'EdgeColor', 'none');
                end
                colorbar(ax);
                set(ax, 'box', 'off');
                axis(ax, 'equal'); axis(ax, 'tight');
                VIZ.apply_standard_3d_camera(ax, verts, obj.camera_opts());
                hold(ax, 'off');
                return;
            end

            verts = (coord_in + U_in)';
            faces = obj.triangulation_faces_2d(surf_in);
            patch(ax, 'Faces', faces, 'Vertices', verts, ...
                'FaceVertexCData', Q_node(:), ...
                'FaceColor', 'interp', 'EdgeColor', 'none');

            colorbar(ax);
            set(ax, 'box', 'off');
            axis(ax, 'equal'); axis(ax, 'tight');
            grid(ax, 'on');
            xlabel(ax, 'x', 'Interpreter', 'none');
            ylabel(ax, 'y', 'Interpreter', 'none');
            if draw_outline_2d
                obj.draw_bounding_edges(ax, faces, verts);
            end
            hold(ax, 'off');
        end

        function opts = camera_opts(obj)
            opts = local_camera_opts_from_tag(obj.camera_view_tag);
        end

        function draw_displacement_on_axes(obj, ax, U, scale_factor)
            if obj.mesh_dim == 3
                U_norm = sqrt(sum(U.^2));
                u_max = max(abs(U(:)));
                c_max = max(abs(obj.coord(:)));
                if u_max > 0 && c_max > 0
                    scale = scale_factor * c_max / u_max;
                else
                    scale = 0;
                end

                elem_values = U_norm(obj.elem);
                [bnd_faces, bnd_coord, bnd_vals] = ...
                    VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, elem_values);

                U_x = U(1, :); U_x_elem = U_x(obj.elem);
                U_y = U(2, :); U_y_elem = U_y(obj.elem);
                U_z = U(3, :); U_z_elem = U_z(obj.elem);
                [~, ~, Ux_face] = VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, U_x_elem);
                [~, ~, Uy_face] = VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, U_y_elem);
                [~, ~, Uz_face] = VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, U_z_elem);

                obj.draw_quantity_on_axes(ax, bnd_coord, bnd_faces, ...
                    [Ux_face; Uy_face; Uz_face] * scale, bnd_vals);
                obj.draw_bounding_edges(ax);
                return;
            end

            U_norm = sqrt(sum(U.^2));
            u_max = max(abs(U(:)));
            c_max = max(abs(obj.coord(:)));
            if u_max > 0 && c_max > 0
                scale = scale_factor * c_max / u_max;
            else
                scale = 0;
            end

            % Plot deformed field without auto-outline first.
            obj.draw_quantity_on_axes(ax, obj.coord, obj.elem, U * scale, U_norm, false);

            % Use original mesh boundary as requested.
            faces_orig = obj.triangulation_faces_2d(obj.elem);
            verts_orig = obj.coord';
            obj.draw_bounding_edges(ax, faces_orig, verts_orig);

            % Keep 2D plot frame anchored to original mesh extents.
            mins = min(obj.coord, [], 2);
            maxs = max(obj.coord, [], 2);
            xlim(ax, [mins(1), maxs(1)]);
            ylim(ax, [mins(2), maxs(2)]);
            axis(ax, 'equal');
        end

        function draw_deviatoric_on_axes(obj, ax, U)
            if obj.mesh_dim == 3
                elem_values = VIZ.get_elem_stress_3D(U, obj.B);
                [bnd_faces, bnd_coord, bnd_vals] = ...
                    VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, elem_values);
                obj.draw_quantity_on_axes(ax, bnd_coord, bnd_faces, ...
                    zeros(size(bnd_coord)), bnd_vals);
                obj.draw_bounding_edges(ax);
                return;
            end

            elem_values = VIZ.get_elem_stress_2D(U, obj.B, obj.elem);
            [coord_new, elem_new, val_new] = obj.expand_element_values_2d(elem_values);
            % For 2D deviatoric strain, disable black outline overlay.
            obj.draw_quantity_on_axes(ax, coord_new, elem_new, zeros(size(coord_new)), val_new, false);
        end

        function draw_element_quantity_on_axes(obj, ax, elem_values)
            assert(numel(elem_values) == size(obj.elem, 2), ...
                'VIZ.SolutionPlotter: elem_values must have one value per element.');

            if obj.mesh_dim == 3
                [bnd_faces, bnd_coord, bnd_vals] = ...
                    VIZ.boundary_stress_from_elements_3D(obj.coord, obj.elem, elem_values);
                obj.draw_quantity_on_axes(ax, bnd_coord, bnd_faces, ...
                    zeros(size(bnd_coord)), bnd_vals);
                obj.draw_bounding_edges(ax);
                return;
            end

            [coord_new, elem_new, val_new] = obj.expand_element_values_2d(elem_values(:).');
            obj.draw_quantity_on_axes(ax, coord_new, elem_new, zeros(size(coord_new)), val_new, false);
        end

        function draw_slice_on_axes(obj, ax, TR2, poly3, transformed_points, ...
                norm_E, plane_id, plane_val, free_axes, h, clim, cmap_name) %#ok<INUSL>
            CL = TR2.ConnectivityList;
            V2 = TR2.Points;

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

            vals_on_nodes = local_interp2d(pts2, vals, V2);

            hold(ax, 'on');
            patch(ax, 'Faces', CL, 'Vertices', V2, ...
                'FaceVertexCData', double(vals_on_nodes), ...
                'FaceColor', 'interp', 'EdgeColor', 'none');
            poly2 = poly3(free_axes, :).';
            plot(ax, poly2(:, 1), poly2(:, 2), 'k-', 'LineWidth', 1.2);
            hold(ax, 'off');

            axis(ax, 'equal'); axis(ax, 'tight');
            grid(ax, 'on'); box(ax, 'on');
            if ~isempty(clim) && isfloat(clim) && numel(clim) == 2
                caxis(ax, clim);
            end
            colormap(ax, cmap_name);
            colorbar(ax);
        end

        function faces = triangulation_faces_2d(obj, elem_in) %#ok<INUSD>
            if isempty(elem_in)
                faces = zeros(0, 3);
                return;
            end

            if size(elem_in, 1) > 6
                tri_pattern = [1 6 4; 4 5 2; 3 5 6; 6 4 5];
                faces = [elem_in(tri_pattern(1, :), :)'; ...
                    elem_in(tri_pattern(2, :), :)'; ...
                    elem_in(tri_pattern(3, :), :)'; ...
                    elem_in(tri_pattern(4, :), :)'];
            elseif size(elem_in, 1) > 3
                tri_pattern = [1 6 5; 4 5 6; 2 4 6; 3 4 5];
                faces = [elem_in(tri_pattern(1, :), :)'; ...
                    elem_in(tri_pattern(2, :), :)'; ...
                    elem_in(tri_pattern(3, :), :)'; ...
                    elem_in(tri_pattern(4, :), :)'];
            else
                faces = elem_in(1:3, :)';
            end
        end

        function edges = boundary_edges_from_faces_2d(obj, faces) %#ok<INUSD>
            if isempty(faces)
                edges = zeros(0, 2);
                return;
            end

            edges_all = [faces(:, [1 2]); faces(:, [2 3]); faces(:, [3 1])];
            edges_sorted = sort(edges_all, 2);
            [uniq_edges, ~, ic] = unique(edges_sorted, 'rows');
            counts = accumarray(ic, 1);
            edges = uniq_edges(counts == 1, :);
        end

        function [coord_new, elem_new, val_new] = expand_element_values_2d(obj, elem_values)
            elem_work = obj.elem;
            if size(elem_work, 1) == 15
                elem_base = elem_work(1:6, :);
                elem_new = [reshape(1:numel(elem_base), size(elem_base)); zeros(9, size(elem_base, 2))];
            else
                elem_base = elem_work;
                elem_new = reshape(1:numel(elem_base), size(elem_base));
            end

            n_local = size(elem_base, 1);
            n_elem = size(elem_base, 2);

            x = obj.coord(1, :);
            y = obj.coord(2, :);
            x_elem = reshape(x(elem_base), 1, []);
            y_elem = reshape(y(elem_base), 1, []);
            coord_new = [x_elem; y_elem];

            % Accept both:
            %   (1 x n_elem) / (n_elem x 1): one value per element
            %   (n_local x n_elem): one value per local plotting node
            if isvector(elem_values)
                v = elem_values(:);
                if numel(v) == n_elem
                    val_mat = repmat(v.', n_local, 1);
                elseif numel(v) == n_local * n_elem
                    val_mat = reshape(v, n_local, n_elem);
                else
                    error(['VIZ.SolutionPlotter: elem_values size mismatch. ', ...
                        'Expected n_elem or n_local*n_elem values.']);
                end
            else
                [r, c] = size(elem_values);
                if r == n_local && c == n_elem
                    val_mat = elem_values;
                elseif r == n_elem && c == n_local
                    val_mat = elem_values.';
                elseif r == 1 && c == n_elem
                    val_mat = repmat(elem_values, n_local, 1);
                elseif c == 1 && r == n_elem
                    val_mat = repmat(elem_values.', n_local, 1);
                elseif numel(elem_values) == n_local * n_elem
                    val_mat = reshape(elem_values, n_local, n_elem);
                else
                    error(['VIZ.SolutionPlotter: elem_values size mismatch. ', ...
                        'Expected [n_local x n_elem] or one value per element.']);
                end
            end

            val_new = reshape(val_mat, [], 1);
        end

    end

    methods (Static, Access = private)

        function surf_out = filter_boundary_surface(coord, elem, surf) %#ok<INUSD>
            if isempty(surf)
                surf_out = VIZ.SolutionPlotter.infer_boundary_surface_from_tetra(elem);
                return;
            end

            v = elem(1:4, :);
            F = [v([2 3 4], :), v([1 3 4], :), v([1 2 4], :), v([1 2 3], :)];
            Fs = sort(F, 1).';
            [~, ~, ic] = unique(Fs, 'rows');
            cnt = accumarray(ic, 1);
            bndFs = Fs(cnt(ic) == 1, :);

            surfV = surf(1:3, :);
            surfKs = sort(surfV, 1).';
            isBnd = ismember(surfKs, bndFs, 'rows');
            surf_out = surf(:, isBnd);
        end

        function surf_out = infer_boundary_surface_from_tetra(elem)
            assert(size(elem, 1) >= 4, ...
                'VIZ.SolutionPlotter: cannot infer 3D boundary surface from element order < 4.');

            v = elem(1:4, :);
            F = [v([2 3 4], :), v([1 3 4], :), v([1 2 4], :), v([1 2 3], :)];
            Fs = sort(F, 1).';
            [uniq_faces, ~, ic] = unique(Fs, 'rows');
            cnt = accumarray(ic, 1);
            bndFs = uniq_faces(cnt == 1, :);
            surf_out = bndFs.';
        end

    end

end

function [ttl, xl, yl, free_axes] = local_plane_labels(plane_id, plane_val)
switch plane_id
    case 1
        ttl = sprintf('$x = %.6g$', plane_val); xl = '$y$'; yl = '$z$'; free_axes = [2 3];
    case 2
        ttl = sprintf('$y = %.6g$', plane_val); xl = '$x$'; yl = '$z$'; free_axes = [1 3];
    case 3
        ttl = sprintf('$z = %.6g$', plane_val); xl = '$x$'; yl = '$y$'; free_axes = [1 2];
    otherwise
        error('Invalid plane_id: %d', plane_id);
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
if size(pts2, 1) < 3
    vals_q = repmat(mean(vals), size(q2, 1), 1);
    return;
end
vals_q = griddata(pts2(:, 1), pts2(:, 2), vals, q2(:, 1), q2(:, 2), 'linear');
bad = isnan(vals_q);
if any(bad)
    vals_q(bad) = griddata(pts2(:, 1), pts2(:, 2), vals, q2(bad, 1), q2(bad, 2), 'nearest');
end
if any(isnan(vals_q))
    vals_q(isnan(vals_q)) = mean(vals);
end
end

function [B, Xi, camera_tag] = local_parse_solutionplotter_args(varargin)
B = [];
Xi = [];
camera_tag = 'mirrored_x';

switch numel(varargin)
    case 0
        return;
    case 1
        if local_is_text(varargin{1})
            camera_tag = local_to_char(varargin{1});
        else
            B = varargin{1};
        end
    case 2
        if local_is_text(varargin{2})
            B = varargin{1};
            camera_tag = local_to_char(varargin{2});
        else
            B = varargin{1};
            Xi = varargin{2};
        end
    case 3
        B = varargin{1};
        Xi = varargin{2};
        camera_tag = local_to_char(varargin{3});
    otherwise
        error(['VIZ.SolutionPlotter: bad optional argument list. ', ...
            'Use (B, Xi, camera_tag), (B, Xi), or (camera_tag).']);
end
end

function tag = local_normalize_camera_tag(tag_in)
tag = lower(strtrim(local_to_char(tag_in)));
switch tag
    case {'', 'default', 'mirror_x', 'mirrored_x', 'other_side_x', 'x_mirror'}
        tag = 'mirrored_x';
    case {'comsol'}
        tag = 'comsol';
    otherwise
        error('Unsupported camera_view_tag "%s". Supported: "mirrored_x", "comsol".', tag);
end
end

function opts = local_camera_opts_from_tag(tag)
switch local_normalize_camera_tag(tag)
    case 'comsol'
        opts = struct('view_dir', [-1.2, 1, -2]);
    case 'mirrored_x'
        opts = struct('view_dir', [1.2, 1, -2]);
    otherwise
        error('Unknown camera tag "%s".', tag);
end
end

function tf = local_is_text(v)
tf = ischar(v) || (isstring(v) && isscalar(v));
end

function out = local_to_char(v)
if ischar(v)
    out = v;
elseif isstring(v) && isscalar(v)
    out = char(v);
else
    error('Camera tag must be a char vector or scalar string.');
end
end

function meta = local_default_curve_meta(name)
meta = struct();
meta.title = sprintf('Continuation: %s', name);
meta.xlabel = 'Control variable - $\omega$';
meta.ylabel = 'Continuation parameter';
meta.marker = '';
end

function meta = local_merge_curve_meta_defaults(meta_in, name)
meta = local_default_curve_meta(name);
if isempty(meta_in)
    return;
end
if ~isstruct(meta_in)
    error('curve_meta must be a struct when provided.');
end
if isfield(meta_in, 'title') && ~isempty(meta_in.title)
    meta.title = meta_in.title;
end
if isfield(meta_in, 'xlabel') && ~isempty(meta_in.xlabel)
    meta.xlabel = meta_in.xlabel;
end
if isfield(meta_in, 'ylabel') && ~isempty(meta_in.ylabel)
    meta.ylabel = meta_in.ylabel;
end
if isfield(meta_in, 'marker') && ~isempty(meta_in.marker)
    meta.marker = meta_in.marker;
end
end
