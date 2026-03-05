function fig = plot_displacements_2D(U, coord, elem, scale_factor)
%PLOT_DISPLACEMENTS_2D Backward-compatible wrapper using VIZ.SolutionPlotter.

if nargin < 4 || isempty(scale_factor)
    scale_factor = 0.05;
end

plotter = VIZ.SolutionPlotter(coord, elem, []);
plotter.add_solution('solution', U, [], [], []);
fig = plotter.plot_displacements(scale_factor);
if numel(fig) > 1
    fig = fig(1);
end

drawnow;
end
