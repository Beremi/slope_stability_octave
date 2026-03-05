function fig = plot_deviatoric_strain_2D(U, coord, elem, B, clim_scale)
%PLOT_DEVIATORIC_STRAIN_2D Backward-compatible wrapper using VIZ.SolutionPlotter.

if nargin < 5
    clim_scale = [];
end

plotter = VIZ.SolutionPlotter(coord, elem, [], B, []);
plotter.add_solution('solution', U, [], [], []);
fig = plotter.plot_deviatoric_strain(clim_scale);
if numel(fig) > 1
    fig = fig(1);
end

drawnow;
end
