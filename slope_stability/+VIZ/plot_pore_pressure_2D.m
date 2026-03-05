function fig = plot_pore_pressure_2D(pw, coord, elem)
%PLOT_PORE_PRESSURE_2D Backward-compatible wrapper using VIZ.SolutionPlotter.

plotter = VIZ.SolutionPlotter(coord, elem, []);
plotter.set_pore_pressure(pw);
fig = plotter.plot_pore_pressure();

drawnow;
end
