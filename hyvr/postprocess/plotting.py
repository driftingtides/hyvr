import sys
import numpy as np
import matplotlib.pyplot as plt

def cross_section_pcolor(model, field, y=None, x=None, log=False, xlim_y=None, ylim_y=None,
        xlim_x=None, ylim_x=None, cmap=None):
    """
    Create a pcolor plot of a cross section of the specified field.

    Parameters
    ----------
    model : Model object
    field : str
        The field to be shown, e.g. 'ae', 'fac', 'hat', 'k_iso', etc.
    y : sequence of values or None, optional (default: None)
        If a sequence of values or a single value is given, cross section at this y-value will be
        plotted.
        If it is None, a y-cross section through the center of the model domain will be plotted.
        If it is an empty list, no y-cross section will be plotted.
    x : sequence of values or None, optional (default: None)
        If a sequence of values or a single value is given, cross section at this x-value will be
        plotted.
        If it is None, a x-cross section through the center of the model domain will be plotted.
        If it is an empty list, no x-cross section will be plotted.
    log : bool, optional (default: False)
        Whether to plot the logarithm of the field.
    xlim_y, ylim_y : lists/tuples, optional (default: None)
        Axis limits for the y-cross-sections
    xlim_x, ylim_x : lists/tuples, optional (default: None)
        Axis limits for the x-cross-sections
    cmap : matplotlib colormap, optional (default: None)
        Colormap to use for the plot, e.g. 'prism'.
    """
    if y is None:
        y = [model.grid.y0 + model.grid.ly/2]
    if x is None:
        x = [model.grid.x0 + model.grid.lx/2]
    if np.isscalar(y):
        y = [y]
    if np.isscalar(x):
        x = [x]

    for yi in y:
        fig, ax = plt.subplots()
        y_index = int(np.round((yi-(model.grid.y0+model.grid.dy/2)/model.grid.dy)))
        if 0 <= y_index and y_index < model.grid.ny:
            im = ax.pcolor(model.grid.X[:,y_index,:].T, model.grid.Z[:,y_index,:].T,
                    model.data[field][:,y_index,:].T, cmap=cmap)
            ax.set_title(field + ', y = {:.2f}'.format(yi))
            ax.set_xlabel('x')
            ax.set_ylabel('z')
            if xlim_y is not None: ax.set_xlim(xlim)
            if ylim_y is not None: ax.set_ylim(ylim)
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax)
            plt.show(fig)
        else:
            print("Warning: y = {:.2f} not in model domain".format(yi), file=sys.stderr)

    for xi in x:
        fig, ax = plt.subplots()
        x_index = int(np.round((xi-(model.grid.x0+model.grid.dx/2)/model.grid.dx)))
        if 0 <= x_index and x_index < model.grid.nx:
            im = ax.pcolor(model.grid.Y[x_index,:,:].T, model.grid.Z[x_index,:,:].T,
                    model.data[field][x_index,:,:].T, cmap=cmap)
            ax.set_title(field + ', x = {:.2f}'.format(xi))
            ax.set_xlabel('y')
            ax.set_ylabel('z')
            if xlim_x is not None: ax.set_xlim(xlim)
            if ylim_x is not None: ax.set_ylim(ylim)
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax)
            plt.show(fig)
        else:
            print("Warning: x = {:.2f} not in model domain".format(xi), file=sys.stderr)




