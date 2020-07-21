import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# Axis labels
axis_labels = dict()
axis_labels["T_av"] = r"T ( K )"
axis_labels["T0"] = r"$T_0$ ( K )"
axis_labels["Tp"] = axis_labels["T_av"]
axis_labels["Tm"] = axis_labels["T_av"]
axis_labels["kxx"] = r"$\kappa_{\rm xx}$ ( W / K m )"
axis_labels["dTx"] = r"$\Delta T_{\rm x}$ ( K )"
axis_labels["dTy"] = r"$\Delta T_{\rm y}$ ( K )"
axis_labels["kxx/T"] = r"$\kappa_{\rm xx}$/T ( W / K$^2$ m )"
axis_labels["dTx/T"] = r"$\Delta T_{\rm x}$/T ( % )"
axis_labels["Resistance"] = r"(T-T$_0$)/$\Delta T_{\rm x}$"
axis_labels["kxy"] = r"$\kappa_{\rm xy}$ ( mW / K cm )"
axis_labels["kxy/kxx"] = r"$\kappa_{\rm xy}/\kappa_{\rm xx}$ ( % )"
axis_labels["dTy/dTx"] = r"$\Delta T_{\rm y}/\Delta T_{\rm x}$ ( % )"
axis_labels["Tp_Tm"] = axis_labels["T_av"]
axis_labels["T0_fit"] = axis_labels["T0"]
axis_labels["I_fit"] = r"I ( mA )"

# Legend labels
legend_labels = dict()
legend_labels["H"] = r"H = %sT"
legend_labels["sample"] = r"Sample: %s"
legend_labels["date"] = r"%s"


def Plot(xdata, ydata, xkey, ykey, *args, **kwargs):
    """
    Plots data corresponding to key.

    Parameters:
    ----------------------------------------------------------------------------
    xdata:  1d array
            The x-axis data to plot
    ydata:  1d array
            The y-axis data to plot

    xkey:    string
            Used to format x-axis label

    ykey:   string
            Used to format y-axis label

    kwargs are all passed to ax.plot from matplotlib except the following:
    ----------------------------------------------------------------------------
    show :  Bool
            Determines if the figure is shown ore closed defaults to True
    parameters : dict
            parameter:value used for legends

    axis_fs:    Int
            The axis labels fontsize

    fig:    matplotlib.figure
            Used to draw on an existing figure, requires ax

    ax:     matplotlib.ax
            Used to draw on an existing figure, requires fig
    """

    # Looks for known kwargs
    if "show" in kwargs:
        show = kwargs["show"]
        kwargs.pop("show")
        if type(show) is bool:
            pass
        else:
            if type(show) is None:
                pass
            else:
                raise TypeError("show must be bool or None")
    else:
        show = True

    if "axis_fs" in kwargs:
        axis_fs = kwargs["axis_fs"]
        kwargs.pop("axis_fs")
    else:
        axis_fs = 16

    if "parameters" in kwargs:
        parameters = kwargs["parameters"]
        kwargs.pop("parameters")
        if type(parameters) is dict:
            labels = []
            for key in parameters:
                if type(parameters[key]) is str:
                    if key in axis_labels:
                        labels.append(axis_labels[key] % (parameters[key]))
                    else:
                        labels.append(parameters[key])
                else:
                    raise TypeError("parameters must be strings")
            label = " ,".join(labels)
            label_size = len(labels)
            if label_size <= 1:
                label_font = axis_fs
            elif label_size == 2:
                label_font = axis_fs-2
            else:
                label_font = 10
        else:
            raise TypeError("parameters must be dict containing strings")

    else:
        label_size = 0
        label = ""

    if "fig" in kwargs and "ax" in kwargs:
        fig, ax = kwargs["fig"], kwargs["ax"]
        kwargs.pop("fig")
        kwargs.pop("ax")
        return_fig = False
    else:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        return_fig = True

    # Actual plotting
    ax.plot(xdata, ydata, label=label, *args, **kwargs)

    if ydata.min()*ydata.max() < 0:
        ax.axhline(0, "--k", lw=2)
    else:
        if ydata.min() >= 0:
            ax.set_ylim(0, ax.get_ylim()[1])
        else:
            ax.set_ylim(ax.get_ylim()[0], 0)

    if xkey in ["T_av", "T0"]:
        ax.set_xlim(0)
    else:
        pass

    # Make it pretty
    ax.set_xlabel(axis_labels[xkey], fontsize=axis_fs)
    ax.set_ylabel(axis_labels[ykey], fontsize=axis_fs)
    ax.tick_params(axis="both", which="both",
                   direction="in", top=True, right=True)

    if label_size != 0:
        ax.legend(label_font)
    else:
        pass

    if len(fig.axes) == 1:
        fig.tight_layout(rect=[0.01, 0.01, 1, 0.95])
    else:
        pass

    if show is True:
        plt.show()
    elif show is False:
        plt.close()
    else:
        pass

    return fig, ax
