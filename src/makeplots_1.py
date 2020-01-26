"""
A Python script to make beautiful plots
"""

#%% Imports
from contextlib import suppress
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib
import tikzplotlib


#%% Import ipython and adjust for qt5 backend, keeping compatibility with
# non-ipython envs
with suppress(ImportError):
    from IPython import get_ipython
with suppress(AttributeError, NameError):
    # List available APIs
    get_ipython().run_line_magic("matplotlib", "-l")
    get_ipython().run_line_magic("matplotlib", "qt5")

print(plt.get_backend())
#%% Matplotlib adjustments
# Adjust matplotlib qt5 backend for high DPI monitors
if plt.get_backend() == "Qt5Agg":
    from matplotlib.backends.qt_compat import QtWidgets  # pylint: disable=C0412
    from sys import argv

    qApp = QtWidgets.QApplication(argv)
    plt.matplotlib.rcParams["figure.dpi"] = qApp.desktop().physicalDpiX()

# Set LaTeX compatibility settings
# matplotlib.use("pgf")
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "sans-serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
    }
)

# plt.rcParams["figure.autolayout"] = True

# Adjust figure default size
plt.rcParams["figure.figsize"] = (6, 4)  # Default fig size (w, h))

# Get the default color scheme
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

awgn_xtickspacing = 5

# plt.ioff()

#%% Get AWGN results from files
# Load the file
results_awgn = loadmat("../results/results1.mat")

#%% Get variables
print(list(results_awgn.keys()))

mustry = results_awgn["mustry"][0]
nstry = results_awgn["nstry"][0]

ds = results_awgn["ds"]
ds[~np.isfinite(ds)] = 200
ds = np.clip(ds, -100, 100)
W = results_awgn["W"]
W[~np.isfinite(W)] = 200
W = np.clip(W, -100, 100)
e = results_awgn["e"]
e[~np.isfinite(e)] = 200
# e = np.clip(e, -100, 100)
ys = results_awgn["ys"]
ys[~np.isfinite(ys)] = 200
ys = np.clip(ys, -100, 100)


#%% Show beautiful figures
lw = 1
for k in range(3):
    for i, n in enumerate(nstry):
        for j, mu in enumerate(mustry):
            num = k * len(nstry) * len(mustry) + i * len(mustry) + j
            fig, axs = plt.subplots(num=num, nrows=3, ncols=1, sharex=True)
            axs[0].set_title(f"Ex. 1, signal {k+1}, $n={n}$, $\\mu = {mu}$")
            # Stuff to plot
            print(num)
            d = ds[k, 0]
            y = ys[k, i, j, 0]
            diff = d - y
            ws = W[k, i, j]
            edb = 10 * np.log10(e[k, i, j])
            # Plot signals
            axs[0].plot(d, label="original", linewidth=lw)
            axs[0].plot(y, label="recovered", linewidth=lw)
            axs[0].plot(diff, label="difference", linewidth=lw)
            axs[0].legend(framealpha=1, fontsize=6, loc="lower right")
            axs[0].set_ylabel("Signal level")
            axs[0].yaxis.set_minor_locator(AutoMinorLocator())
            axs[0].grid(b=True, which="major", color="silver", linestyle="-")
            axs[0].grid(b=True, which="minor", color="lightgray", linestyle="--")
            # Plot the weights
            for v in range(n):
                axs[1].plot(ws[:, v], linewidth=lw)
            axs[1].set_ylabel("Weight value")
            axs[1].yaxis.set_minor_locator(AutoMinorLocator())
            axs[1].grid(b=True, which="major", color="silver", linestyle="-")
            axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
            # Plot the error (dB)
            axs[2].plot(edb, linewidth=lw)

            # Do adjustments
            # Like proper axis settings
            if k == 0:  # First signal
                axs[1].set_ylim(bottom=-2, top=2)
                axs[0].set_ylim(bottom=-10, top=10)
                if n == 3:
                    axs[2].set_ylim(bottom=-15, top=3)
                    axs[1].set_ylim(bottom=-5, top=5)
                if n == 5:
                    axs[2].set_ylim(bottom=-30, top=-3)
                if n == 10:
                    axs[2].set_ylim(bottom=-45, top=0)
            if k == 1:  # Second signal
                axs[1].set_ylim(bottom=-5, top=5)
                axs[0].set_ylim(bottom=-10, top=10)
                if n == 3:
                    axs[2].set_ylim(bottom=-43, top=5)
                if n == 5:
                    axs[2].set_ylim(bottom=-43, top=0)
                if n == 10:
                    axs[2].set_ylim(bottom=-45, top=0)
            if k == 2:  # Third signal
                axs[2].set_ylim(bottom=-40, top=-7)
                axs[1].set_ylim(bottom=-1, top=1)
                axs[0].set_ylim(bottom=-3, top=3)
                # On this signal, plot the theoretical weights
                T = len(d)
                t = np.arange(T)
                axs[1].plot(
                    0.2 + 0.3 * t / T, linestyle="--", color="black", linewidth=lw / 2
                )
                axs[1].plot(
                    0.7 * np.cos(4 * np.pi * t / T),
                    linestyle="--",
                    color="black",
                    linewidth=lw / 2,
                )
                axs[1].plot(
                    0.3 * np.sign(np.sin(10 * np.pi * t / T)) - 0.1,
                    linestyle="--",
                    color="black",
                    linewidth=lw / 2,
                )
                axs[1].plot(
                    -0.2 * np.ones(T), linestyle="--", color="black", linewidth=lw / 2
                )
            axs[0].set_xlim(left=0, right=len(d))
            axs[2].yaxis.set_minor_locator(AutoMinorLocator())
            axs[2].grid(b=True, which="major", color="silver", linestyle="-")
            axs[2].grid(b=True, which="minor", color="lightgray", linestyle="--")
            axs[2].set_ylabel("Error (dB)")
            fig.subplots_adjust(hspace=0.15)
            axs[-1].set_xlabel("Time/iteration")
            plt.savefig(
                f"../results/ex1_l{k+1}_n{n}_mu{int(mu*100):d}.png",
                dpi=600,
                bbox_inches="tight",
            )
            plt.savefig(
                f"../report/figures/pgf/ex1_l{k+1}_n{n}_mu{int(mu*100):d}.pgf",
                bbox_inches="tight",
            )
            plt.savefig(
                f"../report/figures/pdf/ex1_l{k+1}_n{n}_mu{int(mu*100):d}.pdf",
                bbox_inches="tight",
            )
            tikzplotlib.save(
                f"../report/figures/tikz/ex1_l{k+1}_n{n}_mu{int(mu*100):d}.tex"
            )
            plt.close()
