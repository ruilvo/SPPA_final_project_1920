"""
A Python script to make beautiful plots
"""

#%% Imports
import glob
from contextlib import suppress
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.cm as cm
import matplotlib


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
        "pgf.rcfonts": True,
    }
)

# plt.rcParams["figure.autolayout"] = True

# Adjust figure default size
plt.rcParams["figure.figsize"] = (6, 4)  # Default fig size (w, h))

# Get the default color scheme
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

#%%
files = glob.glob("../results/results_ex4*.mat")
print(files)


# %%
for i in range(len(files)):
    fname = files[i]
    results = loadmat(fname)
    print(list(results.keys()))
    W = results["W"]
    error = results["err"][:, 0]
    mu = results["mu"][0, 0]
    n = results["n"][0, 0]
    d = results["d"][0]
    ys = results["ys"][0]

    dname = "\(d_1\)" if "d1" in fname else "\(d_2\)"

    fig, axs = plt.subplots(num=i + 1, nrows=3, ncols=1, sharex=True)
    title = f"Results ex. 5 for {dname}, \(n={n:d}\) \(\mu={mu:.1f}\)"
    print(title)
    axs[0].plot(d, label="Target signal")
    axs[0].plot(ys, label="Recovered signal")
    axs[0].plot(ys - d, label="Difference")
    axs[0].legend(loc="lower right", fontsize=6)
    axs[1].plot(W)
    axs[2].plot(10 * np.log10(error))
    axs[2].set_xlim(left=0, right=len(d))
    axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=2))
    axs[1].yaxis.set_minor_locator(AutoMinorLocator(n=2))
    axs[1].grid(b=True, which="major", color="silver", linestyle="-")
    axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
    axs[2].xaxis.set_minor_locator(AutoMinorLocator(n=2))
    axs[2].yaxis.set_minor_locator(AutoMinorLocator(n=2))
    axs[2].grid(b=True, which="major", color="silver", linestyle="-")
    axs[2].grid(b=True, which="minor", color="lightgray", linestyle="--")
    axs[2].set_xlabel("Iteration")
    axs[2].set_ylabel("Error (dB)")
    axs[1].set_ylabel("Weights value")
    axs[0].set_ylabel("Amplitude")
    axs[0].set_title(title)

    plt.savefig(fname.replace(".mat", ".png"), dpi=600, bbox_inches="tight")
    pdfsavename = fname.replace(".mat", ".pdf").replace(
        "../results", "../report/figures/pdf"
    )
    plt.savefig(pdfsavename, bbox_inches="tight")

    plt.close()


# %%
