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

#%% Get the data

# Load the file
results = loadmat("../results/ex5_1.mat")
print(list(results.keys()))

W = results["W"]
error = results["xi"][:, 0]
niters = error.shape[0]
d = results["d"][0]

#%%
plt.figure(1)
plt.title("Constellation evolution for 8-PSK (ex. 5)")
lim = 2000
plt.scatter(
    np.real(d[:lim]), np.imag(d[:lim]), c=np.arange(1, lim + 1), cmap=cm.rainbow
)
plt.colorbar(label="Iteration")
ax = plt.gca()
ax.set_xlim(left=-5, right=5)
ax.set_ylim(bottom=-5, top=5)
ax.xaxis.set_minor_locator(AutoMinorLocator(n=2))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=2))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
fname = f"ex5_8psk_const"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")


#%%
xstep = 50
xnums = np.arange(1, niters + 1, xstep)
fig, axs = plt.subplots(num=5, nrows=2, ncols=1, sharex=True)
axs[0].plot(xnums, np.abs(W[:-1:xstep]))
axs[0].set_ylim(bottom=0, top=0.25)
axs[0].set_xlim(left=0, right=niters)
axs[1].plot(xnums, 10 * np.log10(error[::xstep]))
axs[1].set_xlabel("Iteration")
axs[0].set_ylabel("Weight value (absolute)")
axs[1].set_ylabel("Error (dB)")
axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[1].yaxis.set_minor_locator(AutoMinorLocator(n=2))
axs[0].yaxis.set_minor_locator(AutoMinorLocator(n=2))
axs[1].grid(b=True, which="major", color="silver", linestyle="-")
axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[0].grid(b=True, which="major", color="silver", linestyle="-")
axs[0].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[0].set_title("Results for 8-PSK (ex. 5)")
fname = f"ex5_8psk_res"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")


# %%
