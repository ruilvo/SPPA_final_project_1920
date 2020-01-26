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

# plt.ioff()

#%% Get AWGN results from files
# Load the file
results1 = loadmat("../results/results22_1.mat")
results2 = loadmat("../results/results22_2.mat")

#%% Get variables
print("results1:")
print(list(results1.keys()))
print("results2:")
print(list(results2.keys()))

#%% Load them
W1 = results1["W"]
W2 = results2["W"]
error1 = results1["EDb"]
error2 = results2["EDb"]
N1 = results1["N"]
N2 = results2["N"]

ns1 = np.arange(1, N1 + 2)
ns2 = np.arange(1, N2 + 2)


#%% Show beautiful figures
plt.figure(1)
plt.plot(ns1, error1[0])
ax = plt.gca()
ax.set_xlim(left=0, right=ns1[-1])
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_ylabel("Error (dB)")
ax.set_xlabel("Number of weights, \(n\)")
ax.set_title("Error in relation to the number of weights for \(x_1\)")

fname = f"ex2_pt2_x1_error"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/" + fname + ".pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
tikzplotlib.save("../report/figures/tikz/" + fname + ".tex")

plt.figure(2)
for n in np.arange(1, N1 + 2):
    plt.plot(np.arange(n, N1 + 2), W1[n - 1 :, n - 1], "-")
ax = plt.gca()
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=3))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_xlim(left=0, right=ns1[-1])
ax.set_ylabel("Weight value")
ax.set_xlabel("Number of weights, \(n\)")
ax.set_title("Evolution of the weights for \(x_1\)")

fname = f"ex2_pt2_x1_weights"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/" + fname + ".pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
tikzplotlib.save("../report/figures/tikz/" + fname + ".tex")

plt.figure(3)
plt.plot(ns2, error2[0])
ax = plt.gca()
ax.set_xlim(left=0, right=ns2[-1])
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_ylabel("Error (dB)")
ax.set_xlabel("Number of weights, \(n\)")
ax.set_title("Error in relation to the number of weights for \(x_2\)")

fname = f"ex2_pt2_x2_error"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/" + fname + ".pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
tikzplotlib.save("../report/figures/tikz/" + fname + ".tex")

plt.figure(4)
for n in np.arange(1, N2 + 2):
    plt.plot(np.arange(n, N2 + 2), W2[n - 1 :, n - 1], "-")
ax = plt.gca()
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=3))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_xlim(left=0, right=ns2[-1])
ax.set_ylabel("Weight value")
ax.set_xlabel("Number of weights, \(n\)")
ax.set_title("Evolution of the weights for \(x_2\)")

fname = f"ex2_pt2_x2_weights"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/" + fname + ".pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
tikzplotlib.save("../report/figures/tikz/" + fname + ".tex")


#%%
fig, axs = plt.subplots(num=5, nrows=2, ncols=1, sharex=True)
axs[0].set_title("Results for ex. 2, \(x_1\)")
axs[1].plot(ns1, error1[0])
axs[1].set_xlim(left=0, right=ns1[-1])
axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=5))
axs[1].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[1].grid(b=True, which="major", color="silver", linestyle="-")
axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[1].set_ylabel("Error (dB)")
axs[1].set_xlabel("Number of weights, \(n\)")
for n in np.arange(1, N1 + 2):
    axs[0].plot(np.arange(n, N1 + 2), W1[n - 1 :, n - 1], "-")
axs[0].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[0].grid(b=True, which="major", color="silver", linestyle="-")
axs[0].grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_ylabel("Weight value")
fig.subplots_adjust(hspace=0.15)

fname = f"ex2_pt2_x1_all"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/" + fname + ".pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
tikzplotlib.save("../report/figures/tikz/" + fname + ".tex")

fig, axs = plt.subplots(num=6, nrows=2, ncols=1, sharex=True)
axs[0].set_title("Results for ex. 2, \(x_2\)")
axs[1].plot(ns2, error2[0])
axs[1].set_xlim(left=0, right=ns2[-1])
axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=5))
axs[1].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[1].grid(b=True, which="major", color="silver", linestyle="-")
axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[1].set_ylabel("Error (dB)")
axs[1].set_xlabel("Number of weights, \(n\)")
for n in np.arange(1, N2 + 2):
    axs[0].plot(np.arange(n, N2 + 2), W2[n - 1 :, n - 1], "-")
axs[0].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[0].grid(b=True, which="major", color="silver", linestyle="-")
axs[0].grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_ylabel("Weight value")
fig.subplots_adjust(hspace=0.15)

fname = f"ex2_pt2_x2_all"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/" + fname + ".pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
tikzplotlib.save("../report/figures/tikz/" + fname + ".tex")

#%%
plt.show()
