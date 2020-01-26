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

# plt.ioff()

#%% For x1

# Load the file
results1 = loadmat("../results/results3_1.mat")

print("results1:")
print(list(results1.keys()))

W1 = results1["W"]
error1 = 10 * np.log10(results1["errav"])
N1 = results1["ns"]
nw1 = results1["n"]

plt.figure(1)
plt.plot(np.arange(0, N1, 200), error1[::200, 0])
ax = plt.gca()
ax.set_xlim(left=0, right=N1)
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_ylabel("Error (dB)")
ax.set_xlabel("Iteration")
ax.set_title("Error evolution for \(x_1\)")
ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
fname = f"ex3_x1_error"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
plt.close()

plt.figure(2)
for i in range(W1.shape[0]):
    plt.plot(np.arange(0, N1 + 1, 200), W1[i, ::200], "-", linewidth=1)
ax = plt.gca()
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=3))
ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_xlim(left=0, right=N1)
ax.set_ylabel("Weight value")
ax.set_xlabel("Iteration")
ax.set_title("Evolution of the weights for \(x_1\)")
fname = f"ex3_x1_weights"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
plt.close()

fig, axs = plt.subplots(num=5, nrows=2, ncols=1, sharex=True)
axs[0].set_title("Results for ex. 3, \(x_1\)")
axs[1].plot(np.arange(0, N1, 200), error1[::200, 0])
axs[1].set_xlim(left=0, right=N1)
axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=5))
axs[1].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[1].grid(b=True, which="major", color="silver", linestyle="-")
axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[1].set_ylabel("Error (dB)")
axs[1].set_xlabel("Iteration")
for i in range(W1.shape[0]):
    axs[0].plot(np.arange(0, N1 + 1, 200), W1[i, ::200], "-", linewidth=1)
axs[0].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[0].grid(b=True, which="major", color="silver", linestyle="-")
axs[0].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
axs[0].set_ylabel("Weight value")
fig.subplots_adjust(hspace=0.15)
fname = f"ex3_x1_all"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
plt.close()

#%% For the part 2
del results1, W1, error1, N1, nw1

results2 = loadmat("../results/results3_2.mat")

print("results2:")
print(list(results2.keys()))

W2 = results2["W"]
error2 = 10 * np.log10(results2["errav"])
N2 = results2["ns"]
nw2 = results2["n"]

plt.figure(3)
plt.plot(np.arange(0, N2, 200), error2[::200, 0])
ax = plt.gca()
ax.set_xlim(left=0, right=N2)
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_ylabel("Error (dB)")
ax.set_xlabel("Iteration")
ax.set_title("Error evolution for \(x_2\)")
ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
fname = f"ex3_x2_error"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
plt.close()

plt.figure(4)
for i in range(W2.shape[0]):
    plt.plot(np.arange(0, N2 + 1, 200), W2[i, ::200], "-", linewidth=1)
ax = plt.gca()
ax.xaxis.set_minor_locator(AutoMinorLocator(n=5))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=3))
ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
ax.grid(b=True, which="major", color="silver", linestyle="-")
ax.grid(b=True, which="minor", color="lightgray", linestyle="--")
ax.set_xlim(left=0, right=N2)
ax.set_ylabel("Weight value")
ax.set_xlabel("Iteration")
ax.set_title("Evolution of the weights for \(x_2\)")
fname = f"ex3_x2_weights"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
plt.close()

fig, axs = plt.subplots(num=6, nrows=2, ncols=1, sharex=True)
axs[0].set_title("Results for ex. 3, \(x_2\)")
axs[1].plot(np.arange(0, N2, 200), error2[::200, 0])
axs[1].set_xlim(left=0, right=N2)
axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=5))
axs[1].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[1].grid(b=True, which="major", color="silver", linestyle="-")
axs[1].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[1].set_ylabel("Error (dB)")
axs[1].set_xlabel("Iteration")
for i in range(W2.shape[0]):
    axs[0].plot(np.arange(0, N2 + 1, 200), W2[i, ::200], "-", linewidth=1)
axs[0].yaxis.set_minor_locator(AutoMinorLocator(n=3))
axs[0].grid(b=True, which="major", color="silver", linestyle="-")
axs[0].grid(b=True, which="minor", color="lightgray", linestyle="--")
axs[1].ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
axs[0].set_ylabel("Weight value")
fig.subplots_adjust(hspace=0.15)
fname = f"ex3_x2_all"
plt.savefig("../results/" + fname + ".png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pdf/" + fname + ".pdf", bbox_inches="tight")
plt.close()
