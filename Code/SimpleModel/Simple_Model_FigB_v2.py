import matplotlib.pyplot as plt
import numpy as np
import numba as nb
import timeit
from matplotlib.ticker import MultipleLocator

# Same as Simple_Model_FigB.py but cleaner code

# Set a random seed for reproducibility
np.random.seed(42)

# Define a function to compute FD, using numba for performance optimization
@nb.njit()
def analytics(t, fs, fd, mus, mud):
    ns, nd = len(fs), len(fd)
    Ps, Pd = np.zeros(ns+1), np.zeros(nd)
    fs[::-1].sort(), fd[::-1].sort()

    # Calculate probabilities for single and double mutations
    for i in range(ns):
        Ps[i] = (1 - np.exp(-t * mus)) * np.exp(-t * mus * i)
    Ps[ns] = np.exp(-t * mus * ns)

    for i in range(nd):
        Pd[i] = (1 - np.exp(-t * mud)) * np.exp(-t * mud * i)

    # Calculate the FD matrix
    H = np.zeros((ns+1, nd))
    for i in range(ns):
        for j in range(nd):
            H[i, j] = fs[i] < fd[j]
    H[ns] = 1

    # Compute the final FD value
    FD = np.dot(Ps, np.dot(H, Pd)) / (1 - np.exp(-t * (ns * mus + nd * mud)))
    return FD

# Initialize parameters
ns, nd = 12, 20  # Number of single and double mutations
alpha = 4*10**(-3)  # Alpha parameter for mutation rate
mus = 1
mud = mus*alpha
mus = mus*3*0.76/5.89
mud = mud*3*0.99*0.52/9.63

landscapes = 10000  # Number of landscapes to simulate
datapoints = 60  # Number of data points
tmin, tmax = 10**(-1.2), 10**(4.2)  # Mutation supply range

# Create a logarithmically spaced array for mutation supply
T = np.logspace(np.log10(tmin), np.log10(tmax), num=datapoints, base=10)

# Initialize an array to store FD values across landscapes
all_y = np.zeros((landscapes, datapoints))

# Time the simulation
start = timeit.default_timer()
for l in range(landscapes):
    fs, fd = np.random.random(size=ns), np.random.random(size=nd)
    y = [analytics(T[i], fs, fd, mus, mud) for i in range(len(T))]
    all_y[l] = y
stop = timeit.default_timer()
print("Time:", stop - start)

# Plotting
fig = plt.figure(figsize=(6.3, 5))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
fig.subplots_adjust(bottom=0.2)

us = ns * (1 - np.exp(-mus * T))  # Average number of single mutations
ud = nd * (1 - np.exp(-mud * T))  # Average number of double mutations

# Plot lines
ax1.plot(T, all_y.mean(axis=0), color="black", lw=2, label="$\\overline{F}_D$")
ax1.plot(T, [nd * (mud) / (nd * mud + ns * mus)] * len(T), "--", lw=1.82, color="black", label="$f_b$")
ax1.plot(T, us, lw=2, label='$\\overline{u}_S$', color="#62AAC5")
ax1.axvline(x=1/mus, color="black", lw=2, linestyle="dotted")
ax1.plot(T, ud, lw=2, label='$\\overline{u}_D$', color="#DD903B")
ax1.plot(T, ud / (ud + us), lw=2, label='$\\overline{u}_D/(\\overline{u}_D + \\overline{u}_S$)', color="#a45e5f")

# Annotations
ax1.annotate(xy=(14000, 20), xytext=(5, 0), textcoords='offset points', text="20", va='center', color="#DD903B")
ax1.annotate(xy=(14000, 12), xytext=(5, 0), textcoords='offset points', text="12", va='center', color="#62AAC5")
ax1.annotate(xy=(14000, 20/32), xytext=(5, 0), textcoords='offset points', text="$\\frac{20}{12+20}$", va='center', color="black")

# Set scales and ticks
ax1.set_xscale("log"), ax1.set_yscale("log")
ax1.set_xticks([1/(mus*ns), 12/(mus*ns), 100/(mus*ns), 10000/(mus*ns), 100000/(mus*ns), 10000000/(mus*ns)])
ax1.set_xticklabels([1, 12, 100, 10000, 100000, 10000000])
ax1.tick_params(axis='x', labelsize=10)
ax1.xaxis.set_minor_locator(MultipleLocator(100000))

# Configure secondary x-axis (ax2)
ax2.set_xscale("log")
ax2.xaxis.set_ticks_position("bottom"), ax2.xaxis.set_label_position("bottom")
ax2.spines["bottom"].set_position(("axes", -0.1))
ax2.set_frame_on(True), ax2.patch.set_visible(False)
for sp in ax2.spines.values(): sp.set_visible(False)
ax2.spines["bottom"].set_visible(True)
ax2.set_xticks([0.01/(mud*nd),1/(mud*nd),20/(mud*nd),100/(mud*nd),10000/(mud*nd),100000/(mud*nd)])
ax2.set_xticklabels([0.01, 1, 20, 100, 10000,100000])
ax2.tick_params(axis='x', labelsize=10)
ax2.xaxis.set_minor_locator(MultipleLocator(100000))

# Set labels and limits
ax1.xaxis.set_label_coords(-.008, -.03), ax2.xaxis.set_label_coords(.000, -.13)
ax1.set_ylabel("$F_{D}$", fontsize=14)
ax1.legend(fontsize=13, framealpha=0.9, loc=(0.667, 0.012), ncol=1, handlelength=1.3)
plt.ylim(0.0005, 28), ax1.set_xlim(tmin, tmax), ax2.set_xlim(tmin, tmax)

# Save and display the plot
plt.savefig('fig_B.pdf')
plt.show()
