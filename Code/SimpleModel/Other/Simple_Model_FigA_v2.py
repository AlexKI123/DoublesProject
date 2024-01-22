import matplotlib.pyplot as plt
import numpy as np
import numba as nb
import timeit
from matplotlib.ticker import MultipleLocator

# Same as Simple_Model_FigA.py but cleaner code

# Function to compute FD
@nb.njit()
def analytics(t, fs, fd, mus, mud):
    ns, nd = len(fs), len(fd)
    Ps, Pd = np.zeros(ns+1), np.zeros(nd)
    fs[::-1].sort(), fd[::-1].sort()

    for i in range(ns):
        Ps[i] = (1 - np.exp(-t * mus)) * np.exp(-t * mus * i)
    Ps[ns] = np.exp(-t * mus * ns)

    for i in range(nd):
        Pd[i] = (1 - np.exp(-t * mud)) * np.exp(-t * mud * i)

    H = np.zeros((ns+1, nd))
    for i in range(ns):
        for j in range(nd):
            H[i, j] = fs[i] < fd[j]
    H[ns] = 1

    FD = np.dot(Ps, np.dot(H, Pd)) / (1 - np.exp(-t * (ns * mus + nd * mud)))
    return FD


def generate_data(ns, nd, mus, mud, landscapes, datapoints, tmin, tmax):
    T = np.logspace(np.log10(tmin), np.log10(tmax), num=datapoints, base=10)
    all_y = np.zeros((landscapes, datapoints))

    for l in range(landscapes):
        fs, fd = np.random.random(size=ns), np.random.random(size=nd)
        y = np.array([analytics(T[i], fs, fd, mus, mud) for i in range(len(T))])
        all_y[l] = y

    return T, all_y


def plot_data(T, all_y, mus, ns, mud, nd):
    fig = plt.figure(figsize=(4.2, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    fig.subplots_adjust(bottom=0.2)

    for i in range(60):
        ax1.plot(T, all_y[i], color="gray", lw=0.5)
    ax1.plot(T, all_y.mean(axis=0), color="black", lw=1.6, label='$\\overline{F}_D$')
    ax1.plot(T, [nd * (mud) / (nd * mud + ns * mus)] * len(T), "--", color="black", lw=1.6, label='$f_b$')

    ax1.set_xscale("log"), ax2.set_xscale("log")
    ax1.set_xticks([1/(mus * ns) * 10**i for i in range(5)])
    ax1.set_xticklabels([10**i for i in range(5)])
    ax1.tick_params(axis='x', labelsize=10)
    ax1.set_xlabel("# Singles", fontsize=12)
    ax1.set_ylabel("$F_{D}$", fontsize=16)

    ax2.xaxis.set_ticks_position("bottom"), ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.1))
    ax2.set_frame_on(True), ax2.patch.set_visible(False)
    for sp in ax2.spines.values(): sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)

    ax2.set_xticks([0.01/(mud * nd) * 10**i for i in range(7)])
    ax2.set_xticklabels([10**i for i in range(-2, 5)])
    ax2.tick_params(axis='x', labelsize=10)
    ax2.set_xlabel("# Doubles", fontsize=12)

    ax1.xaxis.set_minor_locator(MultipleLocator(10000))
    ax2.xaxis.set_minor_locator(MultipleLocator(10000))

    ax1.xaxis.set_label_coords(.069, -.024)
    ax2.xaxis.set_label_coords(.08, -.12)

    ax1.legend(fontsize=12, loc=(0.02, 0.6))
    plt.ylim(0.0, 0.010), ax1.set_xlim(tmin, tmax), ax2.set_xlim(tmin, tmax)
    plt.tight_layout()
    plt.savefig('fig_A.pdf')
    plt.show()


# Parameters
ns, nd = 12, 20
alpha = 4*10**(-3)
mus = 1
mud = mus*alpha
mus = mus*3*0.76/5.89
mud = mud*3*0.99*0.52/9.63
landscapes, datapoints = 100, 60
tmin, tmax = 10**(-1.2), 10**(1.2)

np.random.seed(42)

# Generate data
start = timeit.default_timer()
T, all_y = generate_data(ns, nd, mus, mud, landscapes, datapoints, tmin, tmax)
stop = timeit.default_timer()
print("Time: ", stop - start)

# Plotting
plot_data(T, all_y, mus, ns, mud, nd)
