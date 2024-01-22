# Import required libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


# Parameter initialization for the three types of mutations (S, D, and I)
ns = 10
nd = 100
ni = 1000
mus = 1
mud = 10**(-2)
mui = 10**(-4)


datapoints = 30
tmin = 10**(-3)
tmax = 10**(6)

# Create logarithmically spaced array for mutation supply
T = np.logspace(np.log10(tmin), np.log10(tmax), num=datapoints, base=10)  # mutation supply


us = ns*(1-np.exp(-mus*T))
ud = nd*(1-np.exp(-mud*T))
ui = ni*(1-np.exp(-mui*T))


# Plotting
fig = plt.figure(figsize=(5.7,5.7))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax3 = ax1.twiny()


# Adjusting layout for the twin X axes
fig.subplots_adjust(bottom=0.2)

ax1.plot(T, us, color="#62AAC5", lw=2, label='$\\overline{u}_S$')
ax1.plot(T, ud, color="#DD903B", lw=2, label='$\\overline{u}_D$')
ax1.plot(T, ui, color="#446455", lw=2, label='$\\overline{u}_I$')
ax1.plot(T, ud/(ud+us+ui), lw=2, color="#a45e5f", label='$\\overline{u}_D$/($\\overline{u}_D$+$\\overline{u}_S$+$\\overline{u}_I$)')
ax1.plot(T, ui/(ud+us+ui), lw=2, color="C4", label='$\\overline{u}_I$/($\\overline{u}_D$+$\\overline{u}_S$+$\\overline{u}_I$)')


ax1.set_xscale("log")
ax1.set_yscale("log")

# Setting ticks for type S
# Set the x-ticks
ax1.set_xticks([0.1/(mus*ns), 1/(mus*ns), 10/(mus*ns), 100/(mus*ns), 1000/(mus*ns), 10000/(mus*ns), 100000/(mus*ns), 10**6/(mus*ns),10**7/(mus*ns)])
ax1.set_xticklabels([0.1, 1, 10, "$10^2$", "$10^3$", "$10^4$","$10^5$","$10^6$","$10^7$"])

# Set the fontsize for x-tick labels
ax1.tick_params(axis='x', labelsize=10)


ax1.set_xlabel("# S", fontsize=12)


# Settings for the twin X axis for type D
ax2.set_xscale("log")
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")
ax2.spines["bottom"].set_position(("axes", -0.1))
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
for sp in ax2.spines.values():
    sp.set_visible(False)
ax2.spines["bottom"].set_visible(True)

# Set the x-ticks
ax1.set_xticks([0.01/(mus*nd), 0.1/(mus*nd), 1/(mus*nd), 10/(mus*nd), 100/(mus*nd), 1000/(mus*nd), 10000/(mus*nd), 100000/(mus*ns), 10**6/(mus*ns),10**7/(mus*ns)])
ax1.set_xticklabels([0.01,0.1, 1, 10, "$10^2$", "$10^3$", "$10^4$","$10^5$","$10^6$","$10^7$"])

# Set the fontsize for x-tick labels
ax1.tick_params(axis='x', labelsize=10)


ax2.set_xlabel("# D", fontsize=12)

# Settings for the second twin X axis for type I
ax3.set_xscale("log")
ax3.xaxis.set_ticks_position("bottom")
ax3.xaxis.set_label_position("bottom")
ax3.spines["bottom"].set_position(("axes", -0.2))
ax3.set_frame_on(True)
ax3.patch.set_visible(False)
for sp in ax3.spines.values():
    sp.set_visible(False)
ax3.spines["bottom"].set_visible(True)

# Set the x-ticks
ax1.set_xticks([0.001/(mui*ni),0.01/(mui*ni),0.1/(mui*ni),1/(mui*ni),10/(mui*ni),100/(mui*ni),1000/(mui*ni),10000/(mui*ni),100000/(mui*ni), 10**6/(mui*ni)])
ax1.set_xticklabels([0.001,0.01, 0.1, 1, 10, "$10^2$", "$10^3$", "$10^4$","$10^5$","$10^6$"])

# Set the fontsize for x-tick labels
ax1.tick_params(axis='x', labelsize=10)

ax3.set_xlabel("# I", fontsize=12)


ax1.xaxis.set_minor_locator(MultipleLocator(10000000))
ax2.xaxis.set_minor_locator(MultipleLocator(10000000))
ax3.xaxis.set_minor_locator(MultipleLocator(10000000))
ax1.xaxis.set_label_coords(.00, -.025)
ax2.xaxis.set_label_coords(.000, -.125)
ax3.xaxis.set_label_coords(.000, -.225)


ax1.set_xlim(tmin,tmax)
ax2.set_xlim(tmin,tmax)
ax3.set_xlim(tmin,tmax)
ax1.set_ylim(3*10**(-5),2*10**(3))

plt.title("$n_S={}$, $n_D={}$, $n_I={}$".format(ns,nd,ni), fontsize=14)
ax1.legend(fontsize=13, framealpha=1, loc=4)

#plt.savefig('simpleModelthreeTypes.pdf')
plt.show()