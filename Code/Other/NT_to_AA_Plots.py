import matplotlib.pyplot as plt
import numpy as np


# Plot the result of the transformation simulations

# Constants
CODONS = 30
L = CODONS * 3
LOG_MIN = np.log10(10**(-5))
LOG_MAX = np.log10(10**(-3))
LOG_POINTS = 10

# Plotting function
def plot_graph(x_theory, y_theory, x_simu, y_simu, x_label, y_label, file_name):
    plt.figure(figsize=(3.8, 3.8))
    plt.title("$codons=30$, $generations=10^7$")
    plt.plot(x_theory, y_theory, "-", label="Expectation", color="black")
    plt.plot(x_simu, y_simu, "D", label="Simulation", color="black")
    plt.xlabel(x_label, fontsize=13)
    plt.ylabel(y_label, fontsize=13)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize=13)
    plt.tight_layout()
    plt.savefig(file_name, bbox_inches='tight')
    plt.show()

# (tilde_nu_S -> tilde_mu_S)
tilde_nu_S_simu = [10**(-5), 10**(-4), 10**(-3)]
tilde_mu_S_simu = [2.2163333333333332e-05, 0.00022842666666666667, 0.00228063]
tilde_nu_S_theory = np.logspace(LOG_MIN, LOG_MAX, LOG_POINTS)
tilde_mu_S_theory = tilde_nu_S_theory * 3 * (1 - 0.23958333333333334)
plot_graph(tilde_nu_S_theory, tilde_mu_S_theory, tilde_nu_S_simu, tilde_mu_S_simu, '$\\tilde\\nu_S$', '$\\tilde\\mu_S$', 'Fig1.pdf')

# (tilde_nu_D -> tilde_mu_D)
tilde_nu_D_simu = [10**(-5), 10**(-4), 10**(-3)]
tilde_mu_D_simu = [1.5446666666666666e-05, 0.00015396666666666668, 0.0015293466666666668]
tilde_nu_D_theory = np.logspace(LOG_MIN, LOG_MAX, LOG_POINTS)
tilde_mu_D_theory = (L - 1) / CODONS * tilde_nu_D_theory * (1 - 0.011363636363636364) * 0.5185
plot_graph(tilde_nu_D_theory, tilde_mu_D_theory, tilde_nu_D_simu, tilde_mu_D_simu, '$\\tilde\\nu_D$', '$\\tilde\\mu_D$', 'Fig2.pdf')

# tilde_nu_S -> mu_S
mu_S_simu = [3.729786871270247e-06, 3.756958226768969e-05, 0.0003737703665814152]
mu_S_theory = tilde_nu_S_theory * 3 * (1 - 0.23958333333333334) / 6.11
plot_graph(tilde_nu_S_theory, mu_S_theory, tilde_nu_S_simu, mu_S_simu, '$\\tilde\\nu_S$', '$\\mu_S$', 'Fig3.pdf')

# tilde_nu_D -> mu_D
mu_D_simu = [1.6078878378090972e-06, 1.6051366202544487e-05, 0.00016086807171646484]
mu_D_theory = (L - 1) / CODONS * tilde_nu_D_theory * (1 - 0.011363) * 0.5185 / 10.375
plot_graph(tilde_nu_D_theory, mu_D_theory, tilde_nu_D_simu, mu_D_simu, '$\\tilde\\nu_D$', '$\\mu_D$', 'Fig4.pdf')
