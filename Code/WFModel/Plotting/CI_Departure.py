import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from scipy.special import lambertw

'''
Initial departure of FD from fb for different DBFE in the Wright-Fisher model.
This code compares the results of our simulations with the analytical results described in the SI.
'''

def calc_harmonic(n):
    return sum(1/d for d in range(1, n + 1))


ns = 120
nd = 200
s = 0.1

########################################################
# Exponential DBFE
folder = 'dataRandom/meanS01D008ns120nd200_exp'


ris0 = []
ris = []
f = []

for sfile in sorted(glob.glob(folder + '/*ris.npy')):
    result = re.search('' + folder + '/(.*)_', sfile)
    seed = result.group(1)
    ris.append(np.load(folder + "/" + seed + "_ris.npy"))
    ris0.append(np.load(folder + "/" + seed + "_ris0.npy"))
    f.append(np.load(folder + "/" + seed + "_fit.npy"))

x = np.load(glob.glob(folder + '/*MS.npy')[0])
ris = np.array(ris)
ris0 = np.array(ris0)
f = np.array(f)

print("ris.shape", ris.shape)

landscapes = ris.shape[0]

alpha = 4 * 10 ** (-3)
mus = 10 ** (-10)
mud = 10 ** (-10) * alpha
mus = mus * 3 * 0.76
mud = mud * 3 * 0.99 * 0.52

U = ns * mus / 5.80 + nd * mud / 9.63
Us = ns * mus / 5.80
fractionAAlevel = mud / 9.63 * nd / (mud / 9.63 * nd + mus / 5.80 * ns)

print("mus", mus)
print("fractionAAlevel", fractionAAlevel)

Avg_meanline = 0

for i in range(landscapes):
    w_SN = np.random.exponential(scale=0.1, size=ns)
    w_DN = np.random.exponential(scale=0.1, size=nd)
    meanSNB = np.mean(w_SN)
    meanDNB = np.mean(w_DN)
    meanline = fractionAAlevel * 2 * meanDNB / (fractionAAlevel * 2 * meanDNB + (1 - fractionAAlevel) * 2 * meanSNB)

    Avg_meanline += meanline

fraction_of_time_a_AAsub_only_due_to_doubles = Avg_meanline / landscapes
print("Mean fraction", fraction_of_time_a_AAsub_only_due_to_doubles)
print("ris[:,:,ns+1:].sum(axis=2).mean(axis=0)[0]", ris[:, :, ns + 1:].sum(axis=2).mean(axis=0)[0])

N_list = []
for i in range(len(x)):
    ms = x[i]
    N = np.floor(ms / U)
    N_list.append(N)
N_list = np.array(N_list)



ii = 0
for r in ris:
    if ii > 500 and ii <= 600:
        plt.plot(x, np.sum(r[:, ns + 1:], axis=1), color="grey", alpha=0.3)
    ii += 1

plt.plot(x, ris[:, :, ns + 1:].sum(axis=2).mean(axis=0), color="black", linewidth=3, label="$\\overline{F}_{D}$")
plt.axhline([fraction_of_time_a_AAsub_only_due_to_doubles], linestyle="--", color="black", linewidth=3,
                label="$f_b$")

epsilon = 1
plt.axvline([np.real(epsilon / (4 * lambertw(epsilon / (4 * U))))], linestyle=(0, (3, 7)), color="black",
                linewidth=2.5)
plt.axvline([np.real(epsilon * ns / (4 * lambertw(epsilon / (4 * mus / 5.80))))], linestyle=(0, (1, 5)),
                color="black", linewidth=2.5)

fbtest = fraction_of_time_a_AAsub_only_due_to_doubles
mD = 0.08
mS = 0.1
mustest = mus / 5.80
mudtest = mud / 9.63

FD4 = fbtest * mD / (fbtest * mD + (1 - fbtest) * mS)


ratiptest = mS * (mD - mS) * (4 * mD * mS + mD ** 2 + mS ** 2) / (mD + mS) ** 3
print("ratiptest", ratiptest)


print(mustest * ns * N_list)
print(x)


UD = nd*mudtest
US = ns*mustest


FD6 = 1/(1 + (US*(mS + 1/2*N_list*(-((4*mD**2*(mD + 2*mS)*UD)/(mD + mS)**2) - 3*mS*US)*np.log(N_list)))/(UD*(mD + 1/2*N_list*(-3*mD*UD - (4*mS**2*(2*mD + mS)*US)/(mD + mS)**2)*np.log(N_list))))
plt.plot(x, FD6, label="FD6")

FD7 = FD4 + ((mD - mS)*UD*US*(4*mD**3*UD + 4*mS**3*US + 3*mD**2*mS*(3*UD + US) + 3*mD*mS**2*(UD + 3*US))*N_list*np.log(N_list))/(2*(mD + mS)**2*(mD*UD + mS*US)**2)
plt.plot(x, FD7, label="FD7")


print(fbtest * mD / (fbtest * mD + (1 - fbtest) * mS))
print(UD * mD / (UD*mD+US*mS))


plt.ylim(0.001, 0.01)
plt.yscale("log")

plt.xscale("log")
plt.legend()

plt.show()