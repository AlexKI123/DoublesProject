import matplotlib.pyplot as plt

# Setting the maximum number of codons to be considered
CodonsMax = 4

# Pre-calculated results from AA_accessible for different number of considered codons
results = [[5.80327868852459, 11.475409836065573, 8.475409836065573, 0.24408014571949, 0.0036429872495446266], [5.80327868852459, 13.228701961838215, 9.05227089492072, 0.24408014571949, 0.008934276926752067], [5.803278688524342, 13.810360044820778, 9.244557914539358, 0.24408014571949, 0.010257099346053929], [5.80327868852459, 14.101189086311189, 9.340701424348293, 0.24408014571949, 0.01085838226391841]]


plt.figure(figsize=(5.5, 4))
ax = plt.subplot(111)

ax.plot(range(1,CodonsMax+1), [item[0] for item in results],".-", color='#3F9DB7', label="S",  markersize=10)
ax.plot(range(1,CodonsMax+1), [item[1] for item in results],".-", color='Red', label="D", markersize=10)
ax.plot(range(1,CodonsMax+1), [item[2] for item in results],".-", color='#EA8B1C', label="D \ S", markersize=10)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)
plt.xlabel("# Codons",fontsize=12)
plt.ylabel("<# Acc. AA seq. per codon>",fontsize=12)

plt.xticks(range(1,CodonsMax+1))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)


# results AA_accessible_on_ring
plt.axhline(y=5.803278688524342, color='#3F9DB7', linestyle='--')
plt.axhline(y=14.973676210784909, color='Red', linestyle='--')
plt.axhline(y=9.629131953776472, color='#EA8B1C', linestyle='--')

plt.ylim(0, 16)

plt.savefig('Fig.pdf', bbox_inches='tight')

plt.show()
