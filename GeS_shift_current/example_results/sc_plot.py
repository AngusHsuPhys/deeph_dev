import numpy as np
import matplotlib.pyplot as plt


plt.figure(figsize=(8, 4))
plt.rcParams.update({'font.size': 13,
                    'font.family': 'Arial',
                    'mathtext.fontset': 'cm'}
                    )
sc = np.loadtxt(f'shiftcond_111.dat')

plt.plot(sc1[:, 0], 2 * sc1[:, 1] * 15 / 100, label=r'$\sigma^{xxx}$', c='royalblue', linewidth=2, linestyle="--")
plt.xlabel(r'$\omega$ (eV)', fontsize=18)
plt.ylabel(r'Shift current $\sigma^{xxx}$ (10$^{2}$μA∙Å/V$^{2}$)', fontsize=18)
plt.legend(prop={'size': 18}, handlelength=1.5, frameon=False)
plt.tight_layout()

plt.savefig(f"sc.svg", format="svg", transparent=True)
