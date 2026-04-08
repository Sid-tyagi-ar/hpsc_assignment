import numpy as np
import matplotlib.pyplot as plt

files = ["snap_0.txt", "snap_500.txt", "snap_1000.txt", "snap_1500.txt"]

fig, axs = plt.subplots(2, 2, figsize=(8,8))

axs = axs.flatten()

for i, f in enumerate(files):
    data = np.loadtxt(f)
    axs[i].scatter(data[:,0], data[:,1], s=2)
    
    t_val = f.split('_')[1].split('.')[0]
    axs[i].set_title(f"t = {t_val}")
    
    axs[i].set_xticks([])
    axs[i].set_yticks([])

plt.tight_layout()
plt.savefig("snapshots.png", dpi=300)
plt.show()
