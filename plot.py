import numpy as np
import matplotlib.pyplot as plt
import glob

# -------- Free fall --------
data = np.loadtxt("freefall.txt")
t = data[:,0]
z = data[:,1]

g = 9.81
z_analytical = z[0] - 0.5*g*t**2

plt.figure()
plt.plot(t,z,label="Numerical")
plt.plot(t,z_analytical,'--',label="Analytical")
plt.title("Free Fall")
plt.xlabel("Time")
plt.ylabel("Height")
plt.legend()

# -------- Error vs dt --------
data = np.loadtxt("error_dt.txt")
plt.figure()
plt.loglog(data[:,0], data[:,1], marker='o')
plt.xlabel("dt")
plt.ylabel("Error")
plt.title("Error vs Timestep")

# -------- Bounce --------
data = np.loadtxt("bounce.txt")
plt.figure()
plt.plot(data[:,0], data[:,1])
plt.xlabel("Time")
plt.ylabel("Height")
plt.title("Bounce Height")

# -------- Energy --------
data = np.loadtxt("energy.txt")
plt.figure()
plt.plot(data[:,0], data[:,1])
plt.xlabel("Time")
plt.ylabel("Kinetic Energy")
plt.title("Energy vs Time")

# -------- Snapshots --------
files = sorted(glob.glob("snap_*.txt"))
for f in files:
    data = np.loadtxt(f)
    plt.figure()
    plt.scatter(data[:,0], data[:,1], s=2)
    plt.title(f)

plt.show()
