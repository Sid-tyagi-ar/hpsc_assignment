# DEM Particle Simulation with OpenMP

This project implements a **Discrete Element Method (DEM)** based particle simulator in C++, with **OpenMP parallelisation** for improved performance.  

The simulation models:
- Particle motion under gravity  
- Particle–particle interactions  
- Particle–wall interactions  

It also includes a Python script for generating plots to verify correctness and analyze behavior.

## Build and Run

### Compile
```bash
make
make run
```

This generates output data files (.txt) used for plotting.

# Generate Plots
make plot

This produces:
Free fall (numerical vs analytical)
Error vs timestep
Bounce height vs time
Kinetic energy vs time
Particle snapshots
# Clean Build
``` bash
make clean
```

Removes compiled binaries and generated files.


# Requirements
g++ with OpenMP support
Python 3
Python libraries:
numpy
matplotlib

# Install dependencies:
``` python
pip install numpy matplotlib
```

📬 Notes

This project was developed as part of a High Performance Scientific Computing (HPSC) assignment and demonstrates both numerical simulation and parallel computing concepts.
