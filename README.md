# FSCOPE
## Numerical tool for Fluctuation Spectroscopy of superconductors
---
Supplemental code published in

**Fluctuation spectroscopy: From Rayleigh-Jeans waves to Abrikosov vortex clusters**

_A. A. Varlamov, A. Galda, and A. Glatz_

Rev. Mod. Phys. 90, 015009 (2018)

https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.90.015009

---

The source code and make file can be found in src/

- Compilation
-- run make_FSCOPE.sh to compile change compiler and flags as needed


- Usage:
+ run the executable without arguments for usage information

Parameters for the computation can be either defined on the command line directly or in a separate “ini” file.
Example:

> FSCOPE ctype=100 tmin=1.01 dt=0.01 Nt=100 hmin=0.1 dh=0.0 Nh=1 Tc0tau=0.01 Tc0tauphi=1 >> sigma.txt

or 

> FSCOPE sigma.ini >> sigma.txt

with file

sigma.ini
-------------------
ctype=100
tmin=1.01
dt=0.01
Nt=100
hmin=0.1
dh=0.0
Nh=1
Tc0tau=0.01
Tc0tauphi=1
-------------------

This example calculates 100 values for the fluctuation correction to conductivity at fixed field h=0.1 as function of temperature t from 1.01 to 2.00 in 0.01 steps.
