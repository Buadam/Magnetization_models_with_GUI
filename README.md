# Mean Field magnetization models with FM or AFM interactions

## Description
This GUI was created to fit the magnetization m(H) and m(T) data measured on an anisotropic magnetic system with S=5/2 spin. The program offers 6 different mean field models for the curve fitting:


## Magnetization models:
Model functions located in the CODE\Model folder

Isotropic paramagnet: simple Brillouin's function with S=5/2. 

Isotropic ferromagnet: Self-consistent mean-field solution using Brillouin's function.

Isotropic antiferromagnet: Self-consistent mean-field solution for two sublattices. 

Anisotropic paramagnet with an axial anisotropy parameter, D. Magnetization curve obtained from the numerical solution of the spin-Hamiltonian. The axial symmetry of the system allows for a 2d solution without the loss of generality.

Anisotropic paramagnet with an axial an an in-plane anisotropy parameter, D and E, respectively. Magnetization curve obtained from the numerical solution of the spin-Hamiltonian. A full 3d solution of the model is required.

Anisotropic ferromagnet with a single axial anisotropy parameter, D: Magnetization solved self-consistently numerically from the spin-Hamiltonian. 

## Input data
Located in the DATA folder. Measured on malariapigment crystals in a suspension with and without magnetic-field cooling. See details in the paper below:

Á. Butykai, Á. Orbán, V. Kocsis, D. Szaller, S. Bordács, E. Tátrai-szekeres, L. F. Kiss, A. Bóta, B. G. Vértessy, T. Zelles, I. Kézsmárki, Malaria pigment crystals as magnetic micro-rotors: key for high-sensitivity diagnosis. Scientific Reports 3, 1431, (2013)

DOI: https://doi.org/10.1038/srep01431
