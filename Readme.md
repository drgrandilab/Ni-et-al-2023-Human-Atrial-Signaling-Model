### Source code for single cell, popoulations, and tissue modeling

H Ni, S Morotti, X Zhang, D Dobrev, E Grandi. Integrative human atrial modeling unravels interactive PKA and CaMKII signaling as key determinant of atrial arrhythmogenesis. bioRxiv, 2022

- Compiler: intel icpc
- Dependencies: CVODE and zlib
- Data analysis: Python & Matlab
- Data visulization: Paraview




#### Dependencies and complilers
* compiler: C++ with CXX=mpiicpc (tested with intel mpi for parallel computing on clusters)
* CVODE with Sundials (tested with Version 5.7)
* zlib (for tissue geometry input)





#### Folders

* Single_Cells: source code for single cell modeling
* Populations: source code for population simulation, DAD analysis, and logistic regression analysis
* 2D_Tissue_modeling: source code for 2D simulations
* 3D_Spatial_Myocyte_modeling: source code for running 3D spatial cell modeling
* PV-like_Populations: source code for simulating PV-like myocytes
