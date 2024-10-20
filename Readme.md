### Source code for single cell, popoulations, and tissue modeling

H Ni, S Morotti, X Zhang, D Dobrev, E Grandi. Integrative human atrial modelling unravels interactive protein kinase A and Ca2+/calmodulin-dependent protein kinase II signalling as key determinants of atrial arrhythmogenesis. 
Cardiovasc Res. 2023 Oct 24;119(13):2294-2311. doi: 10.1093/cvr/cvad118. 

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
* Matlab_Single_Cell: a copy of single cell model implemented in Matlab 
