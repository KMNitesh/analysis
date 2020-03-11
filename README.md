## Molecular dynamics trajectory analysis package for Tinker & Gromacs & Amber software

#### Architecture
<img src="doc/Architecture.png" alt="Software Architecture" title="Software Architecture" width="500" height="359" />


#### Typical analysis functions
 
1.   trajectory format conversion between Tinker ARC, Gromacs Trr, Gromacs Xtc, Gromacs Gro and Amber Netcdf
2.   Diffusion coefficient based on Einstein and Green-Kubo equation
3.   Hydrogen bond analysis
4.   Rotational autocorrelation function
5.   NMR NOE calculation
6.   Water residence time
7.   Cluster analysis based on linkage method
8.   Aggregation volume calculation for specific molecule
9.   IR & Raman Spectrum calculation
10.  Density of State (DOS) Spectrum
11.  Hydrogen bond lifetime analysis
12.  Various angle related analysis and can convert probability distribution to Gibbs free energy plot
13.  Orientation-Resolved Radial Distribution Functions (JCTC 2019, 15, 803−812)
14.  Conditional Time Correlation Function (JCTC 2019, 15, 803−812)
15.  Radius of gyration (mass-weighted)
16.  Biphase Mix Index
17.  Co-plane Index (CPI)
18.  RMSD & RMSF with outputting superposed structures
19.  Use AmberMask for selecting Residues, Atoms in topology structure
20.  Use Domain Specific Language(DSL) script to drive analysis process 
21.  Other practical utilities for using Tinker

#### Typical QM analysis functions

1.   NBO spin summary
2.   NAO orbital contribution by driving Multiwfn
3.   Delocalization Index (DI) by driving Multiwfn
4.   QTAIM analysis by driving Multiwfn
5.   ADCH Charge

#### Build Requirements
- Language : C++20 ( GCC 9.x or above )
- Build System :  CMake 3.15
- Third-party libraries :  Boost 1.72,
                           Intel Threading Building Blocks(TBB) 2020.1 ( for Multi-Core Parallelism ), 
                           NetCDF,  FFTW3, 
                           GROMACS library 5.1.4 ( for reading and writting Gromacs topology and trajectory file ),
                           Google Test 1.10 ( for unit test ),
                           Eigen 3.3 ( C++ template library for linear algebra ),
                           Pugixml 1.10 ( XML format process ),
                           nlohmann/json 3.7.3,
                           GNU readline
                           
                          
                           
#### Execution Flow
<img src="doc/ExecutionFlow.png" alt="Software Execution Flow" title="Software Execution Flow" width="700" height="285" />


                           

 
