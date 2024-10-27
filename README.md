# Modified Voronoi and Multi-Scale Tessellation

Modified Voronoi tessellation and Multi-Scale Tessellation for computing differential operators of particle velocity.



## Description

The objective of this project is to implement a method for computing differential operators based on the variation of the volume of modified Voronoi tessellations, as described by **Maurel–Oujia, Matsuda, et al. (2024)**. Additionally, it aims to port the Python code developed during the CTR Summer Program 2022, used for multiresolution analysis as described by **Matsuda, Schneider, Oujia, et al. (2022)**. A significant portion of this code was implemented during the JSPS Short-Term Program in 2023 at JAMSTEC in Japan.


## Licensing and Contributions

The code is managed under the **GNU General Public License v3.0**. Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed. Further details can be found in the LICENSE file.

The original version of the C++ code was developed by **Thibault Maurel Oujia** at **I2M, Aix-Marseille University** to compute the volume cells in parallel using OpenMP for shared memory and to compute the divergence and curl of the particle velocity.

During the JSPS Short-Term Program in 2023, Thibault implemented parallel computation of the volume using MPI for distributed memory and wavelets on graphs using MPI. The implementation of wavelets on graphs is a C++ version of the code originally developed in Python by **Keigo Matsuda** and Thibault Maurel Oujia during the CTR Summer Program 2022.


## Prerequisites

This code is developed to be built and run on **Linux** and **macOS** machines. The following dependencies are required:

- **C++ Compiler**: Supporting C++17 or higher (e.g., GCC 9.0+, Clang 10.0+).
- **CMake**: Version 3.1 or higher.
- **CGAL**: Computational Geometry Algorithms Library.
- **TBB**: Threading Building Blocks for parallelism.
- **MPI**: Message Passing Interface for distributed memory parallelism.
- **GMP**: GNU Multiple Precision Arithmetic Library.
- **Qhull**: For computing convex hulls.
- **Boost**: Boost C++ Libraries.
- **OpenMP**: For shared-memory parallelism.



## References

- Thibault Maurel Oujia, “On the Particle Dynamics in Fully Developed Turbulence: Tessellation, Multiresolution and Machine Learning Methods.”
- K. Matsuda et al., “Multiresolution analysis of inertial particle tessellations for clustering dynamics,” Proceedings of the Summer Program, Center for Turbulence Research, Stanford University, 2022.
- Thibault Maurel–Oujia, Keigo Matsuda, and Kai Schneider, “Computing differential operators of the particle velocity in moving particle clouds using tessellations,” Journal of Computational Physics, vol. 498, 112658, 2024.















