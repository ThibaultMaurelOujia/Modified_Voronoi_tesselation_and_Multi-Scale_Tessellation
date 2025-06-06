# Modified Voronoi and Multi-Scale Tessellation

Modified Voronoi tessellation and Multi-Scale Tessellation for computing differential operators of particle velocity.



## Description

The objective of this project is to implement a method for computing differential operators based on the variation of the volume of modified Voronoi tessellations, as described by **Maurel–Oujia, Matsuda, et al. (2024)**. Additionally, it aims to port the Python code developed during the CTR Summer Program 2022, used for multiresolution analysis as described by **Matsuda, Schneider, Oujia, et al. (2022)**. A significant portion of this code was implemented during the JSPS Short-Term Program in 2023 at JAMSTEC in Japan.


## Licensing and Contributions

The code is managed under the **GNU General Public License v3.0**. Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed. Further details can be found in the LICENSE file.

The original version of the C++ code was developed by **Thibault Maurel Oujia** at **I2M, Aix-Marseille University** to compute the volume cells in parallel using OpenMP for shared memory and to compute the divergence and curl of the particle velocity.

During the JSPS Short-Term Program in 2023, Thibault Maurel Oujia implemented parallel computation of the volume using MPI for distributed memory and wavelets on graphs using MPI. The implementation of wavelets on graphs is a C++ version of the code originally developed in Python by **Keigo Matsuda** and Thibault Maurel Oujia during the CTR Summer Program 2022.


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



## Compilation and Execution

To compile and run the code in **Release** mode, follow these steps:

### Step 1: Clone the Repository

First, clone the repository to your local machine using Git:

```bash
git clone git@github.com:ThibaultMaurelOujia/Modified_Voronoi_tesselation_and_Multi-Scale_Tessellation.git
```

### Step 2: Build the Project

Navigate to the project directory and create a new build directory:

```bash
cd Modified_Voronoi_tesselation_and_Multi-Scale_Tessellation
rm -fr build
mkdir build
cd build
```

Configure the project with CMake for a Release build:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

Compile the project:

```bash
cmake --build . --target clean
cmake --build . --config Release
```

### Step 3: Run the Application

After building the project, you can run the executable in different ways:

- **Local Execution:** `./parallel_delaunay_multiscale`
- **Using MPI on multiple processors:** `mpiexec -n 8 ./parallel_delaunay_multiscale ../config_example.txt`
- **Rebuild and Run for Quick Testing:** `clear && cmake --build . --config Release && mpiexec -n 1 ./parallel_delaunay_multiscale ../config_example.txt`

Replace `config_example.txt` with your specific configuration file as needed.











# References

- Thibault Maurel Oujia, *On the Particle Dynamics in Fully Developed Turbulence: Tessellation, Multiresolution and Machine Learning Methods*. Ph.D. Thesis.

- Keigo Matsuda, Kai Schneider, Thibault Maurel Oujia, Jacob West, Suhas Jain, and Kazuki Maeda, "Multiresolution analysis of inertial particle tessellations for clustering dynamics," in *Proceedings of the Summer Program*, Center for Turbulence Research, Stanford University, 2022.

- Thibault Maurel–Oujia, Keigo Matsuda, and Kai Schneider, "Computing differential operators of the particle velocity in moving particle clouds using tessellations," *Journal of Computational Physics*, vol. 498, pp. 112658, 2024.













