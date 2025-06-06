#ifndef PARALLEL_DELAUNAY_H
#define PARALLEL_DELAUNAY_H


#include <mpi.h>
#include <vector>
#include <fstream>
#include <iostream>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>



#include <tbb/global_control.h>


#include "utile.h"
#include "print_data.h"
#include "test_function.h"
#include "modified_voronoi.h"
#include "multiscale_tessellation.h"


bool parallel_delaunay_multiscale(const ArgConfig& ARG_CONFIG, int subdomainNumber, int totalSubdomains);



bool updateGlobalIndicesGlobalProcess(std::vector<int16_t>& GlobalIndices_GlobalProcess,
                                      const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized_subcube,
                                      std::size_t Np, std::size_t Np_Domain, int subdomainNumber, int MPI_rank, int MPI_size, const ArgConfig& ARG_CONFIG);
void updatePeriodizedSubcube(std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized_subcube,
                             std::vector<Point>& V_with_indices_periodized_subcube,
                             const ArgConfig& ARG_CONFIG, std::size_t Np_Domain);
void ifLOAD_GRAPH(const ArgConfig& ARG_CONFIG);






#endif // PARALLEL_DELAUNAY_H








