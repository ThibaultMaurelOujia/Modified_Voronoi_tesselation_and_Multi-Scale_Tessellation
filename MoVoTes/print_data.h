#ifndef PRINT_DATA_H
#define PRINT_DATA_H


#include <vector>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// Delaunay T3
typedef CGAL::Triangulation_data_structure_3<
    // CGAL::Triangulation_vertex_base_3<K>,
    // CGAL::Triangulation_vertex_base_with_info_3<std::size_t, K>,
    CGAL::Triangulation_vertex_base_with_info_3<std::pair<std::size_t, std::size_t>, K>,
    CGAL::Delaunay_triangulation_cell_base_3<K>,
    CGAL::Parallel_tag>                               Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>              Triangulation;
typedef K::Point_3 Point;

void print_points(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& V_with_indices, int max_count);
void print_vertices(Triangulation& T, int max_count);
void print_vertices_with_indices(Triangulation& T, int max_count);
void print_edges(Triangulation& T, int max_count);
void print_edges_with_indices(Triangulation& T, const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& V_with_indices, int max_count);
void print_cells(Triangulation& T, int max_count);
void print_cells_with_indices(Triangulation& T, int max_count);



#endif // PRINT_DATA_H





