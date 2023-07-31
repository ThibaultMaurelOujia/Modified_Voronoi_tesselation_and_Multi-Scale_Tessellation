#ifndef MODIFIED_VORONOI_H
#define MODIFIED_VORONOI_H


#include <omp.h>
#include <vector>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>


#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullPoint.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>


#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_data_structure_3<
    // CGAL::Triangulation_vertex_base_with_info_3<std::size_t, K>,
    CGAL::Triangulation_vertex_base_with_info_3<std::pair<std::size_t, std::size_t>, K>,
    CGAL::Delaunay_triangulation_cell_base_3<K>,
    CGAL::Parallel_tag>                               Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>              Triangulation;
typedef Triangulation::Point                                Point;
typedef K::Point_3 CgalPoint;


class ModifiedVoronoi {
public:
    class Point {
    public:
        double x, y, z;

        Point() : x(0), y(0), z(0) {}
        Point(double x, double y, double z) : x(x), y(y), z(z) {}
    };

    class Cell {
    public:
        std::vector<size_t> cell_indices_vertices;
        std::vector<std::array<uint16_t, 3>> cell_facets_indices_from_cell_indices_vertices;
        double volume;
    };


    void print_cells(int cell_index = -1) const;

    void create_empty_cells(size_t num_cells);
    void create_empty_vertices(Triangulation& T);
    std::vector<double> get_cells_volumes(int Np = -1) const;

    void store_gravity_centers_and_indices(Triangulation& T, const std::vector<std::pair<CGAL::Point_3<K>, std::pair<std::size_t, std::size_t>>>& P_with_indices);

    // Create facet of the modified Voronoi tessellation
    void create_facet_modified_voronoi_cells(bool loading_bar = false, int Np = -1);

    void compute_volume_cell_exact(Cell& cell);
    void compute_volume_convex_hull(Cell& cell);

    void compute_volume_modified_voronoi_cells(int cell_index = -1, bool loading_bar = true, bool exact_volume = true, int Np = -1);
    

private:
    std::vector<Point> vertices;
    std::vector<Cell> cells;
    // std::deque<Point> vertices;
    // std::deque<Cell> cells;
};


#endif // MODIFIED_VORONOI_H


