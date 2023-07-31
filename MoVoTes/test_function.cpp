#include <iostream>

#include "print_data.h"
#include "test_function.h"
#include "modified_voronoi.h"


void test_all(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& V_with_indices, Triangulation& T, ModifiedVoronoi& tessellation) {


    /*******************************************************************
    *					   Delaunay tessellation					   *
    *******************************************************************/


    print_points(V_with_indices, 3);
    print_vertices(T, 3);
    print_vertices_with_indices(T, 3);
    print_edges(T, 3);
    print_edges_with_indices(T, V_with_indices, 3);
    print_cells(T, 3);
    print_cells_with_indices(T, 3);


    /*******************************************************************
    *                    Modified Voronoi tessellation                 *
    *******************************************************************/



    tessellation.print_cells(0);

    
    tessellation.print_cells();













}

