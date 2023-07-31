#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H

#include "print_data.h"
#include "modified_voronoi.h"

void test_all(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& V_with_indices, Triangulation& T, ModifiedVoronoi& tessellation);

#endif // TEST_FUNCTION_H
