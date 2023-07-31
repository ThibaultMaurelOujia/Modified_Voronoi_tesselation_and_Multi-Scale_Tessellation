#include "print_data.h"




void print_points(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& V_with_indices, int max_count) {
    int count = 0;

    std::cout << "Total number of points: " << V_with_indices.size() << std::endl;
    count = 0;
    for (const auto& point_with_index : V_with_indices) {
        if (count >= max_count) {
            break;
        }
        const Point& p = point_with_index.first;
        std::size_t index = point_with_index.second.first;
        // std::cout << "Point: " << p << ", Index: " << index << std::endl;
        std::cout << "P[" << index << "]: " << p << std::endl; 
        ++count;
    }
    std::cout << "_____________________________" << std::endl;
}


void print_vertices(Triangulation& T, int max_count) {
    int count = 0;

    // Pour parcourir les sommets (auto vertex_itr = T.finite_vertices_begin(); vertex_itr != T.finite_vertices_end(); ++vertex_itr) 
    std::cout << "Total number of vertices: " << T.number_of_vertices() << std::endl;
    count = 0;
    for (auto vertex_itr = T.finite_vertices_begin(); vertex_itr != T.finite_vertices_end() && count < max_count; ++vertex_itr, ++count) {
        Point p = vertex_itr->point();
        std::cout << "Vertex: " << p << std::endl;
    }
    std::cout << "_____________________________" << std::endl;
}


void print_vertices_with_indices(Triangulation& T, int max_count) {
    int count = 0;

    std::cout << "Total number of vertices: " << T.number_of_vertices() << std::endl;
    count = 0;
    for (auto vertex_itr = T.finite_vertices_begin(); vertex_itr != T.finite_vertices_end() && count < max_count; ++vertex_itr, ++count) {
        Point p = vertex_itr->point();
        std::pair<std::size_t, std::size_t> idx = vertex_itr->info();
        std::cout << "Vertex: " << p << " (Index: " << idx.first << ")" << std::endl;
    }
    std::cout << "_____________________________" << std::endl;
}


void print_edges(Triangulation& T, int max_count) {
    int count = 0;

    // Pour parcourir les ar�tes (auto edge_itr = T.finite_edges_begin(); edge_itr != T.finite_edges_end(); ++edge_itr)
    std::cout << "Total number of finite edges: " << T.number_of_finite_edges() << std::endl;
    count = 0;
    for (auto edge_itr = T.finite_edges_begin(); edge_itr != T.finite_edges_end() && count < max_count; ++edge_itr, ++count) {
        auto v1 = edge_itr->first->vertex(T.vertex_triple_index(edge_itr->second, 0))->point();
        auto v2 = edge_itr->first->vertex(T.vertex_triple_index(edge_itr->second, 1))->point();
        std::cout << "Edge: " << v1 << " - " << v2 << std::endl;
    }
    std::cout << "_____________________________" << std::endl;
}


void print_edges_with_indices(Triangulation& T, const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& V_with_indices, int max_count) {
    int count = 0;

    std::cout << "Total number of finite edges: " << T.number_of_finite_edges() << std::endl;
    count = 0;
    for (auto edge_itr = T.finite_edges_begin(); edge_itr != T.finite_edges_end() && count < max_count; ++edge_itr, ++count) {
        Triangulation::Edge e = *edge_itr;
        Triangulation::Vertex_handle v1 = e.first->vertex(T.vertex_triple_index(e.second, 0));
        Triangulation::Vertex_handle v2 = e.first->vertex(T.vertex_triple_index(e.second, 1));
        std::pair<std::size_t, std::size_t> idx1 = v1->info();
        std::pair<std::size_t, std::size_t> idx2 = v2->info();

        std::cout << "Edge (indices): " << idx1.first << " | " << idx2.first << std::endl;
        std::cout << "Edge: "
            << V_with_indices[idx1.first].first << " (Index: " << V_with_indices[idx1.first].second.first << ") | "
            << V_with_indices[idx2.first].first << " (Index: " << V_with_indices[idx2.first].second.first << ")"
            << std::endl;
    }
    std::cout << "_____________________________" << std::endl;
}



void print_cells(Triangulation& T, int max_count) {
    int count = 0;

    // Pour parcourir les cellules (t�tra�dres) (auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end(); ++cell_itr)
    std::cout << "Total number of finite cells: " << T.number_of_finite_cells() << std::endl;
    count = 0;
    for (auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end() && count < max_count; ++cell_itr, ++count) {
        Point p1 = cell_itr->vertex(0)->point();
        Point p2 = cell_itr->vertex(1)->point();
        Point p3 = cell_itr->vertex(2)->point();
        Point p4 = cell_itr->vertex(3)->point();
        std::cout << "Cell (tetrahedron): " << p1 << " | " << p2 << " | " << p3 << " | " << p4 << std::endl;
    }
    std::cout << "_____________________________" << std::endl;
}


void print_cells_with_indices(Triangulation& T, int max_count) {
    int count = 0;

    // Pour parcourir les indices qui composent les cellules (t�tra�dres) (auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end(); ++cell_itr)
    std::cout << "Total number of finite cells: " << T.number_of_finite_cells() << std::endl;
    count = 0;
    for (auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end() && count < max_count; ++cell_itr, ++count) {
        std::pair<std::size_t, std::size_t> idx1 = cell_itr->vertex(0)->info();
        std::pair<std::size_t, std::size_t> idx2 = cell_itr->vertex(1)->info();
        std::pair<std::size_t, std::size_t> idx3 = cell_itr->vertex(2)->info();
        std::pair<std::size_t, std::size_t> idx4 = cell_itr->vertex(3)->info();
        std::cout << "Cell (indices): " << idx1.first << " | " << idx2.first << " | " << idx3.first << " | " << idx4.first << std::endl;
        Point p1 = cell_itr->vertex(0)->point();
        Point p2 = cell_itr->vertex(1)->point();
        Point p3 = cell_itr->vertex(2)->point();
        Point p4 = cell_itr->vertex(3)->point();
        std::cout << "Cell (tetrahedron): " << p1 << " | " << p2 << " | " << p3 << " | " << p4 << std::endl;
    }
    std::cout << "_____________________________" << std::endl;
}










