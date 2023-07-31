#include "modified_voronoi.h"



void ModifiedVoronoi::create_empty_cells(size_t num_cells) {
    std::cout << "CELLS ALLOC:      " << std::flush;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_cells; ++i) {
        Cell empty_cell;
        cells.push_back(empty_cell); // !!!!!!!! changer ca pour d'abord defiler pour allouer
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|0|10|20|30|40|50|60|70|80|90|100|  Time: " << duration/1000. << " s" << std::endl;
}


void ModifiedVoronoi::create_empty_vertices(Triangulation& T) {
    std::cout << "VERTICES ALLOC:   " << std::flush;
    auto start_time = std::chrono::high_resolution_clock::now();
    for (auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end(); ++cell_itr) {
        Point empty_vertex;
        vertices.push_back(empty_vertex); // !!!!!!!! changer ca pour d'abord defiler pour allouer
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|0|10|20|30|40|50|60|70|80|90|100|  Time: " << duration/1000. << " s" << std::endl;
}


std::vector<double> ModifiedVoronoi::get_cells_volumes(int Np) const {
    std::vector<double> volumes;

    if (Np == -1) {
        Np = cells.size();
    }

    volumes.reserve(Np);

    int count = 0;
    for (const auto& cell : cells) {
        if (count < Np) {
            volumes.push_back(cell.volume);
            ++count;
        }
        else {
            break;
        }
    }

    return volumes;
}


void ModifiedVoronoi::print_cells(int cell_index) const {
    // tessellation.print_cells(0);
    if (cell_index == -1) {
        // Imprimer toutes les cellules
        for (size_t i = 0; i < cells.size(); ++i) {
            const Cell& cell = cells[i];
            std::cout << "Cell " << i << ":" << std::endl;
            std::cout << "  Cell indices vertices: ";
            for (size_t index : cell.cell_indices_vertices) {
                std::cout << index << " ";
            }
            std::cout << std::endl;

            std::cout << "  Cell facets indices from cell indices vertices: " << std::endl;
            for (const auto& facet : cell.cell_facets_indices_from_cell_indices_vertices) {
                std::cout << "    ";
                for (size_t index : facet) {
                    std::cout << index << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "  Volume: " << cell.volume << std::endl;
        }
    }
    else if (cell_index >= 0 && cell_index < cells.size()) {
        // Imprimer la cellule sp�cifi�e
        const Cell& cell = cells[cell_index];
        std::cout << "Cell " << cell_index << ":" << std::endl;
        std::cout << "  Cell indices vertices: ";
        for (size_t index : cell.cell_indices_vertices) {
            std::cout << index << " ";
        }
        std::cout << std::endl;

        std::cout << "  Cell facets indices from cell indices vertices: " << std::endl;
        for (const auto& facet : cell.cell_facets_indices_from_cell_indices_vertices) {
            std::cout << "    ";
            for (size_t index : facet) {
                std::cout << index << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "  Volume: " << cell.volume << std::endl;
    }
    else {
        std::cerr << "Invalid cell index" << std::endl;
    }
}


// create_vertices_modified_voronoi_cells
void ModifiedVoronoi::store_gravity_centers_and_indices(Triangulation& T, const std::vector<std::pair<CGAL::Point_3<K>, std::pair<std::size_t, std::size_t>>>& P_with_indices) {

    std::cout << "COMPUTE VERTICES: " << std::flush;
    auto start_time = std::chrono::high_resolution_clock::now();

    size_t vertex_index = 0;
    for (auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end(); ++cell_itr) {
        std::array<std::size_t, 4> indices;
        std::array<CGAL::Point_3<K>, 4> points;
        for (int i = 0; i < 4; ++i) {
            std::pair<std::size_t, std::size_t> idx = cell_itr->vertex(i)->info();
            indices[i] = idx.first;
            points[i] = P_with_indices[indices[i]].first;
        }

        // Calculate the center of gravity manually
        double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
        for (const auto& point : points) {
            x_sum += point.x();
            y_sum += point.y();
            z_sum += point.z();
        }

        ModifiedVoronoi::Point gravity_center_mv(x_sum / 4.0, y_sum / 4.0, z_sum / 4.0);

        // Update the gravity center in vertices
        vertices[vertex_index] = gravity_center_mv;

        // Add the gravity center index to each cell in ModifiedVoronoi
        for (std::size_t index : indices) {
            cells[index].cell_indices_vertices.push_back(vertex_index);
        }

        vertex_index++;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|0|10|20|30|40|50|60|70|80|90|100|  Time: " << duration/1000. << " s" << std::endl;
}


void ModifiedVoronoi::create_facet_modified_voronoi_cells(bool loading_bar, int Np) {
    // Calculer l'enveloppe convexe et le volume pour toutes les cellules
    int num_cells = static_cast<int>(cells.size());

    if (Np != -1) {
        num_cells = Np;
    }

    std::atomic<int> progress_count(0);
    std::vector<bool> progress_markers(11, false);

    auto start_time = std::chrono::high_resolution_clock::now();
    if (loading_bar) {
        std::cout << "CREATE FACET:     |0" << std::flush;
        progress_markers[0] = true;
    }


    #pragma omp parallel
    {
        // Copie et normalisation des points de cell_vertices
        std::vector<CGAL::Point_3<K>> cell_vertices;
        std::vector<double> points_3D;
        std::vector<size_t> original_indices;

        cell_vertices.reserve(20);
        points_3D.reserve(3 * 20);
        original_indices.reserve(4);
        
        #pragma omp for
        for (int i = 0; i < num_cells; ++i) {
            Cell& cell = cells[i];

            cell_vertices.clear();
            points_3D.clear();

            for (size_t idx : cell.cell_indices_vertices) {
                const Point& point = vertices[idx];
                cell_vertices.push_back(CGAL::Point_3<K>(point.x, point.y, point.z));
            }

            // Calcul du centre de gravite
            double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
            for (const auto& vertex : cell_vertices) {
                x_sum += vertex.x();
                y_sum += vertex.y();
                z_sum += vertex.z();
            }
            size_t num_vertices = cell_vertices.size();
            CGAL::Point_3<K> center_of_gravity(x_sum / num_vertices, y_sum / num_vertices, z_sum / num_vertices);

            // Normalisation des points
            for (const auto& vertex : cell_vertices) {
                CGAL::Vector_3<K> centered_vertex = vertex - center_of_gravity;
                double norm = std::sqrt(centered_vertex.squared_length());
                CGAL::Point_3<K> normed_vertex = center_of_gravity + (centered_vertex / norm);
                points_3D.push_back(normed_vertex.x());
                points_3D.push_back(normed_vertex.y());
                points_3D.push_back(normed_vertex.z());
            }

            // Calcul de l'enveloppe convexe avec Qhull
            int ndim = 3;
            std::string comment = "";
            std::string qhull_command = "Qt";

            try {
                orgQhull::Qhull qhull(comment.c_str(), ndim, num_vertices, points_3D.data(), qhull_command.c_str());

                // Calcul du volume en parcourant les faces de l'enveloppe convexe des points normalis�s
                orgQhull::QhullFacetList facets = qhull.facetList();
                for (const auto& facet : facets) {
                    original_indices.clear();
                    orgQhull::QhullVertexSet vertices = facet.vertices();
                    for (const auto& vertex : vertices) {
                        int original_index = vertex.point().id();
                        //original_indices.push_back(cell.cell_indices_vertices[original_index]);
                        original_indices.push_back(original_index);
                    }

                    std::array<uint16_t, 3> facet_indices = { static_cast<uint16_t>(original_indices[0]), static_cast<uint16_t>(original_indices[1]), static_cast<uint16_t>(original_indices[2]) };
                    cell.cell_facets_indices_from_cell_indices_vertices.push_back(facet_indices);
                }
            }
            catch (orgQhull::QhullError& e) {
                std::cerr << e.what() << std::endl;
                throw std::runtime_error("QhullError in compute_volume_cell");
            }

            size_t count = ++progress_count;
            if (loading_bar) {
                size_t progress = count * 100 / num_cells;
                size_t index = progress / 10;
                if (index < progress_markers.size()) {
                    if (!progress_markers[index]) {
                        #pragma omp critical
                        {
                            if (!progress_markers[index]) {
                                std::cout << "|" << progress << std::flush;
                                progress_markers[index] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    if (loading_bar) {
        std::cout << "|  Time: " << duration/1000. << " s" << std::endl;
    }
}


void ModifiedVoronoi::compute_volume_cell_exact(Cell& cell) {
    // Copie des points de cell_vertices
    std::vector<CGAL::Point_3<K>> cell_vertices;
    for (size_t idx : cell.cell_indices_vertices) {
        const Point& point = vertices[idx];
        cell_vertices.push_back(CGAL::Point_3<K>(point.x, point.y, point.z));
    }

    // Calcul du centre de gravit�
    double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
    for (const auto& vertex : cell_vertices) {
        x_sum += vertex.x();
        y_sum += vertex.y();
        z_sum += vertex.z();
    }
    size_t num_vertices = cell_vertices.size();
    CGAL::Point_3<K> center_of_gravity(x_sum / num_vertices, y_sum / num_vertices, z_sum / num_vertices);

    // Calcul du volume en utilisant les indices des facettes stock�s dans cell_facets_indices_from_cell_indices_vertices
    double total_volume = 0.0;
    for (const auto& facet_indices : cell.cell_facets_indices_from_cell_indices_vertices) {
        CGAL::Point_3<K> p1 = cell_vertices[facet_indices[0]];
        CGAL::Point_3<K> p2 = cell_vertices[facet_indices[1]];
        CGAL::Point_3<K> p3 = cell_vertices[facet_indices[2]];

        CGAL::Vector_3<K> v1 = p2 - p1;
        CGAL::Vector_3<K> v2 = p3 - p1;
        CGAL::Vector_3<K> normal = CGAL::cross_product(v1, v2);

        // Calcul du volume du t�tra�dre form� par la face et le centre de gravit�
        double volume = std::abs(normal * (center_of_gravity - p1)) / 6.0;
        // Ajout du volume au volume total
        total_volume += volume;
    }
    // Mise � jour du volume de la cellule
    cell.volume = total_volume;
}


void ModifiedVoronoi::compute_volume_convex_hull(Cell& cell) {
    std::vector<CGAL::Point_3<K>> cell_vertices;
    CGAL::Polyhedron_3<K> convex_hull;
    for (size_t idx : cell.cell_indices_vertices) {
        const Point& point = vertices[idx];
        cell_vertices.push_back(CGAL::Point_3<K>(point.x, point.y, point.z));
    }
    CGAL::convex_hull_3(cell_vertices.begin(), cell_vertices.end(), convex_hull);

    if (convex_hull.size_of_vertices() < 4) {
        cell.volume = 0.0;
        return;
    }

    // 1. Calculate the center of gravity of the convex hull
    double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
    for (auto vertex_itr = convex_hull.vertices_begin(); vertex_itr != convex_hull.vertices_end(); ++vertex_itr) {
        const CGAL::Point_3<K>& point = vertex_itr->point();
        x_sum += point.x();
        y_sum += point.y();
        z_sum += point.z();
    }
    size_t num_vertices = convex_hull.size_of_vertices();
    CGAL::Point_3<K> center_of_gravity(x_sum / num_vertices, y_sum / num_vertices, z_sum / num_vertices);


    // 2-4. Calculate the volume for each face
    double total_volume = 0.0;
    for (auto facet_itr = convex_hull.facets_begin(); facet_itr != convex_hull.facets_end(); ++facet_itr) {
        auto halfedge_itr = facet_itr->facet_begin();
        CGAL::Point_3<K> p1 = halfedge_itr->vertex()->point();
        ++halfedge_itr;
        CGAL::Point_3<K> p2 = halfedge_itr->vertex()->point();
        ++halfedge_itr;
        CGAL::Point_3<K> p3 = halfedge_itr->vertex()->point();

        CGAL::Vector_3<K> v1 = p2 - p1;
        CGAL::Vector_3<K> v2 = p3 - p1;
        CGAL::Vector_3<K> normal = CGAL::cross_product(v1, v2);

        // Calculate the volume of the tetrahedron formed by the face and the center of gravity
        double volume = std::abs(normal * (center_of_gravity - p1)) / 6.0;

        // 3. Add the volume to the total volume
        total_volume += volume;
    }

    // 5. Update the cell volume
    cell.volume = total_volume;
}


void ModifiedVoronoi::compute_volume_modified_voronoi_cells(int cell_index, bool loading_bar, bool exact_volume, int Np) {
    if (cell_index == -1) {
        // Calculer l'enveloppe convexe et le volume pour toutes les cellules
        int num_cells = static_cast<int>(cells.size());

        if (Np != -1) {
            num_cells = Np;
        }

        std::atomic<int> progress_count(0);
        std::vector<bool> progress_markers(11, false);

        auto start_time = std::chrono::high_resolution_clock::now();
        if (loading_bar) {
            std::cout << "COMPUTE VOLUME:   |0";
            progress_markers[0] = true;
        }

        #pragma omp parallel for
        for (int i = 0; i < num_cells; ++i) {
            if (exact_volume) {
                compute_volume_cell_exact(cells[i]);
            }
            else {
                compute_volume_convex_hull(cells[i]);
            }
            size_t count = ++progress_count;
            if (loading_bar) {
                size_t progress = count * 100 / num_cells;
                size_t index = progress / 10;
                if (index < progress_markers.size()) {
                    if (!progress_markers[index]) {
                        #pragma omp critical
                        {
                            if (!progress_markers[index]) {
                                std::cout << "|" << progress << std::flush;
                                progress_markers[index] = true;
                            }
                        }
                    }
                }
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        if (loading_bar) {
            std::cout << "|  Time: " << duration / 1000. << " s" << std::endl;
        }
    }
    else {
        // Calculer l'enveloppe convexe et le volume pour une cellule sp�cifique
        if (static_cast<size_t>(cell_index) >= cells.size()) {
            throw std::out_of_range("Cell index out of range");
        }
        if (exact_volume) {
            compute_volume_cell_exact(cells[cell_index]);
        }
        else {
            compute_volume_convex_hull(cells[cell_index]);
        }
    }
}









