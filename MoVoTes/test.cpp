/*






find_package(UPCXX REQUIRED)
target_link_libraries(<your_target> ${UPCXX_LIBRARIES})
target_include_directories(<your_target> PUBLIC ${UPCXX_INCLUDE_DIRS})




#include <iostream>
#include <upcxx/upcxx.hpp>

int main(int argc, char *argv[]) {
    upcxx::init(); // Initialiser UPC++

    int rank = upcxx::rank_me();
    int total_ranks = upcxx::rank_n();

    std::cout << "Hello World from rank " << rank << " out of " << total_ranks << " ranks." << std::endl;

    upcxx::finalize(); // Terminer UPC++

    return 0;
}



















#include <iostream>
#include <upcxx/upcxx.hpp>

constexpr size_t ARRAY_SIZE = 200ULL * 1024ULL * 1024ULL * 1024ULL; // 200 Go

int main(int argc, char *argv[]) {
    upcxx::init();

    if (upcxx::rank_me() == 0) {
        std::cout << "Allocating a distributed array of size 200 GB..." << std::endl;
    }

    // Allouer un tableau distribué de 200 Go
    upcxx::global_ptr<int64_t> global_array = upcxx::new_array<int64_t>(ARRAY_SIZE / sizeof(int64_t));

    // ...

    upcxx::finalize();
    return 0;
}




// Lecture d'un élément à distance
upcxx::future<int64_t> read_future = upcxx::rget(global_array + remote_index);
int64_t remote_value = read_future.wait();

// Écriture d'un élément à distance
int64_t new_value = 42;
upcxx::future<> write_future = upcxx::rput(new_value, global_array + remote_index);
write_future.wait();

// Opération atomique à distance (ajout)
int64_t add_value = 5;
upcxx::future<int64_t> atomic_future = upcxx::atomic_fetch_add(global_array + remote_index, add_value);
int64_t old_value = atomic_future.wait();











upcxx::init();

// Votre code de définition des types et de génération des points ici...

if (upcxx::rank_me() == 0) {
    // Exécuter le code pour la première triangulation sur le nœud 0
    Triangulation T1(P_with_indices.begin(), P_with_indices.end(), &locking_ds);
} else if (upcxx::rank_me() == 1) {
    // Exécuter le code pour la deuxième triangulation sur le nœud 1
    Triangulation T2(P_with_indices.begin(), P_with_indices.end(), &locking_ds);
}

upcxx::finalize();


















#include <upcxx/upcxx.hpp>
#include <vector>
// Inclure les headers nécessaires pour CGAL et les autres structures de données

int main() {
    upcxx::init();

    // Votre code de définition des types et de génération des points ici...

    int num_nodes = upcxx::rank_n(); // Obtenir le nombre total de nœuds
    std::vector<Triangulation*> triangulations(num_nodes, nullptr); // Créer un tableau de pointeurs de triangulation

    // Créer une triangulation pour chaque nœud et stocker l'adresse dans le tableau
    for (int i = 0; i < num_nodes; i++) {
        if (upcxx::rank_me() == i) {
            Triangulation* T = new Triangulation(P_with_indices.begin(), P_with_indices.end(), &locking_ds);
            triangulations[i] = T;
        }
    }

    // Effectuer les calculs nécessaires sur les triangulations
    // ...

    // Nettoyer la mémoire allouée pour les triangulations
    for (Triangulation* T : triangulations) {
        if (T) {
            delete T;
        }
    }

    upcxx::finalize();
    return 0;
}

















#include <upcxx/upcxx.hpp>
#include <vector>
// Inclure les headers nécessaires pour CGAL et les autres structures de données

// Fonction qui sera appelée à distance pour obtenir la position du premier point de la première cellule
upcxx::future<Point> get_first_point_position(Triangulation* T) {
    auto cell_itr = T->finite_cells_begin();
    Point first_point_position = cell_itr->vertex(0)->point();
    return upcxx::make_future(first_point_position);
}

int main() {
    upcxx::init();

    // Votre code de définition des types et de génération des points ici...
    // Votre code pour créer les triangulations ici...

    if (upcxx::rank_me() == 0) {
        // Le nœud 0 demande l'information au nœud 1 (supposant que les données sont sur le nœud 1)
        int target_node = 1;
        Triangulation* target_triangulation = triangulations[target_node];
        upcxx::future<Point> point_future = upcxx::rpc(target_node, get_first_point_position, target_triangulation);

        // Attendre que la réponse soit disponible et afficher la position du point
        point_future.wait();
        Point first_point_position = point_future.result();
        std::cout << "Position du premier point de la première cellule sur le nœud 1 : " << first_point_position << std::endl;
    }

    // Votre code pour nettoyer la mémoire allouée pour les triangulations ici...

    upcxx::finalize();
    return 0;
}


























#include <iostream>
#include <Eigen/Dense>

Eigen::Matrix3d rotationMatrix(const Eigen::Vector3d& axis, double angle) {
    Eigen::Matrix3d rotation;
    double c = cos(angle);
    double s = sin(angle);
    double one_minus_c = 1.0 - c;
    double x = axis.x();
    double y = axis.y();
    double z = axis.z();

    rotation(0, 0) = x * x * one_minus_c + c;
    rotation(0, 1) = x * y * one_minus_c - z * s;
    rotation(0, 2) = x * z * one_minus_c + y * s;

    rotation(1, 0) = y * x * one_minus_c + z * s;
    rotation(1, 1) = y * y * one_minus_c + c;
    rotation(1, 2) = y * z * one_minus_c - x * s;

    rotation(2, 0) = z * x * one_minus_c - y * s;
    rotation(2, 1) = z * y * one_minus_c + x * s;
    rotation(2, 2) = z * z * one_minus_c + c;

    return rotation;
}

Eigen::Matrix3d alignNormalToZAxis(const Eigen::Vector3d& normal) {
    Eigen::Vector3d z_axis(0, 0, 1);

    // Find the axis of rotation (cross product of normal and z_axis)
    Eigen::Vector3d axis_of_rotation = normal.cross(z_axis).normalized();

    // Find the angle of rotation (dot product of normal and z_axis)
    double angle_of_rotation = acos(normal.dot(z_axis));

    // Calculate the rotation matrix
    Eigen::Matrix3d rotation = rotationMatrix(axis_of_rotation, angle_of_rotation);

    return rotation;
}

int main() {
    // Define the normal vector of the face
    Eigen::Vector3d normal(1, 1, 1);
    normal.normalize();

    // Calculate the rotation matrix to align the normal with the z-axis
    Eigen::Matrix3d rotation = alignNormalToZAxis(normal);

    // Apply the rotation to the points of the face
    Eigen::Vector3d point(1, 2, 3);
    Eigen::Vector3d rotated_point = rotation * point;

    std::cout << "Rotated point: " << rotated_point.transpose() << std::endl;

    return 0;
}











*/