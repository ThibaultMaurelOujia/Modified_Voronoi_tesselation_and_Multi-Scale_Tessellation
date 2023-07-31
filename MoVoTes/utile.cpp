#include <fstream>
#include <vector>
#include <iostream>


#include "utile.h"



/*******************************************************************
*					       CONFIGURATION      					   *
*******************************************************************/


Config read_config_file(const std::string& filename) {
    Config config;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Erreur: impossible d'ouvrir le fichier de configuration.");
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "DATA_PATHS") {
            while (std::getline(file, line) && line != "END_DATA_PATHS") {
                std::istringstream path_iss(line);
                std::string input_path, output_path;
                path_iss >> input_path >> output_path;
                config.DATA_PATHS.push_back({input_path, output_path});
            }
        } else if (key == "FILENAMES") {
            while (std::getline(file, line) && line != "END_FILENAMES") {
                std::istringstream filenames_iss(line);
                std::string filename;
                std::vector<std::string> filenames_in_directory;
                while (filenames_iss >> filename) {
                    filenames_in_directory.push_back(filename);
                }
                config.FILENAMES.push_back(filenames_in_directory);
            }
        } else if (key == "SUFFIXES") {
            iss >> config.POSITION_SUFFIX >> config.VELOCITY_SUFFIX;
        } else if (key == "VOLUME") {
            iss >> std::boolalpha >> config.BOOL_VOLUME;
        } else if (key == "DIVERGENCE") {
            iss >> std::boolalpha >> config.BOOL_DIVERGENCE;
        } else if (key == "CURL") {
            iss >> std::boolalpha >> config.BOOL_CURL;
        } else if (key == "GRAPH") {
            iss >> std::boolalpha >> config.BOOL_GRAPH;
        } else if (key == "WAVELET") {
            iss >> std::boolalpha >> config.BOOL_WAVELET;
        } else if (key == "SKIP_DELAUNAY_VOLUME_GRAPH") {
            iss >> std::boolalpha >> config.BOOL_SKIP_DELAUNAY_VOLUME_GRAPH;
        } else if (key == "VOLUME_MERGE") {
            iss >> std::boolalpha >> config.BOOL_VOLUME_MERGE;
        } else if (key == "DIVERGENCE_MERGE") {
            iss >> std::boolalpha >> config.BOOL_DIVERGENCE_MERGE;
        } else if (key == "CURL_MERGE") {
            iss >> std::boolalpha >> config.BOOL_CURL_MERGE;
        } else if (key == "GRAPH_MERGE") {
            iss >> std::boolalpha >> config.BOOL_GRAPH_MERGE;
        } else if (key == "DIVERGENCE_WRITE") {
            iss >> std::boolalpha >> config.BOOL_DIVERGENCE_WRITE;
        } else if (key == "CURL_WRITE") {
            iss >> std::boolalpha >> config.BOOL_CURL_WRITE;
        } else if (key == "HELICITY_WRITE") {
            iss >> std::boolalpha >> config.BOOL_HELICITY_WRITE;
        } else if (key == "SUBDOMAINS") {
            iss >> config.SUBDOMAINS;
        } else if (key == "LOAD_GRAPH") {
            iss >> std::boolalpha >> config.LOAD_GRAPH;
        } else if (key == "DELTA_T") {
            iss >> config.DELTA_T;
        } else if (key == "SLICE_SIZE") {
            iss >> config.SLICE_SIZE;
        } else if (key == "DOMAIN_SIZE") {
            double values[6];
            std::string value_str;
            for(int i = 0; i < 6; i++){
                if (!(iss >> value_str)) {
                    throw std::runtime_error("Not enough values for DOMAIN_SIZE");
                }
                if (value_str == "TWO_PI") values[i] = TWO_PI;
                else if (value_str == "PI") values[i] = PI;
                else values[i] = std::stod(value_str);
            }
            config.DOMAIN_SIZE = std::make_tuple(values[0], values[1], values[2], values[3], values[4], values[5]);
            if (iss >> value_str) {
                throw std::runtime_error("Too many values for DOMAIN_SIZE");
            }
        } else if (key == "SUBCUBES") {
            std::string value_str;
            while(iss >> value_str) {
                std::size_t dash_pos = value_str.find("-");
                if (dash_pos == std::string::npos) { // pas de tiret
                    int value = std::stoi(value_str);
                    config.SUBCUBES.push_back({value, value});
                } else { // il y a un tiret
                    int value1 = std::stoi(value_str.substr(0, dash_pos));
                    int value2 = std::stoi(value_str.substr(dash_pos + 1));
                    config.SUBCUBES.push_back({value1, value2});
                }
            }
        } else if (key == "BOOL_TEST_RANDOM") {
            iss >> std::boolalpha >> config.BOOL_TEST_RANDOM;
        } else if (key == "SEQUENTIAL") {
            iss >> std::boolalpha >> config.BOOL_SEQUENTIAL;
        }
    }

    return config;
}


/*******************************************************************
*					             WRITE      					   *
*******************************************************************/


void write_points_to_binary_file(const std::vector<Point>& points, const std::string& filename) {
    std::ofstream output_file(filename, std::ios::binary);

    std::cout << "WRITE DATA IN : " << filename << std::endl;

    if (!output_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (const auto& point : points) {
        double x = point.x();
        double y = point.y();
        double z = point.z();
        output_file.write(reinterpret_cast<const char*>(&x), sizeof(double));
        output_file.write(reinterpret_cast<const char*>(&y), sizeof(double));
        output_file.write(reinterpret_cast<const char*>(&z), sizeof(double));
    }

    output_file.close();
}


void write_points_with_indices_to_binary_file(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices, const std::string& filename) {
    std::ofstream output_file(filename, std::ios::binary);

    std::cout << "WRITE DATA IN : " << filename << std::endl;

    if (!output_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (const auto& point_with_index : points_with_indices) {
        const Point& point = point_with_index.first;
        double x = point.x();
        double y = point.y();
        double z = point.z();
        output_file.write(reinterpret_cast<const char*>(&x), sizeof(double));
        output_file.write(reinterpret_cast<const char*>(&y), sizeof(double));
        output_file.write(reinterpret_cast<const char*>(&z), sizeof(double));
    }

    output_file.close();
}


void write_double_array_to_binary_file(const std::vector<double>& double_array, const std::string& filename) {
    std::ofstream output_file(filename, std::ios::binary);

    std::cout << "WRITE DATA IN : " << filename << std::endl;

    if (!output_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (const double& value : double_array) {
        output_file.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }

    output_file.close();
}


void write_int_array_to_binary_file(const std::vector<std::size_t>& int_array, const std::string& filename) {
    std::ofstream output_file(filename, std::ios::binary);

    std::cout << "WRITE DATA IN : " << filename << std::endl;

    if (!output_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (const std::size_t& value : int_array) {
        output_file.write(reinterpret_cast<const char*>(&value), sizeof(std::size_t));
    }

    output_file.close();
}



void write_original_indicies_to_binary_file(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points, std::size_t Np, const std::string& filename) {
    std::ofstream output_file(filename, std::ios::binary);

    std::cout << "WRITE DATA IN : " << filename << std::endl;
    
    if (!output_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    for (std::size_t i = 0; i < Np; ++i) {
        const std::size_t index = points[i].second.second;
        output_file.write(reinterpret_cast<const char*>(&index), sizeof(std::size_t));
    }

    output_file.close();
}


/*******************************************************************
*					             READ       					   *
*******************************************************************/


void read_points_from_binary_file(const std::string& filename, std::vector<Point>& points) {
    std::ifstream input_file(filename, std::ios::binary);

    std::cout << "READ DATA FROM : " << filename << std::endl;

    if (!input_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return;
    }

    while (input_file) {
        double x, y, z;
        input_file.read(reinterpret_cast<char*>(&x), sizeof(double));
        input_file.read(reinterpret_cast<char*>(&y), sizeof(double));
        input_file.read(reinterpret_cast<char*>(&z), sizeof(double));

        if (input_file) {
            points.emplace_back(x, y, z);
        }
    }

    input_file.close();
}


bool read_points_with_indices_from_binary_file(const std::string& filename, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices) {
    std::ifstream input_file(filename, std::ios::binary);

    std::cout << "READ DATA WITH INDICES FROM : " << filename << std::endl;

    if (!input_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return false;
    }

    std::size_t index = 0;
    while (input_file) {
        double x, y, z;
        input_file.read(reinterpret_cast<char*>(&x), sizeof(double));
        input_file.read(reinterpret_cast<char*>(&y), sizeof(double));
        input_file.read(reinterpret_cast<char*>(&z), sizeof(double));

        if (input_file) {
            points_with_indices.emplace_back(Point(x, y, z), std::make_pair(index, index));
            ++index;
        }
    }

    input_file.close();

    return true;
}


bool read_double_array_from_binary_file(const std::string& filename, std::vector<double>& double_array) {
    std::ifstream input_file(filename, std::ios::binary);

    std::cout << "READ DATA FROM : " << filename << std::endl;

    if (!input_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return false;
    }

    while (input_file) {
        double value;
        input_file.read(reinterpret_cast<char*>(&value), sizeof(double));

        if (input_file) {
            double_array.push_back(value);
        }
    }

    input_file.close();

    return true;
}


bool read_int_from_binary_file(const std::string& filename, std::vector<std::size_t>& indices) {
    std::ifstream input_file(filename, std::ios::binary);

    std::cout << "READ DATA FROM : " << filename << std::endl;

    if (!input_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return false;
    }

    while (input_file) {
        std::size_t index;
        input_file.read(reinterpret_cast<char*>(&index), sizeof(std::size_t));

        if (input_file) {
            indices.push_back(index);
        }
    }

    input_file.close();

    return true;
}


/*******************************************************************
*					             MERGE       					   *
*******************************************************************/


void merge_data(const ArgConfig& ARG_CONFIG, int totalSubdomains, std::size_t max_index) {

    if (!ARG_CONFIG.BOOL_VOLUME_MERGE && !ARG_CONFIG.BOOL_DIVERGENCE_MERGE && !ARG_CONFIG.BOOL_CURL_MERGE && !ARG_CONFIG.BOOL_GRAPH_MERGE)
        return;

    std::cout << "MERGE DATA" << std::endl;

    // Determine the maximum index.
    if (max_index == -1){
        max_index = 0;
        for (int i = 0; i < totalSubdomains; ++i) {
            std::vector<std::size_t> index_subcube;
            read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".index_" + std::to_string(i), index_subcube);
            for (const auto& index : index_subcube) {
                if (index > max_index) {
                    max_index = index;
                }
            }
        }
    }
    std::cout << "Np = " << max_index + 1 << std::endl;

    // Create vectors of the correct size.
    std::vector<double> cell_volumes(max_index + 1, 0.0);
    std::vector<double> cell_volumes_x_y_z(max_index + 1, 0.0);
    std::vector<double> cell_volumes_x(max_index + 1, 0.0);
    std::vector<double> cell_volumes_y(max_index + 1, 0.0);
    std::vector<double> cell_volumes_z(max_index + 1, 0.0);

    std::vector<std::size_t> level(max_index + 1, 0.0);
    std::vector<std::size_t> merge(max_index + 1, 0.0);


    // Proceed as before.
    for (int i = 0; i < totalSubdomains; ++i) {
        std::vector<std::size_t> index_subcube;
        std::vector<double> cell_volumes_subcube, cell_volumes_x_y_z_subcube, cell_volumes_x_subcube, cell_volumes_y_subcube, cell_volumes_z_subcube;
        std::vector<std::size_t> key_subcube, level_subcube, merge_subcube;

        bool SUCCESS_VOL;
        bool SUCCESS_DIVERGENCE;
        bool SUCCESS_CURL;
        read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".index_" + std::to_string(i), index_subcube);
        if (ARG_CONFIG.BOOL_VOLUME_MERGE){
            SUCCESS_VOL = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".vol_" + std::to_string(i), cell_volumes_subcube);
        }
        if (ARG_CONFIG.BOOL_DIVERGENCE_MERGE){
            SUCCESS_DIVERGENCE = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".vol_x_y_z_" + std::to_string(i), cell_volumes_x_y_z_subcube);
        }
        if (ARG_CONFIG.BOOL_CURL_MERGE){
            SUCCESS_CURL = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".vol_0_z_-y_" + std::to_string(i), cell_volumes_x_subcube);
            read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".vol_-z_0_x_" + std::to_string(i), cell_volumes_y_subcube);
            read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + ".vol_y_-x_0_" + std::to_string(i), cell_volumes_z_subcube);
        }
        if (ARG_CONFIG.BOOL_GRAPH_MERGE){
            std::string subdomainNumber = "";
            if (totalSubdomains != 1) {
                subdomainNumber = "_" + std::to_string(i);
            }
            std::string subdirectory = ARG_CONFIG.DATA_PATH.second + "/GW_" + ARG_CONFIG.FILENAME;
            read_int_from_binary_file(subdirectory + "/" + ARG_CONFIG.FILENAME + subdomainNumber + ".key", key_subcube);
            read_int_from_binary_file(subdirectory + "/" + ARG_CONFIG.FILENAME + subdomainNumber + ".level", level_subcube);
            read_int_from_binary_file(subdirectory + "/" + ARG_CONFIG.FILENAME + subdomainNumber + ".merge", merge_subcube);
        }

        for (std::size_t j = 0; j < index_subcube.size(); ++j) {
            if (ARG_CONFIG.BOOL_VOLUME_MERGE && SUCCESS_VOL){
                cell_volumes[index_subcube[j]] = cell_volumes_subcube[j];
            }
            if (ARG_CONFIG.BOOL_DIVERGENCE_MERGE && SUCCESS_DIVERGENCE){
                cell_volumes_x_y_z[index_subcube[j]] = cell_volumes_x_y_z_subcube[j];
            }
            if (ARG_CONFIG.BOOL_CURL_MERGE && SUCCESS_CURL){
                cell_volumes_x[index_subcube[j]] = cell_volumes_x_subcube[j];
                cell_volumes_y[index_subcube[j]] = cell_volumes_y_subcube[j];
                cell_volumes_z[index_subcube[j]] = cell_volumes_z_subcube[j];
            }
            if (ARG_CONFIG.BOOL_GRAPH_MERGE){
                if (key_subcube[j] < max_index + 1){
                    level[key_subcube[j]] = level_subcube[j];
                    merge[key_subcube[j]] = merge_subcube[j];
                }
            }
        }
    }

    // Write the merged data to a file.
    if (ARG_CONFIG.BOOL_VOLUME_MERGE){
        write_double_array_to_binary_file(cell_volumes, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol");
    }
    if (ARG_CONFIG.BOOL_DIVERGENCE_MERGE){
        write_double_array_to_binary_file(cell_volumes_x_y_z, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_x_y_z");
    }
    if (ARG_CONFIG.BOOL_CURL_MERGE){
        write_double_array_to_binary_file(cell_volumes_x, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_0_z_-y");
        write_double_array_to_binary_file(cell_volumes_y, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_-z_0_x");
        write_double_array_to_binary_file(cell_volumes_z, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_y_-x_0");
    }
    if (ARG_CONFIG.BOOL_GRAPH_MERGE){
        write_int_array_to_binary_file(level, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".level");
        write_int_array_to_binary_file(merge, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".merge");
        // write_double_array_to_binary_file(gvol, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".gvol");
    }
}


void computes_quantities(const ArgConfig& ARG_CONFIG) {

    if (!ARG_CONFIG.BOOL_DIVERGENCE_WRITE && !ARG_CONFIG.BOOL_CURL_WRITE && !ARG_CONFIG.BOOL_HELICITY_WRITE)
        return;

    std::cout << "COMPUTES QUANTITIES" << std::endl;

    // Create vectors of the correct size.
    std::vector<double> cell_volumes;

    bool SUCCESS;
    SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol", cell_volumes);
    if (!SUCCESS) return;
    if (ARG_CONFIG.BOOL_DIVERGENCE_WRITE){
        std::vector<double> cell_volumes_x_y_z;
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_x_y_z", cell_volumes_x_y_z);
        if (!SUCCESS) return;
        std::vector<double> DIV = compute_discrete_velocity_divergence(cell_volumes, cell_volumes_x_y_z, ARG_CONFIG.DELTA_T);
        write_double_array_to_binary_file(DIV, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".div");
    }
    if (ARG_CONFIG.BOOL_CURL_WRITE){
        std::vector<double> cell_volumes_x;
        std::vector<double> cell_volumes_y;
        std::vector<double> cell_volumes_z;
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_0_z_-y", cell_volumes_x);
        if (!SUCCESS) return;
        std::vector<double> CURL_X = compute_discrete_velocity_divergence(cell_volumes, cell_volumes_x, ARG_CONFIG.DELTA_T);
        write_double_array_to_binary_file(CURL_X, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".curl_x");
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_-z_0_x", cell_volumes_y);
        if (!SUCCESS) return;
        std::vector<double> CURL_Y = compute_discrete_velocity_divergence(cell_volumes, cell_volumes_y, ARG_CONFIG.DELTA_T);
        write_double_array_to_binary_file(CURL_Y, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".curl_y");
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".vol_y_-x_0", cell_volumes_z);
        if (!SUCCESS) return;
        std::vector<double> CURL_Z = compute_discrete_velocity_divergence(cell_volumes, cell_volumes_z, ARG_CONFIG.DELTA_T);
        write_double_array_to_binary_file(CURL_Z, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".curl_z");
    }
    if (ARG_CONFIG.BOOL_HELICITY_WRITE){
        std::vector<double> CURL_X;
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".curl_x", CURL_X);
        if (!SUCCESS) return;
        std::vector<double> CURL_Y;
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".curl_y", CURL_Y);
        if (!SUCCESS) return;
        std::vector<double> CURL_Z;
        SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".curl_z", CURL_Z);
        if (!SUCCESS) return;
        std::vector<Point> VELOCITY;
        read_points_from_binary_file(ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME + ARG_CONFIG.VELOCITY_SUFFIX, VELOCITY);
        std::vector<double> HEL = compute_discrete_velocity_helicity(CURL_X, CURL_Y, CURL_Z, VELOCITY); 
        write_double_array_to_binary_file(HEL, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + ".hel");
    }
}


/*******************************************************************
*					            GENERATE       					   *
*******************************************************************/


void generate_points_in_unit_cube(int num_points, std::vector<Point>& points) {
    CGAL::Random_points_in_cube_3<Point> rnd(1.);

    std::cout << "RANDOM POINTS" << std::endl;

    points.reserve(num_points);
    for (int i = 0; i != num_points; ++i) {
        Point p = *rnd++;
        points.push_back(Point(p.x() / 2 + 0.5, p.y() / 2 + 0.5, p.z() / 2 + 0.5));
    }
}


void generate_points_with_indices_in_unit_cube(int num_points, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices) {
    CGAL::Random_points_in_cube_3<Point> rnd(1.);

    std::cout << "RANDOM POINTS WITH INDICES" << std::endl;

    points_with_indices.reserve(num_points);
    for (std::size_t i = 0; i != num_points; ++i) {
        Point p = *rnd++;
        points_with_indices.push_back(std::make_pair(Point(p.x() / 2 + 0.5, p.y() / 2 + 0.5, p.z() / 2 + 0.5), std::make_pair(i, i)));
    }
}


void generate_points_with_indices_in_two_pi_cube(int num_points, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices) {
    CGAL::Random_points_in_cube_3<Point> rnd(1.);

    std::cout << "RANDOM POINTS WITH INDICES IN TWO PI CUBE" << std::endl;

    points_with_indices.reserve(num_points);
    for (std::size_t i = 0; i != num_points; ++i) {
        Point p = *rnd++;
        points_with_indices.push_back(std::make_pair(Point(p.x() * TWO_PI / 2 + TWO_PI / 2, p.y() * TWO_PI / 2 + TWO_PI / 2, p.z() * TWO_PI / 2 + TWO_PI / 2), std::make_pair(i, i)));
    }
}


/*******************************************************************
*					           PERIODIZE       					   *
*******************************************************************/


std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> periodize_points(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices, const double SLICE_SIZE, const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE) {

    std::cout << "PERIODIZE POINT" << std::endl;

    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> periodized_points(P_with_indices);

    std::size_t next_index = P_with_indices.size();

    std::array<int, 3> range = { -1, 0, 1 };
    for (int i : range) {
        for (int j : range) {
            for (int k : range) {
                if (abs(i) + abs(j) + abs(k) != 0) {
                    for (const auto& point_index_pair : P_with_indices) {
                        Point shifted_point = Point(point_index_pair.first.x() + i * TWO_PI, point_index_pair.first.y() + j * TWO_PI, point_index_pair.first.z() + k * TWO_PI);

                        if (shifted_point.x() >= std::get<0>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.x() <= std::get<1>(DOMAIN_SIZE) + SLICE_SIZE &&
                            shifted_point.y() >= std::get<2>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.y() <= std::get<3>(DOMAIN_SIZE) + SLICE_SIZE &&
                            shifted_point.z() >= std::get<4>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.z() <= std::get<5>(DOMAIN_SIZE) + SLICE_SIZE) {
                            // periodized_points.push_back(std::make_pair(shifted_point, next_index));
                            //periodized_points.push_back(std::make_pair(shifted_point, point_index_pair.second));
                            periodized_points.push_back(std::make_pair(shifted_point, std::make_pair(next_index, point_index_pair.second.first)));
                            ++next_index;
                        }
                    }
                }
            }
        }
    }

    return periodized_points;
}


std::vector<Point> periodize_velocity(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices, const std::vector<Point>& V_with_indices, const double SLICE_SIZE, const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE) {

    std::cout << "PERIODIZE VELOCITY" << std::endl;
   
    std::vector<Point> periodized_velocity(V_with_indices);

    std::array<int, 3> range = { -1, 0, 1 };
    for (int i : range) {
        for (int j : range) {
            for (int k : range) {
                if (abs(i) + abs(j) + abs(k) != 0) {
                    for (const auto& point_index_pair : P_with_indices) {
                        Point shifted_point = Point(point_index_pair.first.x() + i * TWO_PI, point_index_pair.first.y() + j * TWO_PI, point_index_pair.first.z() + k * TWO_PI);

                        if (shifted_point.x() >= std::get<0>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.x() <= std::get<1>(DOMAIN_SIZE) + SLICE_SIZE &&
                            shifted_point.y() >= std::get<2>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.y() <= std::get<3>(DOMAIN_SIZE) + SLICE_SIZE &&
                            shifted_point.z() >= std::get<4>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.z() <= std::get<5>(DOMAIN_SIZE) + SLICE_SIZE) {
                            periodized_velocity.push_back(V_with_indices[point_index_pair.second.first]);
                        }
                    }
                }
            }
        }
    }

    return periodized_velocity;
}


std::tuple<std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>, std::vector<Point>, int> extract_subcube_points(
    const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized,
    const std::vector<Point>& V_with_indices_periodized,
    int subdomainNumber,
    int totalSubdomains,
    double SLICE_SIZE,
    const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE) {

    int cubeRoot = std::cbrt(totalSubdomains);
    if (cubeRoot * cubeRoot * cubeRoot != totalSubdomains) {
        throw std::invalid_argument("totalSubdomains must be a perfect cube.");
    }

    if (subdomainNumber < 0 || subdomainNumber >= totalSubdomains) {
        throw std::invalid_argument("subdomainNumber must be between 0 and totalSubdomains - 1.");
    }

    // Compute subdomain coordinates
    int xSub = subdomainNumber % cubeRoot;
    int ySub = (subdomainNumber / cubeRoot) % cubeRoot;
    int zSub = subdomainNumber / (cubeRoot * cubeRoot);

    double subdomainSize = TWO_PI / cubeRoot;

    // Compute subdomain limits
    double xMin = std::get<0>(DOMAIN_SIZE) + xSub * subdomainSize;
    double xMax = xMin + subdomainSize;
    double yMin = std::get<2>(DOMAIN_SIZE) + ySub * subdomainSize;
    double yMax = yMin + subdomainSize;
    double zMin = std::get<4>(DOMAIN_SIZE) + zSub * subdomainSize;
    double zMax = zMin + subdomainSize;

    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices_periodized_subcube;
    std::vector<Point> V_with_indices_periodized_subcube;
    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices_extended;
    std::vector<Point> V_with_indices_extended;

    int pointsInDomain = 0;

    for (std::size_t i = 0; i < P_with_indices_periodized.size(); ++i) {
        const Point& p = P_with_indices_periodized[i].first;
        bool inDomain = p.x() >= xMin && p.x() <= xMax && p.y() >= yMin && p.y() <= yMax && p.z() >= zMin && p.z() <= zMax;
        bool inExtendedDomain = p.x() >= xMin - SLICE_SIZE && p.x() <= xMax + SLICE_SIZE && p.y() >= yMin - SLICE_SIZE && p.y() <= yMax + SLICE_SIZE && p.z() >= zMin - SLICE_SIZE && p.z() <= zMax + SLICE_SIZE;

        if (inDomain) {
            P_with_indices_periodized_subcube.push_back(std::make_pair(P_with_indices_periodized[i].first, std::make_pair(pointsInDomain, P_with_indices_periodized[i].second.second)));
            V_with_indices_periodized_subcube.push_back(V_with_indices_periodized[i]);
            pointsInDomain++;
        } else if (inExtendedDomain) {
            P_with_indices_extended.push_back(std::make_pair(P_with_indices_periodized[i].first, std::make_pair(pointsInDomain, P_with_indices_periodized[i].second.second)));
            V_with_indices_extended.push_back(V_with_indices_periodized[i]);
        }
    }

    // Append extended points after the inDomain points
    P_with_indices_periodized_subcube.insert(P_with_indices_periodized_subcube.end(), P_with_indices_extended.begin(), P_with_indices_extended.end());
    V_with_indices_periodized_subcube.insert(V_with_indices_periodized_subcube.end(), V_with_indices_extended.begin(), V_with_indices_extended.end());

    // Fix indices for the points in the extended domain
    for (std::size_t i = pointsInDomain; i < P_with_indices_periodized_subcube.size(); ++i) {
        P_with_indices_periodized_subcube[i].second.first = i;
    }

    return std::make_tuple(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, pointsInDomain);
}


std::vector<size_t> extractGlobalIndices(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized_subcube, int Np_Domain) {
    std::vector<size_t> global_indices;
    for(int i = 0; i < Np_Domain; i++) {
        // element.second correspond au pair d'indices (local, global)
        // element.second.second correspond à l'indice global
        global_indices.push_back(P_with_indices_periodized_subcube[i].second.second);
    }
    return global_indices;
}


/*******************************************************************
*					       READ PERIODIZE   					   *
*******************************************************************/


std::tuple<std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>, std::vector<Point>, int, int> read_points_velocity_with_indices_from_binary_file_periodize_extract_subcube(
    const std::string& filename, 
    const std::string& position_suffix, 
    const std::string& velicity_suffix, 
    int subdomainNumber,
    int totalSubdomains,
    double SLICE_SIZE,
    const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE
) {

    std::cout << "READ DATA WITH INDICES FROM : " << filename << std::endl;

    std::string posFile = filename + position_suffix;
    std::string velFile = filename + velicity_suffix;
    
    std::ifstream pos_input_file(posFile, std::ios::binary);
    std::ifstream vel_input_file(velFile, std::ios::binary);

    if (!pos_input_file || !vel_input_file) {
        std::cerr << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
        return std::make_tuple(std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>(), std::vector<Point>(), 0, 0);
    }


    int cubeRoot = std::cbrt(totalSubdomains);
    if (cubeRoot * cubeRoot * cubeRoot != totalSubdomains) {
        throw std::invalid_argument("totalSubdomains must be a perfect cube.");
    }

    if (subdomainNumber < 0 || subdomainNumber >= totalSubdomains) {
        throw std::invalid_argument("subdomainNumber must be between 0 and totalSubdomains - 1.");
    }

    // Compute subdomain coordinates
    int xSub = subdomainNumber % cubeRoot;
    int ySub = (subdomainNumber / cubeRoot) % cubeRoot;
    int zSub = subdomainNumber / (cubeRoot * cubeRoot);

    double xsubdomainSize = get<1>(DOMAIN_SIZE) / cubeRoot;
    double ysubdomainSize = get<3>(DOMAIN_SIZE) / cubeRoot;
    double zsubdomainSize = get<5>(DOMAIN_SIZE) / cubeRoot;

    // Compute subdomain limits
    double xMin = std::get<0>(DOMAIN_SIZE) + xSub * xsubdomainSize;
    double xMax = xMin + xsubdomainSize;
    double yMin = std::get<2>(DOMAIN_SIZE) + ySub * ysubdomainSize;
    double yMax = yMin + ysubdomainSize;
    double zMin = std::get<4>(DOMAIN_SIZE) + zSub * zsubdomainSize;
    double zMax = zMin + zsubdomainSize;

    int pointsInDomain = 0;
    std::array<int, 3> range = { -1, 0, 1 };



    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices_periodized_subcube;
    std::vector<Point> V_with_indices_periodized_subcube;
    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices_extended;
    std::vector<Point> V_with_indices_extended;



    std::size_t index = 0;
    while (pos_input_file && vel_input_file) {
        double px, py, pz, vx, vy, vz;
        
        // Read position and velocity
        pos_input_file.read(reinterpret_cast<char*>(&px), sizeof(double));
        pos_input_file.read(reinterpret_cast<char*>(&py), sizeof(double));
        pos_input_file.read(reinterpret_cast<char*>(&pz), sizeof(double));
        
        vel_input_file.read(reinterpret_cast<char*>(&vx), sizeof(double));
        vel_input_file.read(reinterpret_cast<char*>(&vy), sizeof(double));
        vel_input_file.read(reinterpret_cast<char*>(&vz), sizeof(double));

        if (pos_input_file && vel_input_file) {
            bool inDomain;
            px = std::fmod(px, std::get<1>(DOMAIN_SIZE)); // !!!
            py = std::fmod(py, std::get<3>(DOMAIN_SIZE));
            pz = std::fmod(pz, std::get<5>(DOMAIN_SIZE));
            const Point& p = Point(px, py, pz);

            if (xSub == cubeRoot - 1) {                         // This is the last subdomain, include upper boundaries
                inDomain = p.x() >= xMin && p.x() <= xMax;
            } else {                                            // Other subdomains, exclude upper boundaries
                inDomain = p.x() >= xMin && p.x() < xMax;
            }

            if (ySub == cubeRoot - 1) {
                inDomain = inDomain && p.y() >= yMin && p.y() <= yMax;
            } else {
                inDomain = inDomain && p.y() >= yMin && p.y() < yMax;
            }

            if (zSub == cubeRoot - 1) {
                inDomain = inDomain && p.z() >= zMin && p.z() <= zMax;
            } else {
                inDomain = inDomain && p.z() >= zMin && p.z() < zMax;
            }
            // bool inDomain = p.x() >= xMin && p.x() <= xMax && p.y() >= yMin && p.y() <= yMax && p.z() >= zMin && p.z() <= zMax; // !!!
            // bool inDomain = p.x() >= xMin && p.x() < xMax && p.y() >= yMin && p.y() < yMax && p.z() >= zMin && p.z() < zMax; // !!!
            bool inExtendedDomain = p.x() >= xMin - SLICE_SIZE && p.x() <= xMax + SLICE_SIZE && p.y() >= yMin - SLICE_SIZE && p.y() <= yMax + SLICE_SIZE && p.z() >= zMin - SLICE_SIZE && p.z() <= zMax + SLICE_SIZE;

            if (inDomain) {
                P_with_indices_periodized_subcube.push_back(std::make_pair(Point(px, py, pz), std::make_pair(pointsInDomain, index)));
                V_with_indices_periodized_subcube.push_back(Point(vx, vy, vz));
                pointsInDomain++;
                // if (index == 41598607 || index == 831893555 || index == 308190216 || index == 435797158) { // !!!
                //     std::cout << "subdomainNumber: " << subdomainNumber << " (px, py, pz): (" << px << " , " << py << " , " << pz << << ") " std::endl;
                // }
            } else if (inExtendedDomain) {
                P_with_indices_extended.push_back(std::make_pair(Point(px, py, pz), std::make_pair(pointsInDomain, index)));
                V_with_indices_extended.push_back(Point(vx, vy, vz));
            }
            
            for (int i : range) {
                for (int j : range) {
                    for (int k : range) {
                        if (abs(i) + abs(j) + abs(k) != 0) {
                            Point shifted_point = Point(px + i * TWO_PI, py + j * TWO_PI, pz + k * TWO_PI);
                            bool inExtendedDomainPeriodic = shifted_point.x() >= xMin - SLICE_SIZE && shifted_point.x() <= xMax + SLICE_SIZE && shifted_point.y() >= yMin - SLICE_SIZE && shifted_point.y() <= yMax + SLICE_SIZE && shifted_point.z() >= zMin - SLICE_SIZE && shifted_point.z() <= zMax + SLICE_SIZE;

                            if (shifted_point.x() >= std::get<0>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.x() <= std::get<1>(DOMAIN_SIZE) + SLICE_SIZE &&
                                shifted_point.y() >= std::get<2>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.y() <= std::get<3>(DOMAIN_SIZE) + SLICE_SIZE &&
                                shifted_point.z() >= std::get<4>(DOMAIN_SIZE) - SLICE_SIZE && shifted_point.z() <= std::get<5>(DOMAIN_SIZE) + SLICE_SIZE &&
                                inExtendedDomainPeriodic) {
                                
                                P_with_indices_extended.push_back(std::make_pair(shifted_point, std::make_pair(pointsInDomain, index)));
                                V_with_indices_extended.push_back(Point(vx, vy, vz));
                            }
                        }
                    }
                }
            }
            ++index;
        }
    }

    pos_input_file.close();
    vel_input_file.close();
    

    // Append extended points after the inDomain points
    P_with_indices_periodized_subcube.insert(P_with_indices_periodized_subcube.end(), P_with_indices_extended.begin(), P_with_indices_extended.end());
    V_with_indices_periodized_subcube.insert(V_with_indices_periodized_subcube.end(), V_with_indices_extended.begin(), V_with_indices_extended.end());

    // Fix indices for the points in the extended domain
    for (std::size_t i = pointsInDomain; i < P_with_indices_periodized_subcube.size(); ++i) {
        P_with_indices_periodized_subcube[i].second.first = i;
    }


    return std::make_tuple(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, pointsInDomain, index); // !!!!!!!!!!!!!! a changer
}




/*******************************************************************
*					            ADVECT      					   *
*******************************************************************/


void advect_particles(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double Delta_t, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result) {
    if (T_A.size() != T_B.size()) {
        throw std::invalid_argument("Les deux tableaux de points doivent avoir la m�me taille");
    }

    result.clear();
    result.reserve(T_A.size());

    for (std::size_t i = 0; i < T_A.size(); ++i) {
        const Point& point_A = T_A[i].first;
        const Point& point_B = T_B[i];
        std::size_t index_local = T_A[i].second.first;
        std::size_t index_global = T_A[i].second.second;

        double x_new = point_A.x() + Delta_t * point_B.x();
        double y_new = point_A.y() + Delta_t * point_B.y();
        double z_new = point_A.z() + Delta_t * point_B.z();

        Point new_point(x_new, y_new, z_new);
        result.push_back(std::make_pair(new_point, std::make_pair(index_local, index_global)));
    }
}


void advect_particles_v_perp_x(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double Delta_t, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result) {
    if (T_A.size() != T_B.size()) {
        throw std::invalid_argument("Les deux tableaux de points doivent avoir la m�me taille");
    }

    result.clear();
    result.reserve(T_A.size());

    for (std::size_t i = 0; i < T_A.size(); ++i) {
        const Point& point_A = T_A[i].first;
        const Point& point_B = T_B[i];
        std::size_t index_local = T_A[i].second.first;
        std::size_t index_global = T_A[i].second.second;

        double x_new = point_A.x();
        double y_new = point_A.y() + Delta_t * point_B.z();
        double z_new = point_A.z() - Delta_t * point_B.y();

        Point new_point(x_new, y_new, z_new);
        result.push_back(std::make_pair(new_point, std::make_pair(index_local, index_global)));
    }
}


void advect_particles_v_perp_y(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double Delta_t, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result) {
    if (T_A.size() != T_B.size()) {
        throw std::invalid_argument("Les deux tableaux de points doivent avoir la m�me taille");
    }

    result.clear();
    result.reserve(T_A.size());

    for (std::size_t i = 0; i < T_A.size(); ++i) {
        const Point& point_A = T_A[i].first;
        const Point& point_B = T_B[i];
        std::size_t index_local = T_A[i].second.first;
        std::size_t index_global = T_A[i].second.second;

        double x_new = point_A.x() - Delta_t * point_B.z();
        double y_new = point_A.y();
        double z_new = point_A.z() + Delta_t * point_B.x();

        Point new_point(x_new, y_new, z_new);
        result.push_back(std::make_pair(new_point, std::make_pair(index_local, index_global)));
    }
}


void advect_particles_v_perp_z(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double Delta_t, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result) {
    if (T_A.size() != T_B.size()) {
        throw std::invalid_argument("Les deux tableaux de points doivent avoir la m�me taille");
    }

    result.clear();
    result.reserve(T_A.size());

    for (std::size_t i = 0; i < T_A.size(); ++i) {
        const Point& point_A = T_A[i].first;
        const Point& point_B = T_B[i];
        std::size_t index_local = T_A[i].second.first;
        std::size_t index_global = T_A[i].second.second;

        double x_new = point_A.x() + Delta_t * point_B.y();
        double y_new = point_A.y() - Delta_t * point_B.x();
        double z_new = point_A.z();

        Point new_point(x_new, y_new, z_new);
        result.push_back(std::make_pair(new_point, std::make_pair(index_local, index_global)));
    }
}


std::vector<double> compute_discrete_velocity_divergence(const std::vector<double>& cell_volumes, const std::vector<double>& cell_volumes_x_y_z, double Delta_t) {
    std::vector<double> result(cell_volumes.size());

    for (size_t i = 0; i < cell_volumes.size(); ++i) {
        result[i] = 2.0 / Delta_t * (cell_volumes_x_y_z[i] - cell_volumes[i]) / (cell_volumes_x_y_z[i] + cell_volumes[i]);
    }

    return result;
}


std::vector<double> compute_discrete_velocity_helicity(const std::vector<double>& CURL_X, const std::vector<double>& CURL_Y, const std::vector<double>& CURL_Z, const std::vector<Point>& VELOCITY) {
    std::vector<double> result(CURL_X.size());

    for (std::size_t i = 0; i < CURL_X.size(); ++i) {
        // Calcul du produit scalaire entre CURL et VELOCITY
        double dot_product = CURL_X[i] * VELOCITY[i].x() + CURL_Y[i] * VELOCITY[i].y() + CURL_Z[i] * VELOCITY[i].z();

        // Calcul de la norme de CURL
        double norm_curl = std::sqrt(CURL_X[i] * CURL_X[i] + CURL_Y[i] * CURL_Y[i] + CURL_Z[i] * CURL_Z[i]);

        // Calcul de la norme de VELOCITY
        double norm_velocity = std::sqrt(VELOCITY[i].x() * VELOCITY[i].x() + VELOCITY[i].y() * VELOCITY[i].y() + VELOCITY[i].z() * VELOCITY[i].z());

        // Calcul du cosinus de l'angle entre CURL et VELOCITY
        result[i] = dot_product / (norm_curl * norm_velocity);
    }

    // Cos_theta = []
    // for i in range(len(Curl_time)):
    //     Cos_theta.append(np.dot(Curl_time[i], Vel[i]))
    // Cos_theta = np.array(Cos_theta)/(np.linalg.norm(Curl_time, axis = 1) * np.linalg.norm(Vel, axis = 1))

    return result;
}


void compute_particle_velocity_test_function(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices, std::vector<Point>& velocity_with_indices) {
    velocity_with_indices.reserve(points_with_indices.size());

    for (const auto& point_with_index : points_with_indices) {
        const Point& point = point_with_index.first;
        double x = point.x();
        double y = point.y();
        double z = point.z();

        double u_x = std::sin(x) * std::cos(y) * std::cos(z);
        double u_y = std::cos(x) * std::sin(y) * std::cos(z);
        double u_z = std::cos(x) * std::cos(y) * std::sin(z);

        velocity_with_indices.push_back(Point(u_x, u_y, u_z));
    }
}


void compute_velocity_divergence_test_function(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices, std::vector<double>& divergence_values) {
    
    for (const auto& pair : P_with_indices) {
        const Point& p = pair.first;
        double x = p.x();
        double y = p.y();
        double z = p.z();

        double divergence = 3 * std::cos(x) * std::cos(y) * std::cos(z);
        
        divergence_values.push_back(divergence);
    }
    
}


















/*


int get_neighbor_rank(int xSub, int ySub, int zSub, int dx, int dy, int dz, int cubeRoot) {
    // Calculer les coordonnées du voisin
    int xNeighbor = (xSub + dx + cubeRoot) % cubeRoot;
    int yNeighbor = (ySub + dy + cubeRoot) % cubeRoot;
    int zNeighbor = (zSub + dz + cubeRoot) % cubeRoot;

    // Calculer le rank du voisin à partir de ses coordonnées
    int neighborRank = xNeighbor + yNeighbor * cubeRoot + zNeighbor * cubeRoot * cubeRoot;

    return neighborRank;
}




int cubeRoot = std::cbrt(totalSubdomains);
// Compute subdomain coordinates
int xSub = subdomainNumber % cubeRoot;
int ySub = (subdomainNumber / cubeRoot) % cubeRoot;
int zSub = subdomainNumber / (cubeRoot * cubeRoot);

// Get the rank of the neighbor to the right
int dx = 1; // move right
int dy = 0; // no movement in the y direction
int dz = 0; // no movement in the z direction
int rightNeighborRank = get_neighbor_rank(xSub, ySub, zSub, dx, dy, dz, cubeRoot);





*/



































