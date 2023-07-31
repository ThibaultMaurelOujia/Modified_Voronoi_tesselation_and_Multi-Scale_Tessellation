#ifndef UTILE_H
#define UTILE_H


#include <map>
#include <cmath>
#include <tuple>
#include <chrono>
#include <thread>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>

#include <unistd.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


const double PI = 3.141592653589793;
const double TWO_PI = 6.283185307179586;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;



struct Config {
    std::vector<std::pair<std::string, std::string>> DATA_PATHS; // pair of input and output paths
    std::vector<std::vector<std::string>> FILENAMES; // names of files in each directory
    std::string POSITION_SUFFIX;
    std::string VELOCITY_SUFFIX;
    bool BOOL_VOLUME;
    bool BOOL_DIVERGENCE;
    bool BOOL_CURL;
    bool BOOL_GRAPH;
    bool BOOL_WAVELET;
    bool BOOL_SKIP_DELAUNAY_VOLUME_GRAPH;
    bool BOOL_VOLUME_MERGE;
    bool BOOL_DIVERGENCE_MERGE;
    bool BOOL_CURL_MERGE;
    bool BOOL_GRAPH_MERGE;
    bool BOOL_DIVERGENCE_WRITE;
    bool BOOL_CURL_WRITE;
    bool BOOL_HELICITY_WRITE;
    double DELTA_T; 
    double SLICE_SIZE;
    int SUBDOMAINS;
    std::vector<std::pair<int, int>> SUBCUBES;
    bool LOAD_GRAPH;
    std::tuple<double, double, double, double, double, double> DOMAIN_SIZE;
    bool BOOL_TEST_RANDOM;
    bool BOOL_SEQUENTIAL;
};


struct ArgConfig {
    std::pair<std::string, std::string> DATA_PATH; // single pair of input and output paths
    std::string FILENAME; // name of files in the directory
    std::string POSITION_SUFFIX;
    std::string VELOCITY_SUFFIX;
    bool BOOL_VOLUME;
    bool BOOL_DIVERGENCE;
    bool BOOL_CURL;
    bool BOOL_GRAPH;
    bool BOOL_WAVELET;
    bool BOOL_SKIP_DELAUNAY_VOLUME_GRAPH;
    bool BOOL_VOLUME_MERGE;
    bool BOOL_DIVERGENCE_MERGE;
    bool BOOL_CURL_MERGE;
    bool BOOL_GRAPH_MERGE;
    bool BOOL_DIVERGENCE_WRITE;
    bool BOOL_CURL_WRITE;
    bool BOOL_HELICITY_WRITE;
    double DELTA_T; 
    double SLICE_SIZE;
    int SUBDOMAINS;
    std::vector<std::pair<int, int>> SUBCUBES;
    bool LOAD_GRAPH;
    std::tuple<double, double, double, double, double, double> DOMAIN_SIZE;
    bool BOOL_TEST_RANDOM;
    bool BOOL_SEQUENTIAL;
};


Config read_config_file(const std::string& filename);


void write_points_to_binary_file(const std::vector<Point>& points, const std::string& filename);
void write_points_with_indices_to_binary_file(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices, const std::string& filename);
void write_double_array_to_binary_file(const std::vector<double>& double_array, const std::string& filename);
void write_int_array_to_binary_file(const std::vector<std::size_t>& int_array, const std::string& filename);
void write_original_indicies_to_binary_file(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points, std::size_t Np, const std::string& filename);

void read_points_from_binary_file(const std::string& filename, std::vector<Point>& points);
bool read_points_with_indices_from_binary_file(const std::string& filename, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices);
bool read_double_array_from_binary_file(const std::string& filename, std::vector<double>& double_array);
bool read_int_from_binary_file(const std::string& filename, std::vector<std::size_t>& indices);


void merge_data(const ArgConfig& ARG_CONFIG, int totalSubdomains, std::size_t max_index = -1);
void computes_quantities(const ArgConfig& ARG_CONFIG);


void generate_points_in_unit_cube(int num_points, std::vector<Point>& points);
void generate_points_with_indices_in_unit_cube(int num_points, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices);
void generate_points_with_indices_in_two_pi_cube(int num_points, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices);


std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> periodize_points(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices, const double SLICE_SIZE, const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE);
std::vector<Point> periodize_velocity(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices, const std::vector<Point>& V_with_indices, const double SLICE_SIZE, const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE);

std::tuple<std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>, std::vector<Point>, int> extract_subcube_points(
    const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized,
    const std::vector<Point>& V_with_indices_periodized,
    int subdomainNumber,
    int totalSubdomains,
    double SLICE_SIZE,
    const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE);

std::vector<size_t> extractGlobalIndices(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized_subcube, int Np_Domain);


std::tuple<std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>, std::vector<Point>, int, int> read_points_velocity_with_indices_from_binary_file_periodize_extract_subcube(
    const std::string& filename, 
    const std::string& position_suffix, 
    const std::string& velicity_suffix, 
    int subdomainNumber,
    int totalSubdomains,
    double SLICE_SIZE,
    const std::tuple<double, double, double, double, double, double>& DOMAIN_SIZE);


void advect_particles(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double DT, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result);
void advect_particles_v_perp_x(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double DT, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result);
void advect_particles_v_perp_y(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double DT, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result);
void advect_particles_v_perp_z(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& T_A, const std::vector<Point>& T_B, double DT, std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& result);

std::vector<double> compute_discrete_velocity_divergence(const std::vector<double>& cell_volumes, const std::vector<double>& cell_volumes_x_y_z, double DT);
std::vector<double> compute_discrete_velocity_helicity(const std::vector<double>& CURL_X, const std::vector<double>& CURL_Y, const std::vector<double>& CURL_Z, const std::vector<Point>& VELOCITY);



void compute_particle_velocity_test_function(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& points_with_indices, std::vector<Point>& velocity_with_indices);
void compute_velocity_divergence_test_function(const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices, std::vector<double>& divergence_values);




#endif // UTILE_H








