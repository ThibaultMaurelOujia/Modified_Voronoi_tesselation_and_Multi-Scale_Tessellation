#include "parallel_delaunay.h"



bool parallel_delaunay_multiscale(const ArgConfig& ARG_CONFIG, int subdomainNumber, int totalSubdomains) {

    std::cout << "\n\n\n\n";
    auto start_time_total = std::chrono::high_resolution_clock::now();

    
    //omp_set_max_active_levels(1); // Activating interlocking parallelism
    int num_threads;
    #pragma omp parallel
    {
    #pragma omp single
        {
            num_threads = omp_get_num_threads();
            std::cout << "ACTIVE THREADS: " << num_threads << std::endl;
        }
    }
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, num_threads);

    #ifdef CGAL_LINKED_WITH_TBB
        std::cout << "CGAL_LINKED_WITH_TBB is defined." << std::endl;
    #else
        std::cout << "CGAL_LINKED_WITH_TBB is NOT defined." << std::endl;
    #endif // CGAL_LINKED_WITH_TBB


    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);


    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    // Delaunay T3
    typedef CGAL::Triangulation_data_structure_3<
        // CGAL::Triangulation_vertex_base_3<K>,
        // CGAL::Triangulation_vertex_base_with_info_3<std::size_t, K>,
        CGAL::Triangulation_vertex_base_with_info_3<std::pair<std::size_t, std::size_t>, K>,
        CGAL::Delaunay_triangulation_cell_base_3<K>,
        CGAL::Parallel_tag>                               Tds;
    typedef CGAL::Delaunay_triangulation_3<K, Tds>              Triangulation;
    typedef Triangulation::Point                                Point;


    std::cout << "\n\n";


    /*******************************************************************
    *					          Load data      					   *
    *******************************************************************/


    std::cout << "START LOAD DATA" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices_periodized_subcube;
    std::vector<Point> V_with_indices_periodized_subcube;
    int Np = 0;
    int Np_Domain = 0;
    int FILENAME_Np;
    if (ARG_CONFIG.BOOL_TEST_RANDOM) {
        try {
            FILENAME_Np = std::stod(ARG_CONFIG.FILENAME);
        } catch (std::invalid_argument const &e) {
            std::cout << "Bad input: std::invalid_argument thrown" << '\n';
            return false;
        } catch (std::out_of_range const &e) {
            std::cout << "Integer overflow: std::out_of_range thrown" << '\n';
            return false;
        }
        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices;
        std::vector<Point> V_with_indices;
        generate_points_with_indices_in_two_pi_cube(FILENAME_Np, P_with_indices);
        compute_particle_velocity_test_function(P_with_indices, V_with_indices);

        // Write positions and velocities to files
        std::filesystem::create_directories(ARG_CONFIG.DATA_PATH.first);
        std::string pos_filename = ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME + ARG_CONFIG.POSITION_SUFFIX;
        std::string vel_filename = ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME + ARG_CONFIG.VELOCITY_SUFFIX;
        write_points_with_indices_to_binary_file(P_with_indices, pos_filename);
        write_points_to_binary_file(V_with_indices, vel_filename);

        Np = static_cast<int>(P_with_indices.size());

        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_with_indices_periodized = periodize_points(P_with_indices, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE);
        std::vector<Point> V_with_indices_periodized = periodize_velocity(P_with_indices, V_with_indices, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE);
        std::cout << "Np = " << static_cast<int>(P_with_indices.size()) << " | Np periodized = " << static_cast<int>(P_with_indices_periodized.size()) << std::endl;

        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_with_indices);
        std::vector<Point>().swap(V_with_indices);

        std::tie(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, Np_Domain) = extract_subcube_points(P_with_indices_periodized, V_with_indices_periodized, subdomainNumber, totalSubdomains, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE);
        int Np_Domain_extended = static_cast<int>(P_with_indices_periodized_subcube.size());
        std::cout << "Np sub-cube = " << Np_Domain << " | Np sub-cube slice = " << Np_Domain_extended << std::endl;

        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_with_indices_periodized);
        std::vector<Point>().swap(V_with_indices_periodized);
    } else {
        std::tie(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, Np_Domain, Np) = read_points_velocity_with_indices_from_binary_file_periodize_extract_subcube(ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME, ARG_CONFIG.POSITION_SUFFIX, ARG_CONFIG.VELOCITY_SUFFIX, subdomainNumber, totalSubdomains, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE);
        std::cout << "Np = " << Np << " | Np sub-cube = " << Np_Domain << " | Np sub-cube slice = " << static_cast<int>(P_with_indices_periodized_subcube.size()) << std::endl;
        if (Np == 0) return false;
    }
    

    // write_points_with_indices_to_binary_file(P_with_indices_periodized_subcube, "../data/subcube/P_with_indices_periodized_subcube_" + std::to_string(subdomainNumber) + ".pos");
    

    std::vector<int16_t> GlobalIndices_GlobalProcess(Np, -1); // Set to -1 (meaning no associated rank)
    if (ARG_CONFIG.BOOL_GRAPH){
        updateGlobalIndicesGlobalProcess(GlobalIndices_GlobalProcess, P_with_indices_periodized_subcube, Np, Np_Domain, subdomainNumber, MPI_rank, MPI_size, ARG_CONFIG);
    }


    std::cout << "END LOAD DATA" << std::endl;
    std::cout << "_____________________________" << std::endl;


    /*******************************************************************
    *					   Delaunay tessellation					   *
    *******************************************************************/


    std::cout << "START DELAUNAY" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);


    updatePeriodizedSubcube(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, ARG_CONFIG, Np_Domain);


    std::cout << "Np = " << Np << " | Np sub-cube = " << Np_Domain << " | Np sub-cube slice = " << static_cast<int>(P_with_indices_periodized_subcube.size()) << std::endl;


    auto start_time_Delaunay_2 = std::chrono::high_resolution_clock::now();

    // Construct the locking data-structure, using the bounding-box of the points
    Triangulation::Lock_data_structure locking_ds( \
    CGAL::Bbox_3(std::get<0>(ARG_CONFIG.DOMAIN_SIZE) - ARG_CONFIG.SLICE_SIZE, \
                std::get<2>(ARG_CONFIG.DOMAIN_SIZE) - ARG_CONFIG.SLICE_SIZE, \
                std::get<4>(ARG_CONFIG.DOMAIN_SIZE) - ARG_CONFIG.SLICE_SIZE, \
                std::get<1>(ARG_CONFIG.DOMAIN_SIZE) + ARG_CONFIG.SLICE_SIZE, \
                std::get<3>(ARG_CONFIG.DOMAIN_SIZE) + ARG_CONFIG.SLICE_SIZE, \
                std::get<5>(ARG_CONFIG.DOMAIN_SIZE) + ARG_CONFIG.SLICE_SIZE), \
    64);
    
    // Construct the triangulation in parallel
    Triangulation T(P_with_indices_periodized_subcube.begin(), P_with_indices_periodized_subcube.end(), &locking_ds);
    assert(T.is_valid());

    auto end_time_Delaunay_2 = std::chrono::high_resolution_clock::now();
    auto duration_Delaunay_2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_Delaunay_2 - start_time_Delaunay_2).count();
    std::cout << "DELAUNAY TIME: " << duration_Delaunay_2 / 1000. << " s" << std::endl;


    std::cout << "END DELAUNAY" << std::endl;
    std::cout << "_____________________________" << std::endl;


    /*******************************************************************
    *				    Modified Voronoi tessellation 				   *
    *******************************************************************/


    std::cout << "START MODIFIED VORONOI" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);


    std::vector<double> cell_volumes;


    {
        ModifiedVoronoi tessellation;

        bool exact_volume = true; // true false

        std::vector<double> cell_volumes_x_y_z;
        std::vector<double> cell_volumes_x;
        std::vector<double> cell_volumes_y;
        std::vector<double> cell_volumes_z;

        tessellation.create_empty_cells(P_with_indices_periodized_subcube.size());
        tessellation.create_empty_vertices(T);


        // Compute the volume
        if (ARG_CONFIG.BOOL_VOLUME){
            // MPI_Barrier(MPI_COMM_WORLD);
            tessellation.store_gravity_centers_and_indices(T, P_with_indices_periodized_subcube);
            if (exact_volume){
                tessellation.create_facet_modified_voronoi_cells(true, Np_Domain);
            }
            tessellation.compute_volume_modified_voronoi_cells(-1, true, exact_volume, Np_Domain);
            cell_volumes = tessellation.get_cells_volumes(Np_Domain);
            write_original_indicies_to_binary_file(P_with_indices_periodized_subcube, Np_Domain, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "index_" + std::to_string(subdomainNumber));
            write_double_array_to_binary_file(cell_volumes, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol_" + std::to_string(subdomainNumber));
        }
        
        
        if (ARG_CONFIG.BOOL_DIVERGENCE){
            // MPI_Barrier(MPI_COMM_WORLD);
            std::cout << "Delta t = " << ARG_CONFIG.DELTA_T << std::endl;
            
            // Advect particles using P, V and time step Delta_t
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_x_y_z_with_indices_periodized_subcube;
            advect_particles(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, ARG_CONFIG.DELTA_T, P_x_y_z_with_indices_periodized_subcube);
            
            // Compute the new volume
            tessellation.store_gravity_centers_and_indices(T, P_x_y_z_with_indices_periodized_subcube);
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_x_y_z_with_indices_periodized_subcube);
            tessellation.compute_volume_modified_voronoi_cells(-1, true, exact_volume, Np_Domain);
            cell_volumes_x_y_z = tessellation.get_cells_volumes(Np_Domain);
            write_double_array_to_binary_file(cell_volumes_x_y_z, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol_x_y_z_" + std::to_string(subdomainNumber));
            std::vector<double>().swap(cell_volumes_x_y_z);
        }
        
        
        if (ARG_CONFIG.BOOL_CURL){
            // MPI_Barrier(MPI_COMM_WORLD);
            std::cout << "Delta t = " << ARG_CONFIG.DELTA_T << std::endl;


            // Advect particles using orthogonal velocity
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_x_with_indices_periodized_subcube;
            advect_particles_v_perp_x(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, ARG_CONFIG.DELTA_T, P_x_with_indices_periodized_subcube);

            // Compute the new volume
            tessellation.store_gravity_centers_and_indices(T, P_x_with_indices_periodized_subcube);
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_x_with_indices_periodized_subcube);
            tessellation.compute_volume_modified_voronoi_cells(-1, true, exact_volume, Np_Domain);
            cell_volumes_x = tessellation.get_cells_volumes(Np_Domain);
            write_double_array_to_binary_file(cell_volumes_x, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol_0_z_-y_" + std::to_string(subdomainNumber));
            std::vector<double>().swap(cell_volumes_x);


            // Advect particles using orthogonal velocity
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_y_with_indices_periodized_subcube;
            advect_particles_v_perp_y(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, ARG_CONFIG.DELTA_T, P_y_with_indices_periodized_subcube);
            
            // Compute the new volume
            tessellation.store_gravity_centers_and_indices(T, P_y_with_indices_periodized_subcube);
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_y_with_indices_periodized_subcube);
            tessellation.compute_volume_modified_voronoi_cells(-1, true, exact_volume, Np_Domain);
            cell_volumes_y = tessellation.get_cells_volumes(Np_Domain);
            write_double_array_to_binary_file(cell_volumes_y, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol_-z_0_x_" + std::to_string(subdomainNumber));
            std::vector<double>().swap(cell_volumes_y);


            // Advect particles using orthogonal velocity
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_z_with_indices_periodized_subcube;
            advect_particles_v_perp_z(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, ARG_CONFIG.DELTA_T, P_z_with_indices_periodized_subcube);
            
            // Compute the new volume
            tessellation.store_gravity_centers_and_indices(T, P_z_with_indices_periodized_subcube);
            std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_z_with_indices_periodized_subcube);
            tessellation.compute_volume_modified_voronoi_cells(-1, true, exact_volume, Np_Domain);
            cell_volumes_z = tessellation.get_cells_volumes(Np_Domain);
            write_double_array_to_binary_file(cell_volumes_z, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol_y_-x_0_" + std::to_string(subdomainNumber));
            std::vector<double>().swap(cell_volumes_z);
        }

        
        bool has_velocity_gradient_components = false;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (ARG_CONFIG.velocity_gradient_tensor_components[i][j]) {
                    has_velocity_gradient_components = true;
                    break;
                }
            }
            if (has_velocity_gradient_components) break;
        }

        if (has_velocity_gradient_components) {
            std::cout << "Calcul des composantes du tenseur du gradient de vitesse" << std::endl;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    if (ARG_CONFIG.velocity_gradient_tensor_components[i][j]) {
                        std::cout << "Calcul de du" << component_to_string(i) << "/d" << component_to_string(j) << std::endl;

                        // Advect particles accordingly
                        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> P_advected_with_indices;
                        advect_particles_component(
                            P_with_indices_periodized_subcube,
                            V_with_indices_periodized_subcube,
                            ARG_CONFIG.DELTA_T,
                            i,
                            j,
                            P_advected_with_indices
                        );

                        // Compute the new volume
                        tessellation.store_gravity_centers_and_indices(T, P_advected_with_indices);
                        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_advected_with_indices);
                        tessellation.compute_volume_modified_voronoi_cells(-1, true, exact_volume, Np_Domain);
                        std::vector<double> cell_volumes_component = tessellation.get_cells_volumes(Np_Domain);

                        // Write the volumes to a file
                        std::string filename = ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." +
                            ARG_CONFIG.OUT_SUFFIXES + "vol_du" + component_to_string(i) + "_d" + component_to_string(j) +
                            "_" + std::to_string(subdomainNumber);
                        write_double_array_to_binary_file(cell_volumes_component, filename);
                        std::vector<double>().swap(cell_volumes_component);
                    }
                }
            }
        }
        
        

        // if (ARG_CONFIG.BOOL_TEST_RANDOM){
        //     // MPI_Barrier(MPI_COMM_WORLD);
        //     std::vector<double> Divergence_P_discrete = compute_discrete_velocity_divergence(cell_volumes, cell_volumes_x_y_z, ARG_CONFIG.DELTA_T);
        //     write_double_array_to_binary_file(Divergence_P_discrete, "../data/divergence_discrete.div");

        //     std::vector<double> Divergence_P_exact;
        //     compute_velocity_divergence_test_function(P_with_indices_periodized_subcube, Divergence_P_exact);
        //     write_double_array_to_binary_file(Divergence_P_exact, "../data/divergence_exact.div");

        //     write_points_with_indices_to_binary_file(P_with_indices_periodized_subcube, "../data/random_points.pos");
        // }
        

        std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>().swap(P_with_indices_periodized_subcube);
        std::vector<Point>().swap(V_with_indices_periodized_subcube);
    }


    std::cout << "END MODIFIED VORONOI" << std::endl;
    std::cout << "_____________________________" << std::endl;


    /*******************************************************************
    *				            GRAPH-WAVELETS 		        		   *
    *******************************************************************/


    if (ARG_CONFIG.BOOL_GRAPH){
        std::cout << "START GRAPH-WAVELETS" << std::endl;


        int8_t num_levels;
        if (Np != 0){
            num_levels = static_cast<int>(std::floor(std::log2(Np))); 
        } else {
            num_levels = 0;
        }
        if (ARG_CONFIG.SUBDOMAINS != -1 && ARG_CONFIG.SUBDOMAINS != 0 && ARG_CONFIG.SUBDOMAINS != 1){
            num_levels = 0;
        }
        // num_levels = 3; // !!!!!!!!!!!!!!
        std::cout << "NUMBER OF LEVELS: " << static_cast<int>(num_levels) << std::endl;


        Graph_Wavelets graphWavelets(ARG_CONFIG, T, cell_volumes, GlobalIndices_GlobalProcess, Np, Np_Domain, num_levels);


        MPI_Barrier(MPI_COMM_WORLD);
        
        if (MPI_size == 1 && (ARG_CONFIG.SUBDOMAINS == -1 || ARG_CONFIG.SUBDOMAINS == 0 || ARG_CONFIG.SUBDOMAINS == 1)) {
            graphWavelets.writeToBinary(ARG_CONFIG.DATA_PATH.second, ARG_CONFIG.FILENAME, ARG_CONFIG.OUT_SUFFIXES);
        } else {
            std::string filename = ARG_CONFIG.FILENAME + (num_levels == 0 ? "_LVL0" : "");
            graphWavelets.writeToBinary(ARG_CONFIG.DATA_PATH.second, filename, ARG_CONFIG.OUT_SUFFIXES, subdomainNumber);
            // graphWavelets.writeToBinary(ARG_CONFIG.DATA_PATH.second, ARG_CONFIG.FILENAME, subdomainNumber);
        }


        std::cout << "END GRAPH-WAVELETS" << std::endl;
        std::cout << "_____________________________" << std::endl;
    }


    /*******************************************************************
    *				                  END   		        		   *
    *******************************************************************/


    auto end_time_total = std::chrono::high_resolution_clock::now();
    auto duration_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_total - start_time_total).count();
    std::cout << "TOTAL TIME: " << duration_total / 1000. << " s" << std::endl << std::endl;

    return true;
}










bool parallel_delaunay_multiscale_2D(const ArgConfig& ARG_CONFIG, int subdomainNumber, int totalSubdomains) {

    std::cout << "\n\n\n\n";
    auto start_time_total = std::chrono::high_resolution_clock::now();


    // omp_set_max_active_levels(1); // Activating interlocking parallelism
    int num_threads;
    #pragma omp parallel
    {
    #pragma omp single
        {
            num_threads = omp_get_num_threads();
            std::cout << "ACTIVE THREADS: " << num_threads << std::endl;
        }
    }
    tbb::global_control c(tbb::global_control::max_allowed_parallelism, num_threads);

    #ifdef CGAL_LINKED_WITH_TBB
        std::cout << "CGAL_LINKED_WITH_TBB is defined." << std::endl;
    #else
        std::cout << "CGAL_LINKED_WITH_TBB is NOT defined." << std::endl;
    #endif // CGAL_LINKED_WITH_TBB


    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

    // Define types for 2D
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    // Delaunay T2
    typedef CGAL::Triangulation_data_structure_2<
        // CGAL::Triangulation_vertex_base_2<K>,
        // CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K>,
        CGAL::Triangulation_vertex_base_with_info_2<std::pair<std::size_t, std::size_t>, K>,
        CGAL::Triangulation_face_base_2<K>
    >                               Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds>              Triangulation;
    // typedef Triangulation::Point                                Point_2D;


    std::cout << "\n\n";


    /*******************************************************************
    *                            Load data                             *
    *******************************************************************/


    std::cout << "START LOAD DATA" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<std::pair<Point_2D, std::pair<std::size_t, std::size_t>>> P_with_indices_periodized_subcube;
    std::vector<Point_2D> V_with_indices_periodized_subcube;
    int Np = 0;
    int Np_Domain = 0;
    int FILENAME_Np;
    if (ARG_CONFIG.BOOL_TEST_RANDOM) {
        try {
            FILENAME_Np = std::stod(ARG_CONFIG.FILENAME);
        } catch (std::invalid_argument const &e) {
            std::cout << "Bad input: std::invalid_argument thrown" << '\n';
            return false;
        } catch (std::out_of_range const &e) {
            std::cout << "Integer overflow: std::out_of_range thrown" << '\n';
            return false;
        }
        std::vector<std::pair<Point_2D, std::pair<std::size_t, std::size_t>>> P_with_indices;
        std::vector<Point_2D> V_with_indices;
        generate_points_with_indices_in_two_pi_square(FILENAME_Np, P_with_indices); // !!!
        compute_particle_velocity_test_function_2D(P_with_indices, V_with_indices); // !!!

        // Write positions and velocities to files
        std::filesystem::create_directories(ARG_CONFIG.DATA_PATH.first);
        std::string pos_filename = ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME + ARG_CONFIG.POSITION_SUFFIX;
        std::string vel_filename = ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME + ARG_CONFIG.VELOCITY_SUFFIX;
        write_points_with_indices_to_binary_file_2D(P_with_indices, pos_filename);
        write_points_to_binary_file_2D(V_with_indices, vel_filename);

        Np = static_cast<int>(P_with_indices.size());

        std::vector<std::pair<Point_2D, std::pair<std::size_t, std::size_t>>> P_with_indices_periodized = periodize_points_2D(P_with_indices, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE); // !!!
        std::vector<Point_2D> V_with_indices_periodized = periodize_velocity_2D(P_with_indices, V_with_indices, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE); // !!!
        std::cout << "Np = " << static_cast<int>(P_with_indices.size()) << " | Np periodized = " << static_cast<int>(P_with_indices_periodized.size()) << std::endl;

        std::vector<std::pair<Point_2D, std::pair<std::size_t, std::size_t>>>().swap(P_with_indices);
        std::vector<Point_2D>().swap(V_with_indices);

        std::tie(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, Np_Domain) = extract_subsquare_points(P_with_indices_periodized, V_with_indices_periodized, subdomainNumber, totalSubdomains, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE); // !!!
        int Np_Domain_extended = static_cast<int>(P_with_indices_periodized_subcube.size());
        std::cout << "Np sub-cube = " << Np_Domain << " | Np sub-cube slice = " << Np_Domain_extended << std::endl;

        std::vector<std::pair<Point_2D, std::pair<std::size_t, std::size_t>>>().swap(P_with_indices_periodized);
        std::vector<Point_2D>().swap(V_with_indices_periodized);
    } else {
        std::tie(P_with_indices_periodized_subcube, V_with_indices_periodized_subcube, Np_Domain, Np) = read_points_velocity_with_indices_from_binary_file_periodize_extract_subsquare(ARG_CONFIG.DATA_PATH.first + "/" + ARG_CONFIG.FILENAME, ARG_CONFIG.POSITION_SUFFIX, ARG_CONFIG.VELOCITY_SUFFIX, subdomainNumber, totalSubdomains, ARG_CONFIG.SLICE_SIZE, ARG_CONFIG.DOMAIN_SIZE); // !!!
        std::cout << "Np = " << Np << " | Np sub-cube = " << Np_Domain << " | Np sub-cube slice = " << static_cast<int>(P_with_indices_periodized_subcube.size()) << std::endl;
        if (Np == 0) return false;
    }

    // Rest of the code with similar modifications...

    return true;
}






















/*******************************************************************
*				              FUNCTIONS   		        		   *
*******************************************************************/


bool updateGlobalIndicesGlobalProcess(std::vector<int16_t>& GlobalIndices_GlobalProcess,
                                      const std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized_subcube,
                                      std::size_t Np, std::size_t Np_Domain, int subdomainNumber, int MPI_rank, int MPI_size, const ArgConfig& ARG_CONFIG) 
{

    std::vector<size_t> GlobalIndices_LocalProcess = extractGlobalIndices(P_with_indices_periodized_subcube, Np_Domain);

    // We use MPI_Allgather and MPI_Allgatherv to share GlobalIndices_LocalProcess
    std::vector<int> sendcounts(MPI_size, 0), displs(MPI_size, 0);

    // Determines how much data each process should send
    sendcounts[MPI_rank] = GlobalIndices_LocalProcess.size();

    // Gathers the send accounts of all processes
    int MPI_err_size = MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, sendcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    if(MPI_err_size == MPI_SUCCESS) {
        std::cout << "SEND SIZE " << MPI_rank << " : SUCCESS" << std::endl;
    } else {
        std::cout << "SEND SIZE " << MPI_rank << " : FAIL" << std::endl;
    }
    

    // Calculates displacements for MPI_Allgatherv
    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + sendcounts[i-1];
    }

    // Gathers all data in all_data
    std::vector<size_t> all_data(Np);

    {
        // Security for reading data problem
        if(Np != displs[MPI_size-1] + sendcounts[MPI_size-1] && (ARG_CONFIG.SUBDOMAINS == -1 || ARG_CONFIG.SUBDOMAINS == 0 || ARG_CONFIG.SUBDOMAINS == 1)) {
            std::cout << "ERROR POINT:" << Np << "!=" << displs[MPI_size-1] + sendcounts[MPI_size-1] << std::endl;
            write_original_indicies_to_binary_file(P_with_indices_periodized_subcube, Np_Domain, ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "index_error_" + std::to_string(subdomainNumber));
            MPI_Barrier(MPI_COMM_WORLD);

            if (MPI_rank == 0) {
                all_data.clear(); 

                for (int rank = 0; rank < MPI_size; ++rank) {
                    std::vector<size_t> indices_from_file;
                    std::string file_to_read = ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "index_error_" + std::to_string(rank);
                    
                    if(read_int_from_binary_file(file_to_read, indices_from_file)) {
                        all_data.insert(all_data.end(), indices_from_file.begin(), indices_from_file.end());
                    } else {
                        std::cerr << "Erreur lors de la lecture du fichier " << file_to_read << std::endl;
                    }
                }

                // Trouver les valeurs en double
                std::unordered_map<size_t, size_t> counts;
                for (const auto &value : all_data) {
                    counts[value]++;
                }

                // Créer un vecteur pour les valeurs en double
                std::vector<size_t> duplicates;
                for (const auto &[value, count] : counts) {
                    if (count > 1) {
                        duplicates.push_back(value);
                    }
                }

                std::cout << "ERROR POINT:" << std::endl;
                for (const auto &duplicate : duplicates) {
                    std::cout << duplicate << std::endl;
                }
            }
            return false;
        } 

        int MPI_err_indices;
        // We can not send to large data
        if (Np<100000000) {
            std::cout << "USE OF MPI GATHER" << std::endl;
            MPI_err_indices = MPI_Allgatherv(GlobalIndices_LocalProcess.data(), GlobalIndices_LocalProcess.size(), MPI_UNSIGNED_LONG, all_data.data(), sendcounts.data(), displs.data(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
        } else {
            std::cout << "USE OF MPI SEND" << std::endl;
            // Send data from each process to all other processes one by one
            for (int i = 0; i < MPI_size; ++i) {
                if (MPI_rank == i) { // If this is the current sender process
                    // Send data to all other processes
                    for (int j = 0; j < MPI_size; ++j) {
                        if (j != i) { // Don't send to self
                            MPI_Send(GlobalIndices_LocalProcess.data(), GlobalIndices_LocalProcess.size(), MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD);
                        }
                    }
                    // Copy data for the current process
                    std::copy(GlobalIndices_LocalProcess.begin(), GlobalIndices_LocalProcess.end(), all_data.begin() + displs[i]);
                } else { // If this is a receiver process
                    // Calculate the starting point in the receive buffer for the incoming data
                    size_t* recv_buffer_start = all_data.data() + displs[i];
                    MPI_err_indices = MPI_Recv(recv_buffer_start, sendcounts[i], MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }


        if(MPI_err_indices == MPI_SUCCESS) {
            std::cout << "SEND INDICES " << MPI_rank << " : SUCCESS" << std::endl;
        } else {
            std::cout << "SEND INDICES " << MPI_rank << " : FAIL" << std::endl;
        }

        // Creation of GlobalIndices_GlobalProcess
        for (int i = 0; i < MPI_size; i++)
        {
            for (int j = displs[i]; j < displs[i] + sendcounts[i]; j++)
            {
                GlobalIndices_GlobalProcess[all_data[j]] = i;
            }
        }
        
        std::vector<size_t>().swap(GlobalIndices_LocalProcess);
        std::vector<size_t>().swap(all_data);
    }

    return true;
}


void updatePeriodizedSubcube(std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>>& P_with_indices_periodized_subcube,
                             std::vector<Point>& V_with_indices_periodized_subcube,
                             const ArgConfig& ARG_CONFIG, std::size_t Np_Domain) 
{

    std::vector<std::pair<Point, std::pair<std::size_t, std::size_t>>> new_P_with_indices_periodized_subcube;
    std::vector<Point> new_V_with_indices_periodized_subcube;

    std::set<std::size_t> indices_to_keep;

    {
        auto start_time_Delaunay = std::chrono::high_resolution_clock::now();
        
        // Construct the locking data-structure, using the bounding-box of the points
        Triangulation::Lock_data_structure locking_ds( \
            CGAL::Bbox_3(std::get<0>(ARG_CONFIG.DOMAIN_SIZE) - ARG_CONFIG.SLICE_SIZE, \
                        std::get<2>(ARG_CONFIG.DOMAIN_SIZE) - ARG_CONFIG.SLICE_SIZE, \
                        std::get<4>(ARG_CONFIG.DOMAIN_SIZE) - ARG_CONFIG.SLICE_SIZE, \
                        std::get<1>(ARG_CONFIG.DOMAIN_SIZE) + ARG_CONFIG.SLICE_SIZE, \
                        std::get<3>(ARG_CONFIG.DOMAIN_SIZE) + ARG_CONFIG.SLICE_SIZE, \
                        std::get<5>(ARG_CONFIG.DOMAIN_SIZE) + ARG_CONFIG.SLICE_SIZE), \
            64);

        // Construct the triangulation in parallel
        Triangulation T(P_with_indices_periodized_subcube.begin(), P_with_indices_periodized_subcube.end(), &locking_ds);
        assert(T.is_valid());

        auto end_time_Delaunay = std::chrono::high_resolution_clock::now();
        auto duration_Delaunay = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_Delaunay - start_time_Delaunay).count();
        std::cout << "DELAUNAY TIME: " << duration_Delaunay / 1000. << " s" << std::endl;


        auto start_time_remove_halo = std::chrono::high_resolution_clock::now();

        int nCells = 0;
        for(auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end(); ++cell_itr){
            ++nCells;
        }

        std::vector<Triangulation::Cell_handle> cells(nCells);

        int index = 0;
        for(auto cell_itr = T.finite_cells_begin(); cell_itr != T.finite_cells_end(); ++cell_itr){
            cells[index++] = cell_itr;
        }

        #pragma omp parallel
        {
            std::array<std::size_t, 4> indices_array;
            std::set<std::size_t> indices_to_keep_local;

            #pragma omp for
            for(int cell_index = 0; cell_index < nCells; ++cell_index){
                auto cell_itr = cells[cell_index];
                bool has_vertex_with_index_less_than_Np_Domain = false;
                bool has_vertex_with_index_greater_equal_Np_Domain = false;
                
                for(int i=0; i<4; ++i){ // Une cellule a 4 sommets dans une triangulation 3D
                    Triangulation::Vertex_handle vh = cell_itr->vertex(i);
                    std::pair<std::size_t, std::size_t> indices = vh->info();
            
                    indices_array[i] = indices.first;
            
                    if(indices.first < Np_Domain){
                        has_vertex_with_index_less_than_Np_Domain = true;
                    } else {
                        has_vertex_with_index_greater_equal_Np_Domain = true;
                    }
                }
            
                if(has_vertex_with_index_less_than_Np_Domain && has_vertex_with_index_greater_equal_Np_Domain){
                    for(const auto& index : indices_array){
                        if(index >= Np_Domain){
                            indices_to_keep_local.insert(index);
                        }
                    }
                }
            }

            // Merge local results into global set
            #pragma omp critical
            indices_to_keep.insert(indices_to_keep_local.begin(), indices_to_keep_local.end());
        }


        auto end_time_remove_halo = std::chrono::high_resolution_clock::now();
        auto duration_remove_halo = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_remove_halo - start_time_remove_halo).count();
        std::cout << "REMOVE HALO TIME: " << duration_remove_halo / 1000. << " s" << std::endl; // REMOVE HALO TIME: 0.024 s // Np sub-cube slice = 116044
    }

    new_P_with_indices_periodized_subcube.reserve(Np_Domain+indices_to_keep.size());
    new_V_with_indices_periodized_subcube.reserve(Np_Domain+indices_to_keep.size());

    // Copier les Np_Domain premiers points
    for(std::size_t i = 0; i < Np_Domain; ++i){
        new_P_with_indices_periodized_subcube.push_back(P_with_indices_periodized_subcube[i]);
        new_V_with_indices_periodized_subcube.push_back(V_with_indices_periodized_subcube[i]);
    }

    std::size_t current_index = Np_Domain;
    // Ajouter les points qui sont dans indices_to_keep
    for(const auto& idx : indices_to_keep){
        new_P_with_indices_periodized_subcube.push_back(std::make_pair(P_with_indices_periodized_subcube[idx].first, std::make_pair(current_index, P_with_indices_periodized_subcube[idx].second.second)));
        new_V_with_indices_periodized_subcube.push_back(V_with_indices_periodized_subcube[idx]);
        ++current_index;
    }

    // Faire que P_with_indices_periodized_subcube et V_with_indices_periodized_subcube soient égaux à ces nouveaux vecteurs
    P_with_indices_periodized_subcube = new_P_with_indices_periodized_subcube;
    V_with_indices_periodized_subcube = new_V_with_indices_periodized_subcube;
}


void ifLOAD_GRAPH(const ArgConfig& ARG_CONFIG) {

    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

    std::cout << "START GRAPH-WAVELETS" << std::endl;

    std::vector<double> cell_volumes;
    std::vector<double> cell_volumes_rank;

    bool SUCCESS;
    SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol", cell_volumes);
    if (!SUCCESS) return;
    SUCCESS = read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/subcube/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol_" + std::to_string(MPI_rank), cell_volumes_rank);
    if (!SUCCESS) return;
    
    int Np = cell_volumes.size();
    int Np_Domain = cell_volumes_rank.size();
    std::vector<double>().swap(cell_volumes_rank);
    std::cout << "Np = " << Np << " | Np sub-cube = " << Np_Domain << std::endl;

    int8_t num_levels;
    num_levels = static_cast<int>(std::floor(std::log2(Np))); 
    std::cout << "NUMBER OF LEVELS: " << static_cast<int>(num_levels) << std::endl;
    
    Triangulation T;
    std::vector<int16_t> GlobalIndices_GlobalProcess(Np, -1);


    Graph_Wavelets graphWavelets(ARG_CONFIG, T, cell_volumes, GlobalIndices_GlobalProcess, Np, Np_Domain, num_levels);

    MPI_Barrier(MPI_COMM_WORLD);
    
    graphWavelets.writeToBinary(ARG_CONFIG.DATA_PATH.second, ARG_CONFIG.FILENAME, ARG_CONFIG.OUT_SUFFIXES, MPI_rank);

    std::cout << "END GRAPH-WAVELETS" << std::endl;
    std::cout << "_____________________________" << std::endl;
}






















