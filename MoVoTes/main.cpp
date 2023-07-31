/******************************************************************************************************************************
*
*
*
*	Modified Voronoi tesselation and Multi-Scale Tessellation for differential operators of the particle velocity.
*
*
*	Thibault Moli MAUREL OUJIA
*	Original version of the code:
*       - compute the volume cell in parallel using OpenMP for shared memory
*       - compute the divergence and curl of the particle velocity
*
*   
*   Thibault Moli MAUREL OUJIA
*   The following parts of the code were implemented during the JSPS Short-term 2023:
*       - parallel computation of the volume using MPI for distributed memory
*       - wavelets on graphs using MPI for distributed memory
*       
*
*
*   Graph_Wavelets is a C++ version of the developed python code by
*   Keigo MATSUDA and Thibault Moli MAUREL OUJIA during the CTR Summer Program 2022.
*   The implementation of edge_collapse has been adapted to work in parallel with OpenMP.
*
*
*
*   To be done:
*
*   - Added a clear option that deletes all temporary files.
*
*   - Deallocate all vectors after use.
*   - Accept only adjacent domains.
*   - Make it possible again to create level zero in parallel.
*   - Better MPI comunication, send data to ajacent sub-cube.
*
*
*
******************************************************************************************************************************/


#include <mpi.h>
#include <filesystem>


#include "utile.h"
#include "print_data.h"
#include "test_function.h"
#include "modified_voronoi.h"
#include "parallel_delaunay.h"


int main(int argc, char *argv[]) {
    

    // Initialize MPI
    MPI_Init(&argc, &argv);


    // Get the rank of the process and the total number of processes
    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);


    // Read the configuration file
    Config CONFIG;
    std::string config_file_path = "../config.txt";
    if (argc > 1) {
        config_file_path = argv[1];
    }
    std::cout << std::endl << std::endl << "CONFIGURATION FILE" << config_file_path << std::endl << std::endl;
    CONFIG = read_config_file(config_file_path);
    // for (int i = 0; i < MPI_size; i++) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (i == MPI_rank) {
    //         CONFIG = read_config_file("../config.txt");
    //     }
    // }

    

    // Afficher les paramètres de configuration
    std::cout << "Data paths:" << std::endl;
    for (const auto& paths : CONFIG.DATA_PATHS) {
        std::cout << "Input:  " << paths.first << ", Output: " << paths.second << std::endl;
        std::filesystem::create_directories(paths.second);
        std::filesystem::create_directories(paths.second + "/subcube");
    }

    std::cout << "File names:" << std::endl;
    for (const auto& filenames : CONFIG.FILENAMES) {
        for (const auto& filename : filenames) {
            std::cout << filename << ",  ";
        }
            std::cout << std::endl;
    }
    
    std::cout << "Position suffix: " << CONFIG.POSITION_SUFFIX << ", Velicity suffix: " << CONFIG.VELOCITY_SUFFIX << std::endl;

    std::cout << std::endl << "Volume: " << CONFIG.BOOL_VOLUME << std::endl;
    std::cout << "Divergence: " << CONFIG.BOOL_DIVERGENCE << std::endl;
    std::cout << "Curl: " << CONFIG.BOOL_CURL << std::endl;
    std::cout << "Graph: " << CONFIG.BOOL_GRAPH << std::endl;
    std::cout << "Wavelet: " << CONFIG.BOOL_WAVELET << std::endl << std::endl;

    std::cout << "Delta t: " << CONFIG.DELTA_T << std::endl;
    std::cout << "Slice size: " << CONFIG.SLICE_SIZE << std::endl << std::endl;

    std::cout << "Skip Delaunay, Volume and Graph: " << CONFIG.BOOL_SKIP_DELAUNAY_VOLUME_GRAPH << std::endl;

    std::cout << "Volume Merge: " << CONFIG.BOOL_VOLUME_MERGE << std::endl;
    std::cout << "Divergence Merge: " << CONFIG.BOOL_DIVERGENCE_MERGE << std::endl;
    std::cout << "Curl Merge: " << CONFIG.BOOL_CURL_MERGE << std::endl;
    std::cout << "Graph Merge: " << CONFIG.BOOL_GRAPH_MERGE << std::endl;
    std::cout << "Volume Merge: " << CONFIG.BOOL_DIVERGENCE_WRITE << std::endl;
    std::cout << "Divergence Merge: " << CONFIG.BOOL_CURL_WRITE << std::endl;
    std::cout << "Curl Merge: " << CONFIG.BOOL_HELICITY_WRITE << std::endl << std::endl;

    std::cout << "Subdomains: " << CONFIG.SUBDOMAINS << std::endl;
    std::cout << "Load Graph: " << CONFIG.LOAD_GRAPH << std::endl;
    std::cout << "Subcubes: ";
    for(const auto& pair : CONFIG.SUBCUBES) {
        std::cout << pair.first << "-" << pair.second << " ";
    }
    std::cout << std::endl;

    std::cout << "Random: " << CONFIG.BOOL_TEST_RANDOM << std::endl;
    std::cout << "Sequential: " << CONFIG.BOOL_SEQUENTIAL << std::endl << std::endl;

    std::cout << "Domain size: ";
    std::cout << "[" << std::get<0>(CONFIG.DOMAIN_SIZE) << ":";
    std::cout << std::get<1>(CONFIG.DOMAIN_SIZE) << "] x [";
    std::cout << std::get<2>(CONFIG.DOMAIN_SIZE) << ":";
    std::cout << std::get<3>(CONFIG.DOMAIN_SIZE) << "] x [";
    std::cout << std::get<4>(CONFIG.DOMAIN_SIZE) << ":";
    std::cout << std::get<5>(CONFIG.DOMAIN_SIZE) << "]";
    std::cout << std::endl;

    // Largest perfect cube less than or equal to the total number of processes
    int cubeRoot = static_cast<int>(std::cbrt(MPI_size));
    int totalSubdomains = cubeRoot * cubeRoot * cubeRoot;

    // Print the number of processes
    if (MPI_rank == 0) {
        if (CONFIG.SUBDOMAINS == -1 || CONFIG.SUBDOMAINS == 0 || CONFIG.SUBDOMAINS == 1) {
            std::cout << "Number of processes: " << MPI_size << "  ||  Number of sub-domains: " << totalSubdomains << std::endl;
        } else {
            std::cout << "Number of processes: " << MPI_size << "  ||  Number of sub-domains: " << CONFIG.SUBDOMAINS << std::endl;
        }
    }

    // Calculate the number of tasks for each process
    int tasksPerProcess = totalSubdomains / MPI_size;
    int remainingTasks = totalSubdomains % MPI_size;

    // Distribute the tasks among the processes
    int start = MPI_rank * tasksPerProcess + std::min(MPI_rank, remainingTasks);
    int end = start + tasksPerProcess - 1;

    if (MPI_rank < remainingTasks) {
        end++;
    }
    

    #ifdef __APPLE__
        if (MPI_size != 1) omp_set_num_threads(1); // !!!
    #endif


    int num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            std::cout << "TOTAL THREADS: " << num_threads << std::endl;
        }
    }


    ArgConfig ARG_CONFIG;
    ARG_CONFIG.POSITION_SUFFIX = CONFIG.POSITION_SUFFIX;
    ARG_CONFIG.VELOCITY_SUFFIX = CONFIG.VELOCITY_SUFFIX;
    ARG_CONFIG.BOOL_VOLUME = CONFIG.BOOL_VOLUME;
    ARG_CONFIG.BOOL_DIVERGENCE = CONFIG.BOOL_DIVERGENCE;
    ARG_CONFIG.BOOL_CURL = CONFIG.BOOL_CURL;
    ARG_CONFIG.BOOL_GRAPH = CONFIG.BOOL_GRAPH;
    ARG_CONFIG.BOOL_WAVELET = CONFIG.BOOL_WAVELET;

    ARG_CONFIG.DELTA_T = CONFIG.DELTA_T;
    ARG_CONFIG.SLICE_SIZE = CONFIG.SLICE_SIZE;
    
    ARG_CONFIG.BOOL_SKIP_DELAUNAY_VOLUME_GRAPH = CONFIG.BOOL_SKIP_DELAUNAY_VOLUME_GRAPH;
    
    ARG_CONFIG.BOOL_VOLUME_MERGE = CONFIG.BOOL_VOLUME_MERGE;
    ARG_CONFIG.BOOL_DIVERGENCE_MERGE = CONFIG.BOOL_DIVERGENCE_MERGE;
    ARG_CONFIG.BOOL_CURL_MERGE = CONFIG.BOOL_CURL_MERGE;
    ARG_CONFIG.BOOL_GRAPH_MERGE = CONFIG.BOOL_GRAPH_MERGE;
    ARG_CONFIG.LOAD_GRAPH = CONFIG.LOAD_GRAPH;
    
    ARG_CONFIG.BOOL_DIVERGENCE_WRITE = CONFIG.BOOL_DIVERGENCE_WRITE;
    ARG_CONFIG.BOOL_CURL_WRITE = CONFIG.BOOL_CURL_WRITE;
    ARG_CONFIG.BOOL_HELICITY_WRITE = CONFIG.BOOL_HELICITY_WRITE;

    ARG_CONFIG.SUBDOMAINS = CONFIG.SUBDOMAINS;
    ARG_CONFIG.SUBCUBES = CONFIG.SUBCUBES;
    ARG_CONFIG.BOOL_TEST_RANDOM = CONFIG.BOOL_TEST_RANDOM;
    ARG_CONFIG.BOOL_SEQUENTIAL = CONFIG.BOOL_SEQUENTIAL;

    ARG_CONFIG.DOMAIN_SIZE = CONFIG.DOMAIN_SIZE;

    int ID_DATA_PATHS = 0;
    for (const auto& filenames : CONFIG.FILENAMES) {
        ARG_CONFIG.DATA_PATH = CONFIG.DATA_PATHS[ID_DATA_PATHS];
        ID_DATA_PATHS ++;
        for (const auto& filename : filenames) {
            ARG_CONFIG.FILENAME = filename;
            bool SUCCESS;
            if (!ARG_CONFIG.LOAD_GRAPH && !ARG_CONFIG.BOOL_SKIP_DELAUNAY_VOLUME_GRAPH) {
                if (ARG_CONFIG.SUBDOMAINS == -1 || ARG_CONFIG.SUBDOMAINS == 0 || ARG_CONFIG.SUBDOMAINS == 1) {
                    for (int subdomainNumber = start; subdomainNumber <= end; subdomainNumber++) {
                        SUCCESS = parallel_delaunay_multiscale(ARG_CONFIG, subdomainNumber, totalSubdomains);
                        std::this_thread::sleep_for(std::chrono::milliseconds(10));
                    }
                } else {
                    if (MPI_size != 1) {
                        std::cerr << "There has to be 1 process. There are currently " << MPI_size << " processes. "  << std::endl;
                        return 0;
                    }
                    for(const auto& pair : ARG_CONFIG.SUBCUBES) {
                        for(int subdomainNumber = pair.first; subdomainNumber <= pair.second; ++subdomainNumber) {
                            if (subdomainNumber > ARG_CONFIG.SUBDOMAINS-1) {
                                std::cerr << "WARNING: " << subdomainNumber << " > " << ARG_CONFIG.SUBDOMAINS << "-1 --> continue"  << std::endl;
                                continue;
                            } 
                            std::cout << "SUBCUBES NUMBER: " << subdomainNumber << std::flush;
                            parallel_delaunay_multiscale(ARG_CONFIG, subdomainNumber, ARG_CONFIG.SUBDOMAINS);
                        }
                    }
                }
            } else if (!ARG_CONFIG.BOOL_SKIP_DELAUNAY_VOLUME_GRAPH) {
                ifLOAD_GRAPH(ARG_CONFIG);
            }
            
            // Ensure all processes have finished their calculations
            MPI_Barrier(MPI_COMM_WORLD);

            if (MPI_rank == 0) {
                if (ARG_CONFIG.SUBDOMAINS == -1 || ARG_CONFIG.SUBDOMAINS == 0 || ARG_CONFIG.SUBDOMAINS == 1) {
                    merge_data(ARG_CONFIG, totalSubdomains);
                } else {
                    merge_data(ARG_CONFIG, ARG_CONFIG.SUBDOMAINS);
                }
                computes_quantities(ARG_CONFIG);
                if (ARG_CONFIG.BOOL_WAVELET){
                    Graph_Wavelet_Transform(ARG_CONFIG);
                    Graph_Inverse_Wavelet_Transform(ARG_CONFIG);
                    Power_Spectrum(ARG_CONFIG);
                }
            }
        }
    }

    // Finalize MPI
    MPI_Finalize();


    return 0;
}







// Run the code:


// mkdir build
// cd build


// cmake -DCMAKE_BUILD_TYPE=Release .. 
// cmake --build . --config Release 


// cmake -DCMAKE_BUILD_TYPE=Debug ..
// cmake --build . --config Debug 


// mpiexec -n 1 ./parallel_delaunay_multiscale


// clear && cmake --build . --config Release && mpiexec -n 1 ./parallel_delaunay_multiscale







// Profiling:


// advixe-cl --collect=survey --project-dir=./advisor_project -- ./parallel_delaunay_multiscale

// advixe-cl --collect=tripcounts --project-dir=./advisor_project -- ./parallel_delaunay_multiscale

// vtune -collect memory-consumption -result-dir vtune_results -- ./parallel_delaunay_multiscale


// advixe-gui ./advisor_project

// vtune-gui ./vtune_results





// /opt/intel/oneapi/vtune_profiler/latest/MacOS/vtune-gui ./vtune_results








/*


cmake -DCMAKE_BUILD_TYPE=Debug ..


-g : Cette option demande à GCC d'inclure dans le fichier binaire des informations de débogage qui peuvent être utilisées avec un débogueur comme gdb. Ces informations incluent, par exemple, les noms des variables et des fonctions, ainsi que les numéros de ligne.

-O0 : Cette option désactive toutes les optimisations. Cela peut rendre le code plus facile à comprendre pendant le débogage.

-Wall : Cette option demande à GCC de montrer tous les avertissements. C'est très utile pour attraper des problèmes potentiels dans le code.

-Wextra : Cette option active des avertissements supplémentaires qui ne sont pas inclus dans -Wall.

-pedantic : Cette option demande à GCC de respecter strictement la norme C et de signaler toute violation.

-fsanitize=address : C'est une option très utile pour détecter les erreurs d'adressage comme les débordements de tampon et les fuites de mémoire.


*/















    // MPI_Barrier(MPI_COMM_WORLD);
    // int TIMETIMEIME = 1;
    // auto start = std::chrono::high_resolution_clock::now();
    // auto end = start + std::chrono::seconds(TIMETIMEIME);
    // auto duration_totalduration_total = std::chrono::duration_cast<std::chrono::milliseconds>(start - start_time_total).count();
    // std::cout << "TIME: " << duration_totalduration_total / 1000. << " s" << std::endl << std::endl;
    // while (std::chrono::high_resolution_clock::now() < end) {
    //     // Busy wait
    // } // !!!!!!!!!!!!!!!!!!!


    // std::exit(EXIT_FAILURE);






