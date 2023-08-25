#ifndef MULTISCALE_TESSELLATION_H
#define MULTISCALE_TESSELLATION_H


#include <set>
#include <vector>
#include <filesystem>

#include <mpi.h>
#include <stddef.h>

#include "utile.h"
#include "modified_voronoi.h"



class Graph_Wavelets {
public:
    struct Cell {
        double volume = -1.;
        std::set<size_t> adjacent_cells;
        int8_t level = 0;
        int64_t merge = -1;
        bool done = true;
    };
    
    struct MPI_Cell {
        size_t index;
        double volume;
        int64_t merge;
        int done;
        uint64_t level;
        size_t num_adjacent_cells;
        // We cannot transfer a std::set directly.
    };

    Graph_Wavelets(const ArgConfig& ARG_CONFIG, Triangulation& T, const std::vector<double>& cell_volumes, std::vector<int16_t>& _GlobalIndices_GlobalProcess, size_t Np, size_t Np_Domain, int8_t num_levels, size_t none_periodic = -1);

    void writeToBinary(const std::string& directory, const std::string& filename, const std::string& OUT_SUFFIXES, int number = -1);
    void loadFromBinary(const std::string& directory, const std::string& filename, const std::string& OUT_SUFFIXES, size_t Np, size_t Np_Domain, int number = -1);

    std::set<size_t> cells_indices;
    std::set<size_t> send_indices_local;
    std::set<size_t> send_indices_global;
    std::vector<size_t> previous_level_indices;
    std::vector<size_t> current_level_indices;
    std::vector<int16_t> GlobalIndices_GlobalProcess;


private:
    Triangulation T;
    std::unordered_map<size_t, Cell> cells; 

    void createLevelZero_Graph(const std::vector<double>& cell_volumes, size_t Np, size_t Np_Domain, size_t none_periodic = -1);
    void Upscale_Graph(const ArgConfig& ARG_CONFIG, int8_t num_levels);

    void SequentialExecution(std::vector<std::pair<double, size_t>>& volume_index_pairs, int MPI_size, int MPI_rank);
    void MPI_Cells(std::vector<std::pair<double, size_t>>& volume_index_pairs, int MPI_size, int MPI_rank);
    void Merge_Cells(const ArgConfig& ARG_CONFIG, std::vector<std::pair<double, size_t>>& volume_index_pairs, int8_t num_levels, int MPI_size, int MPI_rank); 
    void Link_Cells(std::vector<std::pair<double, size_t>>& volume_index_pairs, int8_t num_levels); 

    void processNeighborCells(Cell& current_cell, std::set<size_t>& new_adjacent_cells, int8_t num_levels);

    void printCell(size_t N) const;
    void printCell_bis(size_t N) const;
    void printFirstNCells(size_t N) const;
    void printAdjacentCellDetails(size_t index) const;

};



void Graph_Wavelet_Transform(const ArgConfig& ARG_CONFIG, int8_t num_levels = -1);
void Graph_Inverse_Wavelet_Transform(const ArgConfig& ARG_CONFIG, int8_t num_levels = -1);
void Power_Spectrum(const ArgConfig& ARG_CONFIG, int num_levels = -1);







#endif // MULTISCALE_TESSELLATION_H






















