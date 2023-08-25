#include "multiscale_tessellation.h"



/*******************************************************************
*				            GRAPH-WAVELETS 		        		   *
*******************************************************************/


Graph_Wavelets::Graph_Wavelets(const ArgConfig& ARG_CONFIG, Triangulation& T, const std::vector<double>& cell_volumes, std::vector<int16_t>& _GlobalIndices_GlobalProcess, size_t Np, size_t Np_Domain, int8_t num_levels, size_t none_periodic) : T(T), GlobalIndices_GlobalProcess(std::move(_GlobalIndices_GlobalProcess)) {
    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
    if (!ARG_CONFIG.LOAD_GRAPH){
        createLevelZero_Graph(cell_volumes, Np, Np_Domain, none_periodic);
        T.clear();
        T = Triangulation();
    } else {
        loadFromBinary(ARG_CONFIG.DATA_PATH.second, ARG_CONFIG.FILENAME, ARG_CONFIG.OUT_SUFFIXES, Np, Np_Domain, MPI_rank);
    }
    
    for (int8_t i = 1; i < num_levels; ++i) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        Upscale_Graph(ARG_CONFIG, i);
    }
}



/*******************************************************************
*				             LEVEL-ZERO 		        		   *
*******************************************************************/


void Graph_Wavelets::createLevelZero_Graph(const std::vector<double>& cell_volumes, size_t Np, size_t Np_Domain, size_t none_periodic) {

    int num_vertices = Np_Domain;
    std::atomic<int> progress_count(0);
    std::vector<bool> progress_markers(11, false);

    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "CREATE LEVEL 0: |0" << std::flush;
    progress_markers[0] = true;


    cells.reserve(2*Np_Domain);
    //omp_set_num_threads(8);


    // int nVertices = 0;
    // for(auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
    //     ++nVertices;
    // }

    // std::vector<Triangulation::Vertex_handle> vertices(nVertices);
    // std::vector<std::vector<Triangulation::Cell_handle>> all_incident_cells(nVertices);

    // int index = 0;
    // for(auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
    //     vertices[index] = vit;
        
    //     std::vector<Triangulation::Cell_handle> incident_cells;
    //     T.finite_incident_cells(vit, std::back_inserter(incident_cells)); // not thread safe, why ?
    //     all_incident_cells[index] = std::move(incident_cells);

    //     ++index;
    // }

    // int main_cells = 0;

    // #pragma omp parallel
    // {
    //     std::vector<size_t> previous_level_indices_local;

    //     #pragma omp for reduction(+:main_cells)
    //     for(int vit_index = 0; vit_index < nVertices; ++vit_index){
    //         auto vit = vertices[vit_index];
            
    //         std::pair<std::size_t, std::size_t> idx1 = vit->info();
    //         size_t index;
    //         size_t index_local;
    //         if (none_periodic==-1){
    //             index = idx1.second;
    //             index_local = idx1.first;
    //         } else {
    //             index = idx1.first;
    //             index_local = idx1.first;
    //         }
            
    //         if (index_local >= Np_Domain || index_local >= cell_volumes.size()) { // Skip indices out of range
    //             continue;  
    //         }

    //         // Cell& cell = cells[index];
    //         Cell& cell = cells.emplace(index, Cell()).first->second; // need prealloc for parallel
            
    //         main_cells ++;
    //         cell.done = false;
    //         cell.volume = cell_volumes[index_local];
    //         previous_level_indices_local.push_back(index);
            
    //         const auto& incident_cells = all_incident_cells[vit_index];
            
    //         for (const auto& incident_cell : incident_cells) {
    //             for (int i = 0; i < 4; ++i) {
    //                 std::pair<std::size_t, std::size_t> idx2 = incident_cell->vertex(i)->info();
    //                 size_t neighbor_index;
    //                 if (none_periodic==-1){
    //                     neighbor_index = idx2.second;
    //                 } else {
    //                     neighbor_index = idx2.first;
    //                 }
    //                 if (neighbor_index != index && neighbor_index < Np) {
    //                     cell.adjacent_cells.insert({neighbor_index});
    //                 }
    //             }
    //         }
            
    //         int count = ++progress_count;
    //         int progress = count * 100 / num_vertices;
    //         int progress_index = progress / 10;
    //         if (progress_index < progress_markers.size()) {
    //             if (!progress_markers[progress_index]) {
    //                 #pragma omp critical
    //                 {
    //                     if (!progress_markers[progress_index]) {
    //                         std::cout << "|" << progress << std::flush;
    //                         progress_markers[progress_index] = true;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     #pragma omp critical
    //     previous_level_indices.insert(previous_level_indices.end(), previous_level_indices_local.begin(), previous_level_indices_local.end());
    // }


    int main_cells = 0;
    
    //#pragma omp parallel for // parallel not possible with iterator
    for (auto it = T.finite_vertices_begin(); it != T.finite_vertices_end(); ++it) {
        
        std::pair<std::size_t, std::size_t> idx1 = it->info();
        size_t index;
        size_t index_local;
        if (none_periodic==-1){
            index = idx1.second;
            index_local = idx1.first;
        } else {
            index = idx1.first;
            index_local = idx1.first;
        }
        
        
        if (index_local >= Np_Domain || index_local >= cell_volumes.size()) { // Skip indices out of range
            continue;  
        }

        // Cell& cell = cells[index];
        Cell& cell = cells.emplace(index, Cell()).first->second;
        cells_indices.insert(index);
        
        main_cells ++;
        cell.done = false;
        cell.volume = cell_volumes[index_local];
        previous_level_indices.push_back(index);
        

        std::vector<Triangulation::Cell_handle> incident_cells;
        T.finite_incident_cells(it, std::back_inserter(incident_cells));
        
        for (const auto& incident_cell : incident_cells) {
            for (int i = 0; i < 4; ++i) {
                std::pair<std::size_t, std::size_t> idx2 = incident_cell->vertex(i)->info();
                size_t neighbor_index;
                if (none_periodic==-1){
                    neighbor_index = idx2.second;
                } else {
                    neighbor_index = idx2.first;
                }
                if (neighbor_index != index && neighbor_index < Np) {
                    cell.adjacent_cells.insert(neighbor_index);
                }
            }
        }
        
        size_t count = ++progress_count;
        size_t progress = count * 100 / num_vertices;
        size_t progress_index = progress / 10;
        if (progress_index < progress_markers.size()) {
            if (!progress_markers[progress_index]) {
                std::cout << "|" << progress << std::flush;
                progress_markers[progress_index] = true;
            }
        }
    }


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|  Time: " << duration / 1000. << " s" << std::endl;
    
    std::cout << "SIZE CELLS: " << cells.size() << "  ||  MAIN CELLS: " << main_cells << std::endl;
}



/*******************************************************************
*				            UPSCALE GRAPH 		        		   *
*******************************************************************/


void Graph_Wavelets::Upscale_Graph(const ArgConfig& ARG_CONFIG, int8_t num_levels) {

    std::cout << "UPSCALE GRAPH TO LEVEL " << static_cast<int>(num_levels) << ":" << std::endl;


    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);
    MPI_Barrier(MPI_COMM_WORLD);


    // Creating a table of volumes and indices
    current_level_indices.clear();
    std::vector<std::pair<double, size_t>> volume_index_pairs;
    for (const auto& i : previous_level_indices) {
        if (cells[i].level == num_levels - 1 && cells[i].level != -1) {
            volume_index_pairs.push_back({cells[i].volume, i});
            current_level_indices.push_back(i); // Add the index to current_level_indices
        }
    }
    
    // int num_cells = volume_index_pairs.size(); // True table size

    // Update previous_level_indices for the next iteration
    previous_level_indices = current_level_indices;

    // Sorting the volume table while retaining the original indices
    std::sort(volume_index_pairs.begin(), volume_index_pairs.end());


    MPI_Cells(volume_index_pairs, MPI_size, MPI_rank);
    

    if (ARG_CONFIG.BOOL_SEQUENTIAL) { 
        SequentialExecution(volume_index_pairs, MPI_size, MPI_rank);
    } 

    
    Merge_Cells(ARG_CONFIG, volume_index_pairs, num_levels, MPI_size, MPI_rank);


    Link_Cells(volume_index_pairs, num_levels);


    MPI_Barrier(MPI_COMM_WORLD);
    
    
    // Sets .done to false
    int main_cells = 0;
    for (const auto& pair : volume_index_pairs) {
        size_t current_index = pair.second;
        if (current_index == -1) {
            continue; 
        }
        cells[current_index].done = false;
        if (cells[current_index].level == num_levels && pair.first != -1){
            main_cells ++;
        } else {
            cells[current_index].adjacent_cells.clear();
        }
    }
    int num_cells = volume_index_pairs.size();
    std::cout << "MAIN CELLS: " << num_cells << " --> " << main_cells << std::endl;

    // Remove useless cells
    for (auto it = cells.begin(); it != cells.end();) {
        if (cells_indices.find(it->first) == cells_indices.end()) {
            // The key is not in cells_indices, so erase this element
            it = cells.erase(it);
        } else {
            ++it;
        }
    }
}


void Graph_Wavelets::SequentialExecution(std::vector<std::pair<double, size_t>>& volume_index_pairs, int MPI_size, int MPI_rank){
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "SEQ:   |0" << std::flush;

    struct PairStruct {
        double first;
        size_t second;
    };

    // Define the custom MPI data type for std::pair<double, size_t>.
    MPI_Datatype MPI_Pair;
    MPI_Datatype types[2] = { MPI_DOUBLE, MPI_UNSIGNED_LONG };
    int blocklengths[2] = { 1, 1 };
    MPI_Aint offsets[2];

    offsets[0] = offsetof(PairStruct, first);
    offsets[1] = offsetof(PairStruct, second);

    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_Pair);
    MPI_Type_commit(&MPI_Pair);

    std::cout << "|10" << std::flush;

    // Calculates the amount of data to be sent for each process
    std::vector<int> sendcounts(MPI_size, 0), displs(MPI_size, 0);
    sendcounts[MPI_rank] = volume_index_pairs.size();

    // Gathers the send accounts of all processes
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, sendcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::cout << "|20" << std::flush;

    // Calculates the displacements for each process
    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + sendcounts[i-1];
    }

    // Calculates the total size for volume_index_pairs_all_process
    size_t total_size = 0;
    for (size_t count : sendcounts)
    {
        total_size += count;
    }

    std::cout << "|30" << std::flush;


    std::vector<std::pair<double, size_t>> volume_index_pairs_all_process(total_size);

    // Gather all data in volume_index_pairs_all_process
    // MPI_Allgatherv(volume_index_pairs.data(), volume_index_pairs.size(), MPI_Pair,
    //             volume_index_pairs_all_process.data(), sendcounts.data(), displs.data(), MPI_Pair, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (total_size < 100000) {
        MPI_Allgatherv(volume_index_pairs.data(), volume_index_pairs.size(), MPI_Pair,
                    volume_index_pairs_all_process.data(), sendcounts.data(), displs.data(), MPI_Pair, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < MPI_size; ++i) {
            if (MPI_rank == i) { // If this is the current sender process
                // Send data to all other processes
                for (int j = 0; j < MPI_size; ++j) {
                    if (j != i) { // Don't send to self
                        MPI_Send(volume_index_pairs.data(), volume_index_pairs.size(), MPI_Pair, j, 0, MPI_COMM_WORLD);
                    }
                }
                // Copy data for the current process
                std::copy(volume_index_pairs.begin(), volume_index_pairs.end(), volume_index_pairs_all_process.begin() + displs[i]);
            } else { // If this is a receiver process
                // Calculate the starting point in the receive buffer for the incoming data
                auto recv_buffer_start = volume_index_pairs_all_process.data() + displs[i];
                MPI_Recv(recv_buffer_start, sendcounts[i], MPI_Pair, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    // Cleaning
    MPI_Type_free(&MPI_Pair);

    std::cout << "|40" << std::flush;
    
    // std::sort(volume_index_pairs_all_process.begin(), volume_index_pairs_all_process.end());

    std::vector<std::pair<double, size_t>> new_volume_index_pairs_all_process;
    new_volume_index_pairs_all_process.reserve(volume_index_pairs_all_process.size());

    for (const auto& pair : volume_index_pairs_all_process) {
        if (pair.second != static_cast<size_t>(-1)) {
            new_volume_index_pairs_all_process.push_back(pair);
        }
    }

    std::cout << "|50" << std::flush;

    volume_index_pairs_all_process = std::move(new_volume_index_pairs_all_process);

    std::cout << "|60" << std::flush;

    std::sort(volume_index_pairs_all_process.begin(), volume_index_pairs_all_process.end());

    std::cout << "|70" << std::flush;

    // Transform volume_index_pairs into a set for quick lookup
    std::unordered_set<size_t> local_indices;
    for (const auto& pair : volume_index_pairs) {
        local_indices.insert(pair.second);
    }

    std::cout << "|80" << std::flush;


    // Check each element of volume_index_pairs_all_process
    for (auto& pair : volume_index_pairs_all_process) {
        // If the second element is NOT in local_indices, set the first element to -1
        if (local_indices.count(pair.second) == 0) {
            pair.first = -1;
        }
    }

    std::cout << "|90" << std::flush;
    
    volume_index_pairs = volume_index_pairs_all_process; 

    
    std::cout << "|100" << std::flush;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|  Time: " << duration / 1000. << " s" << std::endl;
}


void Graph_Wavelets::MPI_Cells(std::vector<std::pair<double, size_t>>& volume_index_pairs, int MPI_size, int MPI_rank) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "MPI:   |0" << std::flush;

    
    // Determine the maximum size across all processes.
    int max_size = 0;
    int local_size = volume_index_pairs.size();
    MPI_Allreduce(&local_size, &max_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // Extend volume_index_pairs to max_size with "empty" pairs.
    if (volume_index_pairs.size() < max_size) {
        // volume_index_pairs.resize(max_size, {-1, 0}); // Don't use 0. The particle ID 0 already exists.
        volume_index_pairs.resize(max_size, {-1, -1});
    }
    
    std::cout << "|10" << std::flush;

    std::vector<size_t> neighbor_data_to_transfer;

    // Volume indices in ascending order
    for (const auto& pair : volume_index_pairs) {
        // Get the index and the current cell
        size_t current_index = pair.second;
        if (pair.first == -1) { // Skip "empty" pairs. 
            continue; 
        }
        Cell& current_cell = cells[current_index];

        // Browse adjacent cells to check their process rank
        for (const auto& neighbor : current_cell.adjacent_cells) {
            if (GlobalIndices_GlobalProcess[neighbor] != MPI_rank) { // If the neighbour's rank is not the same as the current rank
                neighbor_data_to_transfer.push_back(neighbor); // Add the neighbour's index to the list to be transferred
            }
        }
    }

    std::cout << "|20" << std::flush;

    std::vector<int> sendcounts(MPI_size, 0), displs(MPI_size, 0);

    // Determines how much data each process should send
    sendcounts[MPI_rank] = neighbor_data_to_transfer.size();

    // Gathers send accounts from all processes
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, sendcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Calculates the total size needed for neighbor_indices_transfered
    int total_size = 0;
    local_size = neighbor_data_to_transfer.size();
    MPI_Allreduce(&local_size, &total_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + sendcounts[i-1];
    }

    // Gathers all data in neighbor_data_received
    std::vector<size_t> neighbor_data_received(total_size);
    MPI_Allgatherv(neighbor_data_to_transfer.data(), neighbor_data_to_transfer.size(), MPI_UNSIGNED_LONG,
                neighbor_data_received.data(), sendcounts.data(), displs.data(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    // Remove all duplicates
    send_indices_global.clear();
    send_indices_global.insert(neighbor_data_received.begin(), neighbor_data_received.end());
    neighbor_data_received.clear();

    std::cout << "|30" << std::flush;


    MPI_Datatype MPI_CellData;
    MPI_Datatype types_cell[6] = { MPI_UNSIGNED_LONG, MPI_DOUBLE, MPI_LONG, MPI_INT, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG };
    int blocklengths_cell[6] = { 1, 1, 1, 1, 1, 1 };
    MPI_Aint offsets_cell[6];

    offsets_cell[0] = offsetof(MPI_Cell, index);
    offsets_cell[1] = offsetof(MPI_Cell, volume);
    offsets_cell[2] = offsetof(MPI_Cell, merge);
    offsets_cell[3] = offsetof(MPI_Cell, done);
    offsets_cell[4] = offsetof(MPI_Cell, level);
    offsets_cell[5] = offsetof(MPI_Cell, num_adjacent_cells);

    MPI_Type_create_struct(6, blocklengths_cell, offsets_cell, types_cell, &MPI_CellData);
    MPI_Type_commit(&MPI_CellData);

    std::vector<MPI_Cell> cells_to_transfer;
    for (const auto& data : send_indices_global) {
        if (GlobalIndices_GlobalProcess[data] == MPI_rank) {
            Cell& current_cell = cells[data];
            
            MPI_Cell mpi_cell;
            mpi_cell.index = data;
            mpi_cell.volume = current_cell.volume;
            mpi_cell.level = current_cell.level;
            mpi_cell.merge = current_cell.merge;
            mpi_cell.done = current_cell.done ? 1 : 0;
            mpi_cell.num_adjacent_cells = current_cell.adjacent_cells.size();
            cells_to_transfer.push_back(mpi_cell);
        }
    }

    std::cout << "|40" << std::flush;

    sendcounts[MPI_rank] = cells_to_transfer.size();

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, sendcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    total_size = 0;
    local_size = cells_to_transfer.size();
    MPI_Allreduce(&local_size, &total_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + sendcounts[i-1];
    }

    std::cout << "|50" << std::flush;

    std::vector<MPI_Cell> cells_received(total_size);
    MPI_Allgatherv(cells_to_transfer.data(), cells_to_transfer.size(), MPI_CellData, cells_received.data(), sendcounts.data(), displs.data(), MPI_CellData, MPI_COMM_WORLD);

    MPI_Type_free(&MPI_CellData);
    
    std::cout << "|60" << std::flush;


    std::vector<size_t> adj_cells_to_transfer;

    for (const auto& cell : cells_to_transfer) {
        for (const auto& adj_cell : cells[cell.index].adjacent_cells) { 
            adj_cells_to_transfer.push_back(adj_cell);
        }
    }
    cells_to_transfer.clear();

    std::cout << "|70" << std::flush;

    sendcounts[MPI_rank] = adj_cells_to_transfer.size();

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, sendcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    total_size = 0;
    local_size = adj_cells_to_transfer.size();
    MPI_Allreduce(&local_size, &total_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 1; i < MPI_size; i++)
    {
        displs[i] = displs[i-1] + sendcounts[i-1];
    }

    // std::cout << " total_size " << total_size << std::flush;
    std::vector<size_t> adj_cells_received(total_size); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // std::cout << " total_size " << total_size << std::flush;
    // std::vector<size_t> adj_cells_received;
    // adj_cells_received.reserve(total_size);
    // size_t logStep = 1;
    // try {
    //     for (size_t i = 0; i < total_size; ++i) {
    //         adj_cells_received.push_back(0);

    //         // Check if i is a power of 10 and print it
    //         if (i == logStep) {
    //             std::cout << "Reached index: " << i << std::endl;
    //             logStep *= 10;
    //         }
    //     }
    // } catch (std::exception& e) {
    //     std::cerr << "Error occurred at index: " << adj_cells_received.size() << "\nError message: " << e.what() << std::endl;
    // }
    // std::cout << " total_size " << total_size << std::flush;


    MPI_Barrier(MPI_COMM_WORLD);
    if (total_size < 100000000) {
        std::cout << "|80" << std::flush;
        MPI_Allgatherv(adj_cells_to_transfer.data(), adj_cells_to_transfer.size(), MPI_UNSIGNED_LONG, adj_cells_received.data(), sendcounts.data(), displs.data(), MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    } else {
        std::cout << "|8O" << std::flush;
        for (int i = 0; i < MPI_size; ++i) {
            // std::cout << "   i   " << i << std::flush; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (MPI_rank == i) { // If this is the current sender process
                // Send data to all other processes
                for (int j = 0; j < MPI_size; ++j) {
                    if (j != i) { // Don't send to self
                        MPI_Send(adj_cells_to_transfer.data(), adj_cells_to_transfer.size(), MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD);
                    }
                }
                // Copy data for the current process
                std::copy(adj_cells_to_transfer.begin(), adj_cells_to_transfer.end(), adj_cells_received.begin() + displs[i]);
            } else { // If this is a receiver process
                // Calculate the starting point in the receive buffer for the incoming data
                size_t* recv_buffer_start = adj_cells_received.data() + displs[i];
                MPI_Recv(recv_buffer_start, sendcounts[i], MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    adj_cells_to_transfer.clear();

    std::cout << "|90" << std::flush;


    send_indices_local.clear();
    send_indices_local.insert(neighbor_data_to_transfer.begin(), neighbor_data_to_transfer.end());
    size_t current_adj_cell_index = 0;
    for (size_t i = 0; i < cells_received.size(); i++) {
        const auto& received_cell = cells_received[i];

        // Check if received_cell.index is in send_indices_local
        if(send_indices_local.find(received_cell.index) != send_indices_local.end()) {
            // Update cell information
            Cell& cell = cells.emplace(received_cell.index, Cell()).first->second;

            cell.volume = received_cell.volume;
            cell.merge = received_cell.merge;
            cell.level = received_cell.level;
            cell.done = received_cell.done == 1;
            cell.adjacent_cells.clear(); 

            // Adds cells adjacent to the cell
            for (size_t j = 0; j < received_cell.num_adjacent_cells; j++) {
                size_t adj_cell_index = adj_cells_received[current_adj_cell_index++];
                cell.adjacent_cells.insert(adj_cell_index);
            }
        } else {
            // If received_cell.index is not in send_indices_local, just advance current_adj_cell_index
            current_adj_cell_index += received_cell.num_adjacent_cells;
        }
    }
    cells_received.clear();
    adj_cells_received.clear();

    std::cout << "|100" << std::flush;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|  Time: " << duration / 1000. << " s" << std::endl;
}


void Graph_Wavelets::Merge_Cells(const ArgConfig& ARG_CONFIG, std::vector<std::pair<double, size_t>>& volume_index_pairs, int8_t num_levels, int MPI_size, int MPI_rank) {
    int num_cells = volume_index_pairs.size();
    
    // bool detailed_progress = true; // true // false // !!! 
    bool detailed_progress = (num_cells > 1000000000) ? true : false;
    std::atomic<int> progress_count(0);
    int num_markers = detailed_progress ? 101 : 11; // We need 101 markers for 1% increments, 11 for 10%.
    std::vector<bool> progress_markers(num_markers, false);

    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "MERGE: |0" << std::flush;
    progress_markers[0] = true;

    // Volume indices in ascending order
    for (const auto& pair : volume_index_pairs) {
        
        // Initialized with size 2*MPI_size, all elements set to -1.
        std::vector<int> merged_cells(2 * MPI_size, -1);
        std::vector<double> merged_volumes(2 * MPI_size, -1);

        // Get the index and the current cell
        size_t current_index = pair.second;
        
        auto it = cells.find(current_index); // Find cell in the map
        if (it != cells.end()) { // If the cell exists
            Cell& current_cell = cells[current_index];

            // If the cell has not been merged
            if (pair.first != -1 && !current_cell.done) { 
                // Initialisation of the minimum neighbouring cell
                std::pair<double, size_t> min_neighbor = {std::numeric_limits<double>::max(), 0};

                // Browse adjacent cells to find the one with the smallest volume
                for (const auto& neighbor : current_cell.adjacent_cells) {
                    
                    size_t neighbor_index = neighbor;
                    const Cell& neighbor_cell = cells[neighbor_index];

                    // Check that the neighbouring cell has not already been merged and find the smallest volume 
                    if (!neighbor_cell.done && neighbor_cell.volume < min_neighbor.first) { 
                        min_neighbor = {neighbor_cell.volume, neighbor_index};
                    }
                }

                if (min_neighbor.first != std::numeric_limits<double>::max()) {
                    // Check that a minimum neighbouring cell has been found

                    // Selection of the index of the largest cell by volume between current_cell and min_neighbor
                    size_t even = (current_cell.volume > min_neighbor.first) ? current_index : min_neighbor.second; // even (max_cell_index) : large cell
                    size_t odd = (current_cell.volume > min_neighbor.first) ? min_neighbor.second : current_index;  // odd  (min_cell_index) : small cell


                    // We check if we need to transfer this datas
                    bool condition_neighbor_of_neighbor = false;

                    // Test for cells[odd]
                    for(const auto& adj : cells[odd].adjacent_cells){
                        if(send_indices_global.count(adj) > 0){
                            // Information to be shared among all processes 
                            merged_cells[2*MPI_rank] = static_cast<int>(even); // Set the 2*MPI_rank position to even. 
                            merged_cells[2*MPI_rank + 1] = static_cast<int>(odd); // Set the 2*MPI_rank+1 position to odd. 
                            merged_volumes[2*MPI_rank] = cells[even].volume; 
                            merged_volumes[2*MPI_rank + 1] = cells[odd].volume; 
                            condition_neighbor_of_neighbor = true;
                            break;
                        }
                    }

                    // If condition is not met in cells[odd], check cells[even]
                    if(!condition_neighbor_of_neighbor) {
                        for(const auto& adj : cells[even].adjacent_cells){
                            if(send_indices_global.count(adj) > 0){
                                // Information to be shared among all processes 
                                merged_cells[2*MPI_rank] = static_cast<int>(even); // Set the 2*MPI_rank position to even. 
                                merged_cells[2*MPI_rank + 1] = static_cast<int>(odd); // Set the 2*MPI_rank+1 position to odd. 
                                break;
                            }
                        }
                    }


                    // Update of .done and .merge for the two cells in cells
                    cells[odd].done = true;
                    cells[odd].merge = even;

                    cells[even].level = num_levels;
                    cells[even].done = true;
                    cells[even].merge = odd;
                    cells[even].volume = cells[even].volume + cells[odd].volume; 
                }
                else {
                    // All neighbours have already been merged

                    // We check if we need to transfer this datas
                    for(const auto& adj : current_cell.adjacent_cells){
                        if(send_indices_global.count(adj) > 0){
                            merged_cells[2*MPI_rank] = static_cast<int>(current_index); // Set the 2*MPI_rank position to current_index. 
                            merged_cells[2*MPI_rank + 1] = static_cast<int>(current_index); // Set the 2*MPI_rank+1 position to current_index. 
                            merged_volumes[2*MPI_rank] = cells[current_index].volume; 
                            merged_volumes[2*MPI_rank + 1] = cells[current_index].volume; 
                            break;
                        }
                    }

                    current_cell.level = num_levels;
                    current_cell.done = true;
                    current_cell.merge = current_index;
                }
            }
        }
        // Sending and receiving merge information using MPI_Allgather with MPI_IN_PLACE
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, merged_cells.data(), 2, MPI_INT, MPI_COMM_WORLD); 
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, merged_volumes.data(), 2, MPI_DOUBLE, MPI_COMM_WORLD); 

        // Cells updated on the basis of merge information received
        for (int i = 0; i < MPI_size; ++i) {
            int sent_cell_1 = merged_cells[2 * i];
            int sent_cell_2 = merged_cells[2 * i + 1];

            if (sent_cell_1 != -1 && sent_cell_2 != -1 && i != MPI_rank) {  // Ignore if the value is -1, i.e. no merge, or if it is the current process, so as not to add the same volume twice.
                cells[sent_cell_1].level = num_levels;
                cells[sent_cell_1].done = true;
                cells[sent_cell_1].merge = sent_cell_2;
                // cells[sent_cell_1].volume = cells[sent_cell_1].volume + cells[sent_cell_2].volume; 
                cells[sent_cell_1].volume = merged_volumes[2 * i] + merged_volumes[2 * i + 1]; 
                
                cells[sent_cell_2].done = true;
                cells[sent_cell_2].merge = sent_cell_1;


                // We could look at all the neighbouring domains (at the same time as neighbor_data_to_transfer) and accept only the cells that belong to the neighbouring domain.
                // Check if sent_cell_1 or sent_cell_2 are in send_indices_local // Wrong, but it's an idea to reduce the number of cells in memory.
                // if (send_indices_local.find(sent_cell_1) != send_indices_local.end() || send_indices_local.find(sent_cell_2) != send_indices_local.end()) {
                //     ...
                // }
            }
        }


        if (!ARG_CONFIG.BOOL_SEQUENTIAL) { 
            // Créer un map pour stocker la liste des processus pour chaque indice.
            std::map<int, std::vector<int>> index_process_map;
            for (int i = 0; i < MPI_size; ++i) {
                if (merged_cells[2 * i] != merged_cells[2 * i + 1]) { 
                    index_process_map[merged_cells[2 * i]].push_back(i);
                    index_process_map[merged_cells[2 * i + 1]].push_back(i);
                }
            }

            // Parcourir la map et pour chaque indice qui apparaît plusieurs fois, vérifier les processus.
            for (const auto& pair : index_process_map) {
                if (pair.second.size() > 1) { // Si l'indice apparaît plus d'une fois.
                    double min_sum = std::numeric_limits<double>::max();
                    int min_process = -1;
                    // Parcourir tous les processus pour cet indice.
                    for (int i : pair.second) {
                        double sum = merged_volumes[2 * i] + merged_volumes[2 * i + 1];
                        if (sum < min_sum) {
                            min_sum = sum;
                            min_process = i;
                        }
                    }
                    // Comparer min_process avec chaque processus.
                    for (int i = 0; i < MPI_size; ++i) {
                        if (min_process != i) {  // We do nothing for min_process == i 
                            if (merged_cells[2 * min_process] == merged_cells[2 * i] || merged_cells[2 * min_process] == merged_cells[2 * i + 1]
                                || merged_cells[2 * min_process + 1] == merged_cells[2 * i] || merged_cells[2 * min_process + 1] == merged_cells[2 * i + 1]) {
                                
                                cells[merged_cells[2 * i]].level = num_levels;
                                cells[merged_cells[2 * i]].done = true;
                                cells[merged_cells[2 * i]].merge = merged_cells[2 * i];
                                cells[merged_cells[2 * i]].volume = merged_volumes[2 * i];

                                cells[merged_cells[2 * i + 1]].level = num_levels;
                                cells[merged_cells[2 * i + 1]].done = true;
                                cells[merged_cells[2 * i + 1]].merge = merged_cells[2 * i + 1];
                                cells[merged_cells[2 * i + 1]].volume = merged_volumes[2 * i + 1];

                                // std::cout << std::endl << "MPI_rank " << MPI_rank << " [2 * MPI_rank] " << merged_cells[2 * MPI_rank] << " [2 * MPI_rank + 1] " << merged_cells[2 * MPI_rank + 1] << " [2 * min_process] " << merged_cells[2 * min_process] << " [2 * min_process + 1] " << merged_cells[2 * min_process + 1] << std::endl;
                            }
                        }
                    }
                }
            }
        } 
        // Reset merged_cells to -1 after each merge operation
        std::fill(merged_cells.begin(), merged_cells.end(), -1);

        size_t count = ++progress_count;
        size_t progress = count * 100 / num_cells;
        size_t progress_index = detailed_progress ? progress : progress / 10;
        if (progress_index < progress_markers.size()) {
            if (!progress_markers[progress_index]) {
                #pragma omp critical
                {
                    if (!progress_markers[progress_index]) {
                        std::cout << "|" << progress << std::flush;
                        progress_markers[progress_index] = true;
                    }
                }
            }
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "|  Time: " << duration / 1000. << " s" << std::endl;
}


void Graph_Wavelets::Link_Cells(std::vector<std::pair<double, size_t>>& volume_index_pairs, int8_t num_levels) {
    int num_cells = volume_index_pairs.size();

    // bool detailed_progress = false; // true // false
    bool detailed_progress = (num_cells > 1000000000) ? true : false;
    std::atomic<int> progress_count(0);
    int num_markers = detailed_progress ? 101 : 11; // We need 101 markers for 1% increments, 11 for 10%.
    std::vector<bool> progress_markers(num_markers, false);


    auto start_time_ = std::chrono::high_resolution_clock::now();
    std::cout << "LINK:  |0" << std::flush;
    progress_markers[0] = true;


    #pragma omp parallel for
    for (const auto& pair : volume_index_pairs) {

        size_t current_index = pair.second;
        auto it = cells.find(current_index); // Find cell in the map
        if (it != cells.end()) {
            Cell& current_cell = cells[current_index];

            // If the cell is in the new level
            if (pair.first != -1 && current_cell.level == num_levels) { 
                
                // Create a set to store new neighbours
                std::set<size_t> new_adjacent_cells;
                
                // Examine the neighbours of the current cell and the cell with which it has been merged
                processNeighborCells(current_cell, new_adjacent_cells, num_levels);

                // The same is done for the neighbours of the cell with which the current cell has been merged
                processNeighborCells(cells[current_cell.merge], new_adjacent_cells, num_levels);

                // Remove it self and merge
                new_adjacent_cells.erase(current_index);
                new_adjacent_cells.erase(current_cell.merge);
                new_adjacent_cells.erase(static_cast<size_t>(-1)); // !!! probleme non resolu

                // Update the neighbours of the current cell
                current_cell.adjacent_cells = new_adjacent_cells;
            }
        }

        size_t count = ++progress_count;
        size_t progress = count * 100 / num_cells;
        size_t progress_index = detailed_progress ? progress : progress / 10;
        if (progress_index < progress_markers.size()) {
            if (!progress_markers[progress_index]) {
                #pragma omp critical
                {
                    if (!progress_markers[progress_index]) {
                        std::cout << "|" << progress << std::flush;
                        progress_markers[progress_index] = true;
                    }
                }
            }
        }
    }
    auto end_time_ = std::chrono::high_resolution_clock::now();
    auto duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_).count();
    std::cout << "|  Time: " << duration_ / 1000. << " s" << std::endl;
}


void Graph_Wavelets::processNeighborCells(Cell& current_cell, std::set<size_t>& new_adjacent_cells, int8_t num_levels) {
    for (const auto& neighbor : current_cell.adjacent_cells) {
        Cell& neighbor_cell = cells[neighbor];

        if (neighbor_cell.level == num_levels) {
            // If the neighbour is at the same level, it is added
            new_adjacent_cells.insert(neighbor);
        } else {
            new_adjacent_cells.insert(neighbor_cell.merge);
        }
    }
}



/*******************************************************************
*				         WAVELET TRANSFORM 		        		   *
*******************************************************************/


void Graph_Wavelet_Transform(const ArgConfig& ARG_CONFIG, int8_t num_levels) {

    std::cout << "WAVELET TRANSFORM" << std::endl;

    std::vector<double> volume;
    std::vector<std::size_t> level;
    std::vector<std::size_t> merge;

    std::vector<std::vector<double>> signals;
    std::vector<std::string> signal_extensions;

    bool SUCCESS;
    SUCCESS = read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "level", level);
    SUCCESS *= read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "merge", merge);
    if (!SUCCESS) return;


    read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "vol", volume);
    if (ARG_CONFIG.BOOL_DIVERGENCE){
        signals.push_back(std::vector<double>());
        read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "div", signals.back());
        signal_extensions.push_back("." + ARG_CONFIG.OUT_SUFFIXES + "wdiv");
    }
    if (ARG_CONFIG.BOOL_CURL){
        std::vector<std::string> axes{"curl_x", "curl_y", "curl_z"};
        for (const auto& axis : axes) {
            signals.push_back(std::vector<double>());
            read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + axis, signals.back());
            signal_extensions.push_back("." + ARG_CONFIG.OUT_SUFFIXES + "w" + std::string(axis));
        }
    }


    int Np = static_cast<int>(volume.size());
    
    if (num_levels == -1){
        num_levels = static_cast<int>(std::floor(std::log2(Np)));
    }
    std::cout << "NUMBER OF LEVELS: " << static_cast<int>(num_levels) << std::endl;


    for (size_t current_level = 0; current_level < num_levels - 1; current_level++) {
        #pragma omp parallel for
        for(size_t odd_indices = 0; odd_indices < level.size(); odd_indices++){
            if(level[odd_indices] == current_level){
                size_t even_indices = merge[odd_indices];
                double volume_sum = volume[even_indices] + volume[odd_indices]; // volume conservation
                
                for(auto& signal : signals) {
                    double signal_projection = (volume[even_indices] * signal[even_indices] + volume[odd_indices] * signal[odd_indices]) / volume_sum; // signal conservation
                    double signal_prediction = signal_projection;
                    double signal_detail = signal[odd_indices] - signal_prediction;
                    // double signal_detail = (volume[even_indices] * (signal[even_indices] - signal_prediction) + volume[odd_indices] * (signal[odd_indices] - signal_prediction)) / volume_sum; // for symmetry // wrong !
                    signal[even_indices] = signal_projection;
                    signal[odd_indices] = signal_detail;
                    // signal[even_indices] += signal[odd_indices] / 2; // lifting scheme !!!!!!! 
                }

                volume[even_indices] = volume_sum;
            }
        }
    }

    write_double_array_to_binary_file(volume, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "wvol_odd");
    for (size_t s = 0; s < signals.size(); ++s) {
        write_double_array_to_binary_file(signals[s], ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + signal_extensions[s]);
    }
}


void Graph_Inverse_Wavelet_Transform(const ArgConfig& ARG_CONFIG, int8_t num_levels) {

    std::cout << "INVERSE WAVELET TRANSFORM" << std::endl;

    std::vector<double> volume;
    std::vector<std::size_t> level;
    std::vector<std::size_t> merge;
    
    std::vector<std::vector<double>> signals;
    std::vector<std::string> signal_extensions;

    bool SUCCESS;
    SUCCESS = read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "level", level);
    SUCCESS *= read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "merge", merge);
    if (!SUCCESS) return;

    // Read transformed signals from file
    read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "wvol_odd", volume);
    std::vector<double> volume_even(volume.size(), 0.0);

    int Np = static_cast<int>(volume.size());
    if (num_levels == -1){
        num_levels = static_cast<int>(std::floor(std::log2(Np)));
    }
    std::cout << "NUMBER OF LEVELS: " << static_cast<int>(num_levels) << std::endl;

    for (int current_level = num_levels - 2; current_level >= 0; current_level--) {
        #pragma omp parallel for
        for(size_t odd_indices = 0; odd_indices < volume.size(); odd_indices++){
            if(level[odd_indices] == current_level){
                size_t even_indices = merge[odd_indices];
                double volume_sum = volume[even_indices];
                volume[even_indices] = volume_sum - volume[odd_indices];
                volume_even[odd_indices] = volume[even_indices];
                
                for(auto& signal : signals) {
                    double signal_detail = signal[odd_indices];
                    double signal_prediction = signal[even_indices];
                    signal[odd_indices] = signal_detail + signal_prediction;
                    signal[even_indices] = (volume_sum * signal_prediction - volume[odd_indices] * signal[odd_indices]) / volume[even_indices];
                }
            }
        }
    }

    // Write original signals to file
    write_double_array_to_binary_file(volume_even, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "wvol_even");
    for (size_t s = 0; s < signals.size(); ++s) {
        write_double_array_to_binary_file(signals[s], ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + signal_extensions[s]);
    }
}


void Power_Spectrum(const ArgConfig& ARG_CONFIG, int num_levels) {

    std::cout << "POWER SPECTRUM" << std::endl;

    std::vector<double> volume_odd;
    std::vector<double> volume_even;
    std::vector<std::size_t> level;

    std::vector<std::vector<double>> signals;
    std::vector<std::string> signal_extensions;

    bool SUCCESS;
    SUCCESS = read_int_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "level", level);
    if (!SUCCESS) return;

    // Read transformed signals from file
    read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "wvol_odd", volume_odd);
    read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "wvol_even", volume_even);
    if (ARG_CONFIG.BOOL_DIVERGENCE){
        signals.push_back(std::vector<double>());
        read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "wdiv", signals.back());
        signal_extensions.push_back("." + ARG_CONFIG.OUT_SUFFIXES + "wps_div");
    }
    if (ARG_CONFIG.BOOL_CURL){
        for (const auto& axis : {"wcurl_x", "wcurl_y", "wcurl_z"}) {
            signals.push_back(std::vector<double>());
            read_double_array_from_binary_file(ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + axis, signals.back());
            signal_extensions.push_back("." + ARG_CONFIG.OUT_SUFFIXES + "wps_" + std::string(axis));
        }
    }

    int Np = static_cast<int>(volume_odd.size());
    if (num_levels == -1){
        num_levels = static_cast<int>(std::floor(std::log2(Np)));
    }


    double dim = 3.0;
    std::vector<double> mean_volume_odd(num_levels-1, 0.0);
    std::vector<double> mean_volume_even(num_levels-1, 0.0);
    std::vector<double> mean_volume(num_levels-1, 0.0);
    std::vector<double> volume_scale(num_levels-1, 0.0);
    std::vector<double> scalefactor(volume_odd.size(), 0.0);

    for (int lev = 1; lev < num_levels; lev++) {
        int count_odd = 0, count_even = 0;
        double sum_odd = 0.0, sum_even = 0.0;
        for (size_t i = 0; i < volume_odd.size(); i++) {
            if (level[i] == lev - 1) {
                sum_odd += volume_odd[i];
                count_odd++;
                sum_even += volume_even[i];
                count_even++;
            }
        }
        mean_volume_odd[lev - 1] = sum_odd / count_odd;
        mean_volume_even[lev - 1] = sum_even / count_even;
        mean_volume[lev - 1] = (mean_volume_odd[lev - 1] + mean_volume_even[lev - 1]) / 2;
        volume_scale[lev - 1] = std::pow(mean_volume[lev - 1], 1.0 / dim);
    }

    for (size_t i = 0; i < volume_odd.size(); i++) {
        scalefactor[i] = std::sqrt((volume_odd[i] + volume_even[i]) * volume_odd[i] / volume_even[i]); // !!!
        // scalefactor[i] = std::pow(2, level[i]/2) * std::pow(1+std::pow(volume_odd[i] / volume_even[i], 2.0), 1.0/2.0);
        // scalefactor[i] = 1.0; // no scale factor
        // scalefactor[i] = std::pow((volume_odd[i] + volume_even[i]) * volume_odd[i] / volume_even[i], 1.0/3.0);
    }
    
    std::vector<double> kv(num_levels-1, 0.0); 
    std::transform(volume_scale.begin(), volume_scale.end(), kv.begin(), [](double vs) { return PI / vs; });

    // Calculate statistics for each signal
    for (size_t s = 0; s < signals.size(); s++) {
        auto& signal = signals[s];
        std::vector<double> M1_L1(num_levels-1, 0.0);
        std::vector<double> M2_L1(num_levels-1, 0.0);
        std::vector<double> M1_L2(num_levels-1, 0.0);
        std::vector<double> M2_L2(num_levels-1, 0.0);
        std::vector<int> Nw(num_levels-1, 0);

        for (int lev = 1; lev < num_levels; lev++) {
            double sum_L1 = 0.0, sum_L2 = 0.0;
            double sum_square_L1 = 0.0, sum_square_L2 = 0.0;
            int count = 0;
            for (size_t i = 0; i < signal.size(); i++) {
                if (level[i] == lev - 1) {
                    double current_signal = signal[i]; 
                    sum_L1 += current_signal;
                    sum_L2 += scalefactor[i] * current_signal;
                    sum_square_L1 += std::pow(current_signal, 2);
                    sum_square_L2 += std::pow(scalefactor[i] * current_signal, 2);
                    count++;
                }
            }
            M1_L1[lev - 1] = sum_L1 / count;
            M2_L1[lev - 1] = sum_square_L1 / count;
            M1_L2[lev - 1] = sum_L2 / count;
            M2_L2[lev - 1] = sum_square_L2 / count;
        }

        std::vector<double> variance_L1(num_levels-1);
        std::transform(M2_L1.begin(), M2_L1.end(), M1_L1.begin(), variance_L1.begin(), [](double m2, double m1) { return m2 - std::pow(m1, 2); });

        std::vector<double> variance_L2(num_levels-1);
        std::transform(M2_L2.begin(), M2_L2.end(), M1_L2.begin(), variance_L2.begin(), [](double m2, double m1) { return m2 - std::pow(m1, 2); });

        // Compute Nw
        for (int lev = 1; lev < num_levels; lev++) {
            Nw[lev - 1] = std::count(level.begin(), level.end(), lev - 1);
        }

        // Compute E2_L2
        std::vector<double> dkv(num_levels-1);
        std::vector<double> E2_L2(num_levels-1);
        for (int lev = 1; lev < num_levels; lev++) {
            dkv[lev - 1] = (std::pow(2, 1 / dim) - 1) * kv[lev - 1] * std::log(2) / dim;
            E2_L2[lev - 1] = Nw[lev - 1] * M2_L2[lev - 1] / dkv[lev - 1] / std::pow(TWO_PI, dim);
        }
        write_double_array_to_binary_file(E2_L2, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + signal_extensions[s]);
    }
    write_double_array_to_binary_file(kv, ARG_CONFIG.DATA_PATH.second + "/" + ARG_CONFIG.FILENAME + "." + ARG_CONFIG.OUT_SUFFIXES + "kv");
}



/*******************************************************************
*				                 PRINT   		        		   *
*******************************************************************/


void Graph_Wavelets::printCell(size_t N) const {
    // const Cell& cell = cells[N];
    auto cell_iter = cells.find(N);
    if(cell_iter == cells.end()){
        std::cout << "Cell " << N << " does not exist in the map." << std::endl;
        return;
    }

    const Cell& cell = cell_iter->second;
    std::cout << "Cell " << N << ":" << std::endl;
    std::cout << "  Volume: " << cell.volume << std::endl;
    std::cout << "  Level: " << static_cast<int>(cell.level) << std::endl;
    std::cout << "  Merge: " << cell.merge << std::endl;
    std::cout << "  Domain: " << GlobalIndices_GlobalProcess[N] << std::endl;
    std::cout << "  Done: " << (cell.done ? "true" : "false") << std::endl;
    std::cout << "  Adjacent cells: ";
    for (const auto& adj : cell.adjacent_cells) {
        std::cout << adj << ", ";
    }
    std::cout << std::endl;
}


void Graph_Wavelets::printCell_bis(size_t N) const {
    // const Cell& cell = cells[N];
    auto cell_iter = cells.find(N);
    if(cell_iter == cells.end()){
        std::cout << "Cell " << N << " does not exist in the map." << std::endl;
        return;
    }

    const Cell& cell = cell_iter->second;
    std::cout << "Cell " << N << ":" << std::endl;
    std::cout << "  Merge: " << cell.merge << std::endl;
    std::cout << "  Adjacent cells: ";
    for (const auto& adj : cell.adjacent_cells) {
        std::cout << adj << ", ";
    }
    std::cout << std::endl;
}


void Graph_Wavelets::printFirstNCells(size_t N) const {
    for (size_t i = 0; i < N && i < cells.size(); ++i) {
        // const Cell& cell = cells[i];
        auto cell_iter = cells.find(i);
        if(cell_iter == cells.end()){
            std::cout << "Cell " << i << " does not exist in the map." << std::endl;
            continue;
        }

        const Cell& cell = cell_iter->second;
        std::cout << "Cell " << i << ":" << std::endl;
        std::cout << "  Volume: " << cell.volume << std::endl;
        std::cout << "  Level: " << static_cast<int>(cell.level) << std::endl;
        std::cout << "  Merge: " << cell.merge << std::endl;
        std::cout << "  Domain: " << GlobalIndices_GlobalProcess[N] << std::endl;
        std::cout << "  Done: " << (cell.done ? "true" : "false") << std::endl;
        std::cout << "  Adjacent cells: ";
        for (const auto& adj : cell.adjacent_cells) {
            std::cout << static_cast<int>(adj) << ", ";
        }
        std::cout << std::endl;
    }
}


void Graph_Wavelets::printAdjacentCellDetails(size_t index) const {

    // const Cell& cell = cells[index];
    auto cell_iter = cells.find(index);
    if(cell_iter == cells.end()){
        std::cout << "Cell " << index << " does not exist in the map." << std::endl;
        return;
    }

    const Cell& cell = cell_iter->second;

    std::cout << "Details for cell at index " << index << ":\n"
              << "Volume: " << cell.volume << "\n"
              << "Adjacent Cells:\n";

    for (const auto& neighbor : cell.adjacent_cells) {
        size_t neighbor_index = neighbor;
        auto cell_iter_neighbor = cells.find(neighbor_index);
        if(cell_iter == cells.end()){
            std::cout << "Cell " << neighbor_index << " does not exist in the map." << std::endl;
            continue;
        }

        const Cell& cell_neighbor = cell_iter->second;
        if (neighbor_index < cells.size()) {
            std::cout << "Neighbor index: " << neighbor_index 
                      << ", Neighbor volume: " << cell_neighbor.volume
                      << ", Neighbor done: " << cell_neighbor.done << "\n";
        } else {
            std::cout << "Neighbor index: " << neighbor_index << " is out of range.\n";
        }
    }
}



/*******************************************************************
*				                 WRITE   		        		   *
*******************************************************************/


void Graph_Wavelets::writeToBinary(const std::string& directory, const std::string& filename, const std::string& OUT_SUFFIXES, int number) {

    std::cout << "WRITE DATA : " << directory + "/GW_" + filename << std::endl;
    
    std::string subdomainNumber = "";
    if (number != -1) {
        subdomainNumber = "_" + std::to_string(number);
    }

    // Check that the directory exists, otherwise it will be created.
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directory(directory);
    }

    // Create a sub-directory to store the files
    std::string subdirectory = directory + "/GW_" + filename;
    if (!std::filesystem::exists(subdirectory)) {
        std::filesystem::create_directory(subdirectory);
    }

    std::ofstream key_file(subdirectory + "/" + filename + subdomainNumber + "." + OUT_SUFFIXES + "key", std::ios::binary);
    std::ofstream vol_file(subdirectory + "/" + filename + subdomainNumber + "." + OUT_SUFFIXES + "vol", std::ios::binary);
    std::ofstream level_file(subdirectory + "/" + filename + subdomainNumber + "." + OUT_SUFFIXES + "level", std::ios::binary);
    std::ofstream merge_file(subdirectory + "/" + filename + subdomainNumber + "." + OUT_SUFFIXES + "merge", std::ios::binary);
    std::ofstream neigh_val_file(subdirectory + "/" + filename + subdomainNumber + "." + OUT_SUFFIXES + "neigh_val", std::ios::binary);
    std::ofstream neigh_ptr_file(subdirectory + "/" + filename + subdomainNumber + "." + OUT_SUFFIXES + "neigh_ptr", std::ios::binary);

    // Extraire toutes les clés dans un vecteur
    std::vector<size_t> keys;
    for (const auto& cell : cells) {
        keys.push_back(cell.first);
    }

    // Trier le vecteur
    std::sort(keys.begin(), keys.end());

    // Itérer sur les clés dans l'ordre
    size_t ptr = 0;
    for (const auto& key : keys) {
        // Si la clé n'existe pas dans cells, continue
        if(cells.find(key) == cells.end())
            continue;

        const auto& cell = cells[key];

        key_file.write(reinterpret_cast<const char*>(&key), sizeof(size_t));
        
        vol_file.write(reinterpret_cast<const char*>(&cell.volume), sizeof(double));

        uint64_t level_as_uint64 = static_cast<uint64_t>(cell.level);
        level_file.write(reinterpret_cast<const char*>(&level_as_uint64), sizeof(uint64_t));

        merge_file.write(reinterpret_cast<const char*>(&cell.merge), sizeof(int64_t));

        size_t num_neigh = cell.adjacent_cells.size();
        neigh_ptr_file.write(reinterpret_cast<const char*>(&ptr), sizeof(size_t));

        for (const auto& neigh : cell.adjacent_cells) {
            neigh_val_file.write(reinterpret_cast<const char*>(&neigh), sizeof(size_t));
        }
        ptr += num_neigh;
    }
    // Write the final pointer value which also represents the total number of non-zero elements
    neigh_ptr_file.write(reinterpret_cast<const char*>(&ptr), sizeof(size_t));


    key_file.close();
    vol_file.close();
    level_file.close();
    merge_file.close();
    neigh_val_file.close();
    neigh_ptr_file.close();
}


void Graph_Wavelets::loadFromBinary(const std::string& directory, const std::string& filename, const std::string& OUT_SUFFIXES, size_t Np, size_t Np_Domain, int number) {

    std::cout << "READ DATA : " << directory + "/GW_" + filename << std::endl;

    int MPI_rank, MPI_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_size);

    std::string subdomainNumber = "";
    if (number != -1) {
        subdomainNumber = "_" + std::to_string(number);
    }

    // Ouvrir les fichiers binaires pour la lecture
    std::ifstream key_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainNumber + "." + OUT_SUFFIXES + "key", std::ios::binary);
    std::ifstream vol_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainNumber + "." + OUT_SUFFIXES + "vol", std::ios::binary);
    std::ifstream level_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainNumber + "." + OUT_SUFFIXES + "level", std::ios::binary);
    std::ifstream merge_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainNumber + "." + OUT_SUFFIXES + "merge", std::ios::binary);
    std::ifstream neigh_val_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainNumber + "." + OUT_SUFFIXES + "neigh_val", std::ios::binary);
    std::ifstream neigh_ptr_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainNumber + "." + OUT_SUFFIXES + "neigh_ptr", std::ios::binary);

    // Vérifier que les fichiers ont été ouverts correctement
    if (!key_file) {
        std::cerr << "Error opening key_file for reading." << std::endl;
        return;
    }
    if (!vol_file) {
        std::cerr << "Error opening vol_file for reading." << std::endl;
        return;
    }
    if (!level_file) {
        std::cerr << "Error opening level_file for reading." << std::endl;
        return;
    }
    if (!merge_file) {
        std::cerr << "Error opening merge_file for reading." << std::endl;
        return;
    }
    if (!neigh_val_file) {
        std::cerr << "Error opening neigh_val_file for reading." << std::endl;
        return;
    }
    if (!neigh_ptr_file) {
        std::cerr << "Error opening neigh_ptr_file for reading." << std::endl;
        return;
    }

    size_t key;
    double volume;
    uint64_t level_as_uint64;
    int64_t merge;
    size_t ptr, prev_ptr = 0;
    std::vector<size_t> neigh_val;
    neigh_ptr_file.read(reinterpret_cast<char*>(&ptr), sizeof(size_t)); // read the first : 0
    while(key_file.read(reinterpret_cast<char*>(&key), sizeof(size_t))) {

        // Lire les valeurs des autres fichiers
        vol_file.read(reinterpret_cast<char*>(&volume), sizeof(double));
        level_file.read(reinterpret_cast<char*>(&level_as_uint64), sizeof(uint64_t));
        merge_file.read(reinterpret_cast<char*>(&merge), sizeof(int64_t));
        neigh_ptr_file.read(reinterpret_cast<char*>(&ptr), sizeof(size_t));

        int8_t level = static_cast<int8_t>(level_as_uint64);

        // Lire les voisins
        neigh_val.resize(ptr - prev_ptr);
        for(size_t& neigh : neigh_val) {
            neigh_val_file.read(reinterpret_cast<char*>(&neigh), sizeof(size_t));
        }

        // Insérer les valeurs dans la cellule
        auto& cell = cells[key];
        cell.volume = volume;
        cell.level = level;
        cell.merge = merge;
        cell.done = false;
        cell.adjacent_cells.insert(neigh_val.begin(), neigh_val.end());


        prev_ptr = ptr;
    }

    // Fermer les fichiers
    key_file.close();
    vol_file.close();
    level_file.close();
    merge_file.close();
    neigh_val_file.close();
    neigh_ptr_file.close();


    // Initialisation de cells_indices et previous_level_indices
    for (const auto& cell : cells) {
        cells_indices.insert(cell.first);
        previous_level_indices.push_back(cell.first);
    }

    // Initialize GlobalIndices_GlobalProcess
    GlobalIndices_GlobalProcess.resize(Np, -1);  // Init with -1 assuming no process yet assigned
    
    
    for(int subdomainNumber = 0; subdomainNumber < MPI_size; ++subdomainNumber) {
        std::string subdomainStr = "_" + std::to_string(subdomainNumber);
        std::ifstream key_file(directory + "/GW_" + filename + "_LVL0/" + filename + "_LVL0" + subdomainStr + "." + OUT_SUFFIXES + "key", std::ios::binary);

        size_t key;
        while(key_file.read(reinterpret_cast<char*>(&key), sizeof(size_t))) {
            GlobalIndices_GlobalProcess[key] = subdomainNumber;
        }

        key_file.close();
    }
    
    std::cout << "SIZE CELLS: " << cells.size() << "  ||  MAIN CELLS: " << cells.size() << std::endl;
}





















