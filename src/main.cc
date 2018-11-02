#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>

double alpha = 1; // Heat Constant
double r = 0.25; // r = k/(h^2); k = 1; h = 2;

int main(int argc, char **argv) {

    double t1_temp = atof(argv[1]);
    double t2_temp = atof(argv[2]);
    size_t num_grid_points = atof(argv[3]);
    int num_timesteps = atof(argv[4]);
    
    MPI_Init(&argc, &argv);

    int MPI_num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_num_ranks);
    int MPI_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
    int grid_size;
    int grid_prev_size;

    std::vector<double> grid;
    std::vector<double> grid_prev;
    
    // ignore ranks if too many
    if (num_grid_points < MPI_num_ranks) {
        MPI_num_ranks = num_grid_points;
    }

    if (MPI_rank < MPI_num_ranks) {
    // Divide up grip among the ranks last rank gets remaining grids
    // edge case check
        if (MPI_num_ranks > 1) {
            if (MPI_rank != MPI_num_ranks - 1) {
                grid_size = int(num_grid_points / MPI_num_ranks);
            } else {
                grid_size = num_grid_points - int(num_grid_points / MPI_num_ranks) * (MPI_num_ranks - 1);
            }

            // Assign vector used to calculate values edges get one less space since they use t1 and t2
            if (MPI_rank == 0 || MPI_rank == MPI_num_ranks - 1) {
                grid_prev_size = grid_size + 1;
            } else {
                grid_prev_size = grid_size + 2;
            }

            grid = std::vector<double>(grid_size, 0.0);
            grid_prev = std::vector<double>(grid_prev_size, 0.0);

            for (int i = 0; i < num_timesteps; ++i) {
                // Calculate values for each timestep
                for (int k = 0; k < grid.size(); ++k) {
                    if (MPI_rank == 0) {
                        if (k == 0) {
                            // point touching t1_temp
                            grid[k] = (1 - 2 * r) * grid_prev[k] + r * t1_temp + r * grid_prev[k + 1];
                        } else {
                            grid[k] = (1 - 2 * r) * grid_prev[k] + r * grid_prev[k - 1] + r * grid_prev[k + 1];
                        }
                    } else if (MPI_rank == MPI_num_ranks - 1 && k == grid.size() - 1) {
                        // point touching t2_temp
                        //indices for grid and grid_prev are shifted
                        grid[k] = (1 - 2 * r) * grid_prev[k + 1] + r * grid_prev[k] + r * t2_temp; 
                    } else {
                        // indices for grid and grid_prev are shifted
                        grid[k] = (1 - 2 * r) * grid_prev[k + 1] + r * grid_prev[k] + r * grid_prev[k + 2];
                    }
                }

                // Copy over grid values to grid_prev
                if (MPI_rank == 0) {
                    for (int l = 0; l < grid.size(); ++l) {
                        grid_prev[l] = grid[l];
                    }
                } else {
                    for (int l = 0; l < grid.size(); ++l) {
                        grid_prev[l + 1] = grid[l];
                    }
                }
                // get remain grid_prev values
                if (MPI_rank % 2 == 0) {
                    if (MPI_rank == 0) {
                        MPI_Recv(&grid_prev[grid.size()], 1, MPI_DOUBLE, MPI_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(&grid[grid.size() - 1], 1, MPI_DOUBLE, MPI_rank + 1, 0, MPI_COMM_WORLD);
                    } else if (MPI_rank == MPI_num_ranks - 1) {
                        MPI_Recv(&grid_prev[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(&grid[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD);
                    } else {
                        MPI_Recv(&grid_prev[grid.size() + 1], 1, MPI_DOUBLE, MPI_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&grid_prev[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(&grid[grid.size() - 1], 1, MPI_DOUBLE, MPI_rank + 1, 0, MPI_COMM_WORLD);
                        MPI_Send(&grid[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD); 
                    }
                } else {
                    if (MPI_rank == MPI_num_ranks - 1) {
                        MPI_Send(&grid[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&grid_prev[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    } else {
                        MPI_Send(&grid[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD); 
                        MPI_Send(&grid[grid.size() - 1], 1, MPI_DOUBLE, MPI_rank + 1, 0, MPI_COMM_WORLD);
                        MPI_Recv(&grid_prev[0], 1, MPI_DOUBLE, MPI_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&grid_prev[grid.size() + 1], 1, MPI_DOUBLE, MPI_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
        } else {
            // Only one process is called or needed
            grid = std::vector<double>(num_grid_points, 0.0);
            grid_size = grid.size();

            if (num_grid_points > 1) {
                grid_prev = std::vector<double>(num_grid_points, 0.0);
                for (int i = 0; i < num_timesteps; ++i) {
                    for (int k = 0; k < num_grid_points; ++k) {
                        if (k == 0) {
                            grid[k] = (1 - 2 * r) * grid_prev[k] + r * t1_temp + r * grid_prev[k + 1];
                        } else if (k == num_grid_points - 1) {
                            grid[k] = (1 - 2 * r) * grid_prev[k] + r * grid_prev[k - 1] + r * t2_temp;
                        } else {
                            grid[k] = (1 - 2 * r) * grid_prev[k] + r * grid_prev[k - 1] + r * grid_prev[k + 1];
                        }
                    }
                    grid_prev = grid;
                }
            } else {
                for (int i = 0; i < num_timesteps; ++i) {
                    // TODO steady state after once change to if statement check
                    grid[0] = (1 - 2 * r) * grid[0] + r * t1_temp + r * t2_temp;
                }
            }       
        }
        // merge all values into rank 0
        std::vector<double> outputs(num_grid_points);
        if (MPI_rank == 0) {
            // Copy rank 0 grid into outputs
            for (int point_index = 0; point_index < grid.size(); ++point_index) {
                outputs[point_index] = grid[point_index];
            }

            if (MPI_num_ranks > 1) {
                // Collect remaining grid points from other ranks
                int recv_size = grid.size();
                int outputs_index = recv_size;
                for (int m = 1; m < MPI_num_ranks - 1; ++m) {
                    MPI_Recv(&outputs[outputs_index], recv_size, MPI_DOUBLE, m, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    outputs_index += recv_size;
                }
                // last rank might differ in size
                MPI_Recv(&outputs[outputs_index], num_grid_points - recv_size * (MPI_num_ranks - 1), MPI_DOUBLE, MPI_num_ranks - 1,
                    0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }


            // output to file
            std::ofstream outfile;
            outfile.open("heat1Doutput.csv");
            for (int n = 0; n < outputs.size() - 1; ++n) {
                outfile << outputs[n] << ", ";
                // would add a new line character after a certain amount to aid in open file
            }
            outfile << outputs[outputs.size() - 1];
            outfile.close();
        } else {
            MPI_Send(&grid[0], grid.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    return 0;
}
