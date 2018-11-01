#include <iostream>
#include <vector>

double alpha = 1; // Heat Constant
double r = 0.25; // r = k/(h^2); k = 1; h = 2;
int main(int arc, char *argv[]) {
    double t1_temp = atof(argv[1]);
    double t2_temp = atof(argv[2]);
    size_t num_grid_points = atof(argv[3]);
    int num_timesteps = atof(argv[4]);
    std::vector<double> grid_prev(num_grid_points, 0.0);
    std::vector<double> grid(num_grid_points, 0.0);
    if (num_grid_points > 1) {
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
            grid[0] = (1 - 2 * r) * grid_prev[0] + r * t1_temp + r * t2_temp;
        }
    }
    for (int l = 0; l < num_grid_points - 1; ++l) {
        std::cout << grid[l] << ", ";
    }
    std::cout << grid[num_grid_points - 1];
    return 0;
}
