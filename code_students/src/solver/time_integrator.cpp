#include "solver/time_integrator.hpp"
#include "setup/grid.hpp"

#include <iostream>

RungeKutta2::RungeKutta2(const grid_3D &grid, const fluid &fluid3D) : simulation_grid(grid), fluid_storage(fluid3D.get_fluid_type()) {
	fluid_storage.setup(simulation_grid);
}

void RungeKutta2::save_data(const fluid &fluid3D) {
	for (int i_field = 0; i_field < fluid3D.get_num_fields(); ++i_field) {
		fluid_storage.fluid_data[i_field] = fluid3D.fluid_data[i_field];
	}
}

void RungeKutta2::do_sub_step(const grid_3D &grid, const fluid &changes, fluid &fluid3D, double delta_t, int sub_step) {

	int num_fields = fluid3D.get_num_fields();

	std::cout << " Doing Runge-Kutta step " << sub_step << "\n";
	if (sub_step == 0) {
		save_data(fluid3D);
	}

	for (int i_field = 0; i_field < num_fields; ++i_field) {
		if (sub_step == 0) {

			// TBD by students -> do first Runge-Kutta step

		} else if (sub_step == 1) {

			// TBD by students -> do second Runge-Kutta step

		} else {
			std::cerr << " Error no such substep: " << sub_step << "\n";
			exit(3);
		}
	}
}

int RungeKutta2::get_number_substeps() const { return 2; }