/*
  A pseudo-representative application to model NEK.

Modified:
    Nalini Kumar { UF CCMT }


Original:
    Copyright (C) 2016  { Dylan Rudolph, NSF CHREC, UF CCMT }

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARAMS_H_
#define PARAMS_H_


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>


struct paramstype {

/* -------------------------- Machine/Primary Parameters --------------------------- */
  unsigned int PROBED_RANK; 	// The rank which shows its timing output
  unsigned int CARTESIAN_X, CARTESIAN_Y, CARTESIAN_Z;	// Number of processes in each dimension. MPI must be run with (product of these numbers) processes.

/* -------------------------- Physics/Application Parameters --------------------------- */
  unsigned int TIMESTEPS;		// Number of simulation timesteps
  unsigned int RK;			// Number of Runge Kutta stages
  unsigned int ELEMENT_SIZE;		// Size of each element (cubic)
  unsigned int PHYSICAL_PARAMS;		// Number of physics parameters tracked
  unsigned int ELEMENTS_X, ELEMENTS_Y, ELEMENTS_Z;	// Number of elements per process in each dimension
  unsigned int ELEMENTS_PER_PROCESS;
  unsigned int ELEMENTS_ON_X_FACE, ELEMENTS_ON_Y_FACE, ELEMENTS_ON_Z_FACE;
  unsigned int FACE_SIZE;
  
};


/* Set machine & application parameters from user-specified command line arguments*/
void setup_parameters(int argc, char *argv[], int rank, struct paramstype *params);

void print_parameters(struct paramstype *params);

#endif

