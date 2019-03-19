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


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>

#include "params.h"


/* Set machine & application parameters from user-specified command line arguments*/
void setup_parameters(int argc, char *argv[], int rank, struct paramstype *params)
{
  params->PROBED_RANK = 0;	// The rank which shows its timing output
  params->RK = 3;		// Number of Runge Kutta stages

  params->TIMESTEPS = 3;		// Number of simulation timesteps
  params->ELEMENT_SIZE = 5;	// Size of each element (cubic)
  params->PHYSICAL_PARAMS = 5;	// Number of physics parameters tracked
  params->ELEMENTS_X=2; 
  params->ELEMENTS_Y=2; 
  params->ELEMENTS_Z=2;		// Number of elements perprocess
  params->CARTESIAN_X=2; 
  params->CARTESIAN_Y=2; 
  params->CARTESIAN_Z = 2;	

/*  if (rank == params->PROBED_RANK) {
    printf ("Command line arguments are processed in the following order.\nTIMESTEPS, ELEMENT_SIZE, ELEMENTS_X, ELEMENTS_Y, ELEMENTS_Z, CARTESIAN_X, CARTESIAN_Y, CARTESIAN_Z, PHYSICAL_PARAMS.\n\n");
    printf ("Input args = %d\n\n",argc);
  }
*/
  if (argc > 10) {
    if (rank == params->PROBED_RANK) { printf("Too many input arguments specified. Defaulting to preassigned values for all 9 parameters in params.h"); }
  }
 
  if ( argc == 10) {
    params->TIMESTEPS = atoi( argv[1] );
    params->ELEMENT_SIZE = atoi( argv[2] );
    params->ELEMENTS_X = atoi( argv[3] );
    params->ELEMENTS_Y = atoi( argv[4] );
    params->ELEMENTS_Z = atoi( argv[5] );
    params->CARTESIAN_X = atoi( argv[6] );
    params->CARTESIAN_Y = atoi( argv[7] );
    params->CARTESIAN_Z = atoi( argv[8] );
    params->PHYSICAL_PARAMS= atoi( argv[9] );
//    if (rank == params->PROBED_RANK) { printf("Using all 9 user specified parameters. \n"); }
  }

  else if ( argc == 9) {
    params->TIMESTEPS = atoi( argv[1] );
    params->ELEMENT_SIZE = atoi( argv[2] );
    params->ELEMENTS_X = atoi( argv[3] );
    params->ELEMENTS_Y = atoi( argv[4] );
    params->ELEMENTS_Z = atoi( argv[5] );
    params->CARTESIAN_X = atoi( argv[6] );
    params->CARTESIAN_Y = atoi( argv[7] );
    params->CARTESIAN_Z = atoi( argv[8] );
//    if (rank == params->PROBED_RANK) { printf("8 parameters specified. Using default values for others from params.h. \n"); }
  }

  else if ( argc == 6) {
    params->TIMESTEPS = atoi( argv[1] );
    params->ELEMENT_SIZE = atoi( argv[2] );
    params->ELEMENTS_X = atoi( argv[3] );
    params->ELEMENTS_Y = atoi( argv[4] );
    params->ELEMENTS_Z = atoi( argv[5] );
//    if (rank == params->PROBED_RANK) { printf("5 parameters specified. Using default values for others from params.h. \n"); }
  }

  else if ( argc == 3) {
    params->TIMESTEPS = atoi( argv[1] );
    params->ELEMENT_SIZE = atoi( argv[2] );
//    if (rank == params->PROBED_RANK) { printf("2 parameters specified. Using default values for others from params.h. \n"); }
  }

  else if ( argc == 2) {
    params->TIMESTEPS = atoi( argv[1] );
//    if (rank == params->PROBED_RANK) { printf("1 parameter specified. Using default values for others from params.h. \n"); }
  }
  
  else if ( argc == 1) {
    if (rank == params->PROBED_RANK) { printf("No input arguments specified by user. Using defaults from params.h. \n"); }
  }

  else {
    if (rank == params->PROBED_RANK) { printf("Invalid number of arguments. Command line parameter values will not be used. Using defaults from params.h. \n"); }
  }

  params->ELEMENTS_PER_PROCESS = params->ELEMENTS_X * params->ELEMENTS_Y * params->ELEMENTS_Z;
  params->ELEMENTS_ON_X_FACE = params->ELEMENTS_Y * params->ELEMENTS_Z;
  params->ELEMENTS_ON_Y_FACE = params->ELEMENTS_X * params->ELEMENTS_Z;
  params->ELEMENTS_ON_Z_FACE = params->ELEMENTS_X * params->ELEMENTS_Y;
  params->FACE_SIZE = params->ELEMENT_SIZE * params->ELEMENT_SIZE; 

}


void print_parameters(struct paramstype *params) {
/*
    printf ( "\nTIMESTEPS = %d", params->TIMESTEPS );
    printf ( " ELEMENT_SIZE = %d", params->ELEMENT_SIZE );
    printf ( " ELEMENTS_X = %d", params->ELEMENTS_X );
    printf ( " ELEMENTS_Y = %d", params->ELEMENTS_Y );
    printf ( " ELEMENTS_Z = %d", params->ELEMENTS_Z );
    printf ( " CARTESIAN_X = %d", params->CARTESIAN_X );
    printf ( " CARTESIAN_Y = %d", params->CARTESIAN_Y );
    printf ( " CARTESIAN_Z = %d", params->CARTESIAN_Z );
    printf ( " PHYSICAL_PARAMS = %d", params->PHYSICAL_PARAMS );
*/
    printf ( "%d,%d,%d,%d,%d,%d,%d,%d,%d\n", params->TIMESTEPS, params->ELEMENT_SIZE, params->ELEMENTS_X, params->ELEMENTS_Y, params->ELEMENTS_Z, params->CARTESIAN_X, params->CARTESIAN_Y, params->CARTESIAN_Z, params->PHYSICAL_PARAMS );
}


