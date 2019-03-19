/*
  A pseudo-representative application to model NEK.

Modified:
   Nalini Kumar  { UF CCMT }
   Aravind Neelakantan { UF CCMT }

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
#include <math.h>
#include <mpi.h>
#include <assert.h>

#include "params.h"
#include "dstructs.h" 
#include "utils.h"
#include "flux.h"



/* -------------------------- Machine/Primary Parameters --------------------------- */
  #define MPI_DTYPE MPI_DOUBLE		// MPI datatypes
  #define CARTESIAN_DIMENSIONS 3  	// Setup for MPI Cartesian function calls
  #define CARTESIAN_REORDER 0
  #define CARTESIAN_WRAP {0, 0, 0}

/* --------------------------- Timing Parameter selection  -------------------------- */
     /* Use this to print the timing parameter on specific module in the program
        PROFILE == FALSE  ---->   Print each timestep (csv format) and its avg.
        PROFILE == TRUE   ---->   Compute(A), comm and compute(B) per step with avg.  */

  #define PROFILE



/* ------------------------------ Main Loop ----------------------------------------- */

int main (int argc, char *argv[])
{

  /* ------------------------------- MPI Setup------------------------------ */

  MPI_Init(&argc, &argv);
    
  int rank, comrades;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comrades);


  /* ------------------------------- PARAMETER Setup------------------------------ */

  struct paramstype *params = malloc(sizeof(struct paramstype));
  assert(params != NULL);

  setup_parameters( argc, argv, rank, params);
  if (rank == params->PROBED_RANK) { print_parameters(params); }

  int cart_sizes[CARTESIAN_DIMENSIONS] = {params->CARTESIAN_X, params->CARTESIAN_Y, params->CARTESIAN_Z};
  int cart_wrap[CARTESIAN_DIMENSIONS] = CARTESIAN_WRAP;

  MPI_Comm cart_comm;
  MPI_Cart_create( MPI_COMM_WORLD, CARTESIAN_DIMENSIONS, cart_sizes, cart_wrap,
                   CARTESIAN_REORDER, &cart_comm );



  /* ------------------------------ Timing Setup --------------------------- */
#ifndef PROFILE
  struct timespec tA, tB;
  double t_steps[params->TIMESTEPS], t_avg, t_sum = 0;
#endif

#ifdef PROFILE
  struct timespec tcompA_s, tcompA_e;
  struct timespec tcompB_s, tcompB_e;
  struct timespec tcomm_s, tcomm_e;
    
  int TSxRK = params->TIMESTEPS * (params->RK);
  
  double t_steps_compA[TSxRK], t_avg_compA, t_sum_compA = 0;
  double t_steps_compB[TSxRK], t_avg_compB, t_sum_compB = 0;
  double t_steps_comm[TSxRK], t_avg_comm, t_sum_comm = 0;

  int iA, trA = 0;
  int iB, trB = 0;
  int iC, trC = 0;
#endif


  /* ------------------------------ Memory Setup --------------------------- */
  srand( 11 );

  /* Index variables: {generic, timestep, params->RK-index, element, block, axis} */
  int i, t, r, e, b, axis;

  element elements_Q[ params->ELEMENTS_PER_PROCESS ];
  element elements_R[ params->ELEMENTS_PER_PROCESS ];

  for (e = 0; e < params->ELEMENTS_PER_PROCESS; e++) {
    elements_Q[e] = new_random_element(0, 10, params);
    elements_R[e] = new_zero_element(params);
  }

  /* The same kernel is used for everything */
  matrix kernel = new_random_matrix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, -10, 10);

  /* The same transformation ternix (RX) is used for all elements.
     This is an approximation, there should be one for each element. */
  ternix RX[9];

  for (i = 0; i < 9; i++) {
    RX[i] = new_random_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE, -1, 1);
  }

  /* Intermediate 3D structures: used in conv operation */
  ternix Hx = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  ternix Hy = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  ternix Hz = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);

  /* Intermediate 3D structures: outputs of conv operation */
  ternix Ur = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  ternix Us = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  ternix Ut = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);

  /* Intermediate 3D structures: outputs of derivative operations */
  ternix Vr = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  ternix Vs = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  ternix Vt = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);



  /* ----------------------------------------------------------------------- */
  /* ------------------------------- Main Loop ----------------------------- */
  /* ----------------------------------------------------------------------- */

  /* For each timestep: */
  for ( t = 0; t < params->TIMESTEPS; t++ ) {

#ifndef PROFILE
    if (rank == params->PROBED_RANK) tA = now();
#endif

    /* For each of the three 'stages': */
    for (r = 0; r < params->RK; r++) {


      /* --------------------------- Compute (A) --------------------------- */
#ifdef PROFILE
      if (rank == params->PROBED_RANK) { tcompA_s = now(); }
#endif
      /* For each element owned by this rank: */
      for ( e = 0; e < params->ELEMENTS_PER_PROCESS; e++ ) {

        /* For each block in the element: */
        for ( b = 0; b < params->PHYSICAL_PARAMS; b++ ) {

          /* Generate Ur, Us, and Ut. */
          operation_conv(elements_Q[e]->B[b], RX, Hx, Hy, Hz, Ur, Us, Ut, params);

          /* Perform the three derivative computations (R, S, T). */
          operation_dr(kernel, Ur, Vr, params);
          operation_ds(kernel, Us, Vs, params);
          operation_dt(kernel, Ut, Vt, params);

          /* Add Vr, Vs, and Vt to make R. */
          operation_sum( Vr, Vs, Vt, elements_R[e]->B[b], params );

        }
      }
#ifdef PROFILE  
      if (rank == params->PROBED_RANK) { 
		  tcompA_e = now();
		  t_steps_compA[trA] = tdiff(tcompA_s, tcompA_e);
		  t_sum_compA += t_steps_compA[trA];
		  trA = trA + 1;
      }
#endif


      /* --------------------------- Communicate --------------------------- */
      /* above: plus neighbor, below: minus neighbor, index along this axis */
      int above, below, index;

      /* Unused status flag */
      MPI_Status status;

      /* Cartesian coordinates */
      int coords[CARTESIAN_DIMENSIONS];

      /* Determine our location in the cartesian grid. */
      MPI_Cart_coords(cart_comm, rank, CARTESIAN_DIMENSIONS, coords);

      vector above_faces_to_send, above_faces_to_recv;
      vector below_faces_to_send, below_faces_to_recv;

#ifdef PROFILE
      if (rank == params->PROBED_RANK) { tcomm_s = now(); }
#endif
      for ( axis = 0; axis < CARTESIAN_DIMENSIONS; axis++ ) {

        /* Find our index along this axis. */
        index = coords[axis];

        /* Determine our neighbors. */
        MPI_Cart_shift(cart_comm, axis, 1, &below, &above);

        /* --------------------------- Transfers --------------------------- */

        /* Significant operations are given a heading, everything else is just
           instrumentation and logging. */

        /* ------------------------ Even Axis Index ------------------------ */

        if ( (index % 2) == 0 ) {

          /* If my index on this axis is even:
             - SEND  faces to    ABOVE  neighbor  (23)
             - RECV  faces from  ABOVE  neighbor  (47)
             - SEND  faces to    BELOW  neighbor  (61)
             - RECV  faces from  BELOW  neighbor  (73) */

          if ( above != MPI_PROC_NULL ) {

            /* - - - - - - - - - - - - Prepare Faces - - - - - - - - - - - - */
            above_faces_to_send = new_extracted_faces(elements_R, axis, 1, params);
            above_faces_to_recv = new_empty_faces(axis, params);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Send Above  - - - - - - - - - - - - */
            MPI_Send( above_faces_to_send->V, above_faces_to_send->size,
                      MPI_DTYPE, above, 23, cart_comm );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Recv Above  - - - - - - - - - - - - */
            MPI_Recv( above_faces_to_recv->V, above_faces_to_recv->size,
                      MPI_DTYPE, above, 47, cart_comm, &status );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - Cleanup Faces - - - - - - - - - - - - */
            delete_vector(above_faces_to_send);
            delete_vector(above_faces_to_recv);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

          }

          if ( below != MPI_PROC_NULL ) {

            /* - - - - - - - - - - - - Prepare Faces - - - - - - - - - - - - */
            below_faces_to_send = new_extracted_faces(elements_R, axis, -1, params);
            below_faces_to_recv = new_empty_faces(axis, params);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Send Below  - - - - - - - - - - - - */
            MPI_Send( below_faces_to_send->V, below_faces_to_send->size,
                      MPI_DTYPE, below, 61, cart_comm );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Recv Below  - - - - - - - - - - - - */
            MPI_Recv( below_faces_to_recv->V, below_faces_to_recv->size,
                      MPI_DTYPE, below, 73, cart_comm, &status );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - Cleanup Faces - - - - - - - - - - - - */
            delete_vector(below_faces_to_send);
            delete_vector(below_faces_to_recv);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

          }
        }

        /* ------------------------- Odd Axis Index ------------------------ */

        else {

          /* If my index on this axis is odd:
             - RECV  faces from  BELOW  neighbor  (23)
             - SEND  faces from  BELOW  neighbor  (47)
             - RECV  faces to    ABOVE  neighbor  (61)
             - SEND  faces from  ABOVE  neighbor  (73) */

          if ( below != MPI_PROC_NULL ) {

            /* - - - - - - - - - - - - Prepare Faces - - - - - - - - - - - - */
            below_faces_to_send = new_extracted_faces(elements_R, axis, -1, params);
            below_faces_to_recv = new_empty_faces(axis, params);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Recv Below  - - - - - - - - - - - - */
            MPI_Recv( below_faces_to_recv->V, below_faces_to_recv->size,
                      MPI_DTYPE, below, 23, cart_comm, &status );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Send Below  - - - - - - - - - - - - */
            MPI_Send( below_faces_to_send->V, below_faces_to_send->size,
                      MPI_DTYPE, below, 47, cart_comm );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - Cleanup Faces - - - - - - - - - - - - */
            delete_vector(below_faces_to_send);
            delete_vector(below_faces_to_recv);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

          }

          if ( above != MPI_PROC_NULL ) {

            /* - - - - - - - - - - - - Prepare Faces - - - - - - - - - - - - */
            above_faces_to_send = new_extracted_faces(elements_R, axis, 1, params);
            above_faces_to_recv = new_empty_faces(axis, params);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Recv Above  - - - - - - - - - - - - */
            MPI_Recv( above_faces_to_recv->V, above_faces_to_recv->size,
                      MPI_DTYPE, above, 61, cart_comm, &status );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - - Send Above  - - - - - - - - - - - - */
            MPI_Send( above_faces_to_send->V, above_faces_to_send->size,
                      MPI_DTYPE, above, 73, cart_comm );
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

            /* - - - - - - - - - - - - Cleanup Faces - - - - - - - - - - - - */
            delete_vector(above_faces_to_send);
            delete_vector(above_faces_to_recv);
            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

          }

        }

      } /* for each axis ... */

#ifdef PROFILE  
      if (rank == params->PROBED_RANK) { 
		  tcomm_e = now();
		  t_steps_comm[trC] = tdiff(tcomm_s, tcomm_e);
		  t_sum_comm += t_steps_comm[trC];
		  trC = trC + 1;
      }
#endif


      /* --------------------------- Compute (B) --------------------------- */
#ifdef PROFILE
      if (rank == params->PROBED_RANK) { tcompB_s = now(); }
#endif

      /* For each element owned by this rank: */
      for ( e = 0; e < params->ELEMENTS_PER_PROCESS; e++ ) {

        /* For each block in the element: */
        for ( b = 0; b < params->PHYSICAL_PARAMS; b++ ) {

          /* Perform a fake Runge Kutta stage (without R from the last stage)
             to obtain a new value of Q. */
          operation_rk(elements_R[e]->B[b], elements_Q[e]->B[b], params);

        }
      }
#ifdef PROFILE  
      if (rank == params->PROBED_RANK) { 
		  tcompB_e = now();
		  t_steps_compB[trB] = tdiff(tcompB_s, tcompB_e);
		  t_sum_compB += t_steps_compB[trB];
		  trB = trB + 1;
      }
#endif      
      
      
    } /* For each stage ... */

#ifndef PROFILE
    if (rank == params->PROBED_RANK) { 
	  tB = now();
      printf("%.8f,", tdiff(tA, tB));
      t_steps[t] = tdiff(tA, tB);
      t_sum += t_steps[t];
    }
#endif

  } /* for each timestep ... */


  /* ------- Print execution time profiling outputs ------------------------ */
#ifndef PROFILE
  if (rank == params->PROBED_RANK) {
    t_avg = t_sum/(params->TIMESTEPS);
    printf("\nAverage time: %.8f\n", t_avg);
  }
#endif

#ifdef PROFILE
  if ( rank == params->PROBED_RANK ) {

    /* -------- Compute A -------- */
    for (iA = 0; iA < TSxRK; iA++) { 
	  if (iA == (TSxRK-1)) {
        printf("%.8f\n", t_steps_compA[iA]);
      }
      else {	
		printf("%.8f,", t_steps_compA[iA]); 
	  }
	}

    /* -------- Compute B -------- */
    for (iB = 0; iB < TSxRK; iB++) { 
	  if (iB == (TSxRK-1)) {
        printf("%.8f\n", t_steps_compB[iB]);
      }
      else {	
		printf("%.8f,", t_steps_compB[iB]); 
	  }
	}
		
    /* -------- Communication-------- */
    for (iC = 0; iC < TSxRK; iC++) { 
	  if (iC == (TSxRK-1)) {
        printf("%.8f\n", t_steps_comm[iC]);
      }
      else {	
		printf("%.8f,", t_steps_comm[iC]); 
	  }
	}


    t_avg_compA = t_sum_compA/TSxRK;
    t_avg_compB = t_sum_compB/TSxRK;
    t_avg_comm = t_sum_comm/TSxRK;

    printf("Average: %.8f, %.8f, %.8f\n", t_avg_compA, t_avg_compB, t_avg_comm);

    //printf("Average time of compute(B): %.8f\n", t_avg_compB);

    //printf("Average time of communication: %.8f\n", t_avg_comm);

  }
#endif 


  /* ----------------------------------------------------------------------- */
  /* -------------------------------- Cleanup ------------------------------ */
  /* ----------------------------------------------------------------------- */

  for (e = 0; e < params->ELEMENTS_PER_PROCESS; e++) {
    delete_element(elements_Q[e], params);
    delete_element(elements_R[e], params);
  }

  delete_matrix(kernel);

  for (i = 0; i < 9; i++) {
    delete_ternix(RX[i]);
  }

  delete_ternix(Hx);
  delete_ternix(Hy);
  delete_ternix(Hz);
  delete_ternix(Ur);
  delete_ternix(Us);
  delete_ternix(Ut);
  delete_ternix(Vr);
  delete_ternix(Vs);
  delete_ternix(Vt);

  free(params);
  
  MPI_Finalize();

  return 0;
}

