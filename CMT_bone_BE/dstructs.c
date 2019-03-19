/*
  A pseudo-representative application to model NEK.

Modified:
   Nalini Kumar  { UF CCMT }

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

#include "dstructs.h"
#include "params.h"
#include "time.h"


/* -------------------------- Vector Functions ----------------------------- */

/* Make a new 'vector' type and allocate memory for it. */
vector new_vector(int size)
{
  vector X = malloc(sizeof(vectortype));
  X->size = size;
  X->V = malloc(sizeof( dtype * ) * size);
  return X;
}

/* Free up the memory allocated for the vector X. */
void delete_vector(vector X)
{
  free(X->V);
  free(X);
}

/* Fill a vector with random numbers over [lower, upper) */
void random_fill_vector(vector X, dtype lower, dtype upper)
{
  int i;
  for (i = 0; i < (X->size); i++) {
    X->V[i] = ((dtype) rand() / (RAND_MAX)) * (upper - lower + 1) + lower;
  }
}

/* Return a newly-allocated random vector */
vector new_random_vector(int size, dtype lower, dtype upper)
{
  vector X = new_vector(size);
  random_fill_vector(X, lower, upper);
  return X;
}


/* -------------------------- Matrix Functions ----------------------------- */

/* Make a new 'matrix' type and allocate memory. Use: A->M[row][column]. */
matrix new_matrix(int rows, int cols)
{
  int i;
  matrix A = malloc(sizeof(matrixtype));
  A->rows = rows;
  A->cols = cols;
  A->M = malloc(sizeof( dtype * ) * rows);

  for (i = 0; i < rows; i++) {
    A->M[i] = malloc(sizeof( dtype * ) * cols);
  }

  return A;
}

/* Free up the memory allocated for the matrix A. */
void delete_matrix(matrix A)
{
  int row;
  for (row = 0; row<(A->rows); row++) { free(A->M[row]); }
  free(A->M);
  free(A);
}

/* Zero out the matrix A. */
void zero_matrix(matrix A)
{
  int row, col;
  for(row = 0; row<(A->rows); row++) {
    for(col = 0; col<(A->cols); col++) {
      A->M[row][col] = (dtype) 0;
    }
  }
}

/* Fill a matrix with random numbers over [lower, upper). */
void random_fill_matrix(matrix A, dtype lower, dtype upper)
{
  int row, col;
  for (row = 0; row < (A->rows); row++) {
    for (col = 0; col < (A->cols); col++) {
      A->M[row][col] = (dtype) rand() / RAND_MAX * (upper - lower + 1) + lower;
    }
  }
}

/* Return a newly-allocated random matrix. */
matrix new_random_matrix(int rows, int cols, dtype lower, dtype upper)
{
  matrix A = new_matrix(rows, cols);
  random_fill_matrix(A, lower, upper);
  return A;
}


/* -------------------------- Ternix Functions ----------------------------- */

/* Make a new 'ternix' type and allocate memory for it.
  Access is done by: A->T[row][column][layer]. */
ternix new_ternix(int rows, int cols, int layers)
{
  int i, j;
  ternix A = malloc(sizeof(ternixtype));
  A->rows = rows;
  A->cols = cols;
  A->layers = layers;
  A->T = malloc( sizeof( dtype * ) * rows );

  for (i = 0; i<rows; i++) {
    A->T[i] = malloc( sizeof( dtype * ) * cols);
    for (j = 0; j<cols; j++) {
      A->T[i][j] = malloc( sizeof( dtype ) * layers);
    }
  }

  return A;
}


/* Free up the memory allocated for the ternix A. */
void delete_ternix(ternix A)
{
  int row, col;
  for (row = 0; row<(A->rows); row++) {
    for (col = 0; col<(A->cols); col++) {
      free(A->T[row][col]);
    }
    free(A->T[row]);
  }
  free(A->T);
  free(A);
}


/* Zero out the ternix A. */
void zero_ternix(ternix A)
{
  int row, col, layer;
  for(row = 0; row<(A->rows); row++) {
    for(col = 0; col<(A->cols); col++) {
      for(layer = 0; layer<(A->layers); layer++) {
        A->T[row][col][layer] = (dtype) 0;
      }
    }
  }
}


/* Return a newly-allocated zeroed ternix. */
ternix new_zero_ternix(int rows, int cols, int layers)
{
  ternix A = new_ternix(rows, cols, layers);
  zero_ternix(A);
  return A;
}


/* Fill a ternix with random numbers over [lower, upper). */
void random_fill_ternix(ternix A, dtype lower, dtype upper)
{
  int row, col, layer;
  for (row = 0; row<(A->rows); row++) {
    for (col = 0; col<(A->cols); col++) {
      for(layer = 0; layer<(A->layers); layer++) {
        A->T[row][col][layer] = ((dtype) rand() / (RAND_MAX)) *
                                 (upper - lower + 1) + lower;
      }
    }
  }
}


/* Return a random newly-allocated ternix. */
ternix new_random_ternix(int rows, int cols, int layers,
                         dtype lower, dtype upper)
{
  ternix A = new_ternix(rows, cols, layers);
  random_fill_ternix(A, lower, upper);
  return A;
}


/* -------------------------- Element Functions ---------------------------- */

/* Return an element with PHYSICAL_PARAMTERS blocks of ELEMENT_SIZE,
   randomly filled with ternices over [lower, upper). */
element new_random_element(dtype lower, dtype upper, struct paramstype *params)
{
  int i;
  element A = malloc(sizeof(elementtype));

  A->B = malloc(sizeof( ternix * ) * params->PHYSICAL_PARAMS);

  for (i = 0; i < params->PHYSICAL_PARAMS; i++) {
    A->B[i] = new_random_ternix( params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE,
                                 lower, upper );
  }

  return A;
}


/* Return an element with PHYSICAL_PARAMTERS blocks of ELEMENT_SIZE. */
element new_zero_element(struct paramstype *params)
{
  int i;
  element A = malloc( sizeof(elementtype) );

  A->B = malloc( sizeof( ternix * ) * params->PHYSICAL_PARAMS );

  for (i = 0; i < params->PHYSICAL_PARAMS; i++) {
    A->B[i] = new_zero_ternix(params->ELEMENT_SIZE, params->ELEMENT_SIZE, params->ELEMENT_SIZE);
  }

  return A;
}


/* Frees up the memory allocated for the element A. */
void delete_element(element A, struct paramstype *params)
{
  int i;

  for (i = 0; i < params->PHYSICAL_PARAMS; i++) { delete_ternix( A->B[i] ); }

  free(A->B);
  free(A);
}



