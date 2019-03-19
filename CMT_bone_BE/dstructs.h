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

#ifndef DSTRUCTS_H_
#define DSTRUCTS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "params.h"

/* -------------------------- Type Definitions ---------------------------- */

/* NOTE: the above definition does not know the number of physical parameters,
   so it is dependent on the macro PHYSICAL_PARAMS for the size of B. */


typedef double dtype; // dtype: internal data storage type for calculations

typedef struct {
  int size;
  dtype * V;
} vectortype, *vector;

typedef struct {
  int rows;
  int cols;
  dtype ** M;
} matrixtype, *matrix;

typedef struct {
  int rows;
  int cols;
  int layers;
  dtype *** T;
} ternixtype, *ternix;

typedef struct {
  ternix *B;
} elementtype, *element;

/* -------------------------- Vector Functions ----------------------------- */
	vector new_vector(int size);
	void delete_vector(vector X);
	void random_fill_vector(vector X, dtype lower, dtype upper);
	vector new_random_vector(int size, dtype lower, dtype upper);

/* -------------------------- Matrix Functions ----------------------------- */
	matrix new_matrix(int rows, int cols);
	void delete_matrix(matrix A);
	void zero_matrix(matrix A);
	void random_fill_matrix(matrix A, dtype lower, dtype upper);
	matrix new_random_matrix(int rows, int cols, dtype lower, dtype upper);

/* -------------------------- Ternix Functions ----------------------------- */
	ternix new_ternix(int rows, int cols, int layers);
	void delete_ternix(ternix A);
	void zero_ternix(ternix A);
	ternix new_zero_ternix(int rows, int cols, int layers);
	void random_fill_ternix(ternix A, dtype lower, dtype upper);
	ternix new_random_ternix(int rows, int cols, int layers,
                         dtype lower, dtype upper);

/* -------------------------- Element Functions ---------------------------- */
	element new_random_element(dtype lower, dtype upper, struct paramstype *params);
	element new_zero_element(struct paramstype *params);
	void delete_element(element A, struct paramstype *params);

#endif

