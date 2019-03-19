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

#include "params.h"
#include "dstructs.h"


/* ------------------------------------------------------------------------- */
/* ---------------------------- Face Functions ----------------------------- */
/* ------------------------------------------------------------------------- */

vector new_extracted_faces(element *elements, int axis, int sign, struct paramstype *params)
/* Return a collection of faces from a set of elements, where the faces
   for each physical parameter have been clumped together in anticipation
   of a transfer operation. Possible values for arguments:

     axis:  {0, 1, 2}  |  X, Y, or Z
     sign:  {-1, 1}    |  Minus or Plus Face

   The resulting vector output is of size: PHYSICAL_PARAMS * FACE_SIZE
   multiplied by the number of elements on the face of interest. */
{
  int i, b, e, row, col, layer, plane, EoF;

  EoF=0;
  switch (axis) { /* EoF: elements on face */
  case 0: EoF = params->ELEMENTS_ON_X_FACE; break;
  case 1: EoF = params->ELEMENTS_ON_Y_FACE; break;
  case 2: EoF = params->ELEMENTS_ON_Z_FACE; break; }

  vector faces = new_vector(EoF * params->PHYSICAL_PARAMS * params->FACE_SIZE);

  /* ---------------------------- Extraction ------------------------------- */

  /* The plane is the index of the lower or upper face. */
  plane = (sign > 0) ? params->ELEMENT_SIZE - 1 : 0;

  /* Index in the output vector */
  i = 0;

  /* For each element owned by this rank: */
  for (e = 0; e < EoF; e++) {

    /* For each block in the element: */
    for (b = 0; b < params->PHYSICAL_PARAMS; b++) {

      if ( axis == 0 ) {
        /* If this is the X axis, the plane is on the row dimension. */

        for (col = 0; col < params->ELEMENT_SIZE; col++) {
          for (layer = 0; layer < params->ELEMENT_SIZE; layer++) {
            faces->V[i] = elements[e]->B[b]->T[plane][col][layer]; i++; } }

      } else if ( axis == 1 ) {
        /* If this is the Y axis, the plane is on the column dimension. */

        for (row = 0; row < params->ELEMENT_SIZE; row++) {
          for (layer = 0; layer < params->ELEMENT_SIZE; layer++) {
            faces->V[i] = elements[e]->B[b]->T[row][plane][layer]; i++; } }

      } else if ( axis == 2 ) {
        /* If this is the Z axis, the plane is on the layer dimension. */

        for (row = 0; row < params->ELEMENT_SIZE; row++) {
          for (col = 0; col < params->ELEMENT_SIZE; col++) {
            faces->V[i] = elements[e]->B[b]->T[row][col][plane]; i++; } }
      }
    }
  }
  return faces;
}



vector new_empty_faces(int axis, struct paramstype *params)
/* Same as above, but intended for the recv side, so not initialized. */
{
  int EoF=0;

  switch (axis) { /* EoF: elements on face */
  case 0: EoF = params->ELEMENTS_ON_X_FACE; break;
  case 1: EoF = params->ELEMENTS_ON_Y_FACE; break;
  case 2: EoF = params->ELEMENTS_ON_Z_FACE; break; }

  return new_vector(EoF * params->PHYSICAL_PARAMS * params->FACE_SIZE);
}

/* ------------------------------------------------------------------------- */
/* ------------------------ Faked CMT-Nek Operations ----------------------- */
/* ------------------------------------------------------------------------- */

void operation_dr(matrix A, ternix B, ternix C, struct paramstype *params)
/* Perform the R axis derivative operation, with kernel A and result C. */
{
  zero_ternix(C);

  int k, j, i, g;

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {
        for (g = 0; g < params->ELEMENT_SIZE; g++) {
          C->T[i][j][k] += A->M[i][g] * B->T[g][j][k]; } } } }
}

void operation_ds(matrix A, ternix B, ternix C, struct paramstype *params)
/* Perform the S axis derivative operation, with kernel A and result C. */
{
  zero_ternix(C);

  int k, j, i, g;

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {
        for (g = 0; g < params->ELEMENT_SIZE; g++) {
          C->T[i][j][k] += A->M[j][g] * B->T[i][g][k]; } } } }
}

void operation_dt(matrix A, ternix B, ternix C, struct paramstype *params)
/* Perform the T axis derivative operation, with kernel A and result C. */
{
  zero_ternix(C);

  int k, j, i, g;

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {
        for (g = 0; g < params->ELEMENT_SIZE; g++) {
          C->T[i][j][k] += A->M[k][g] * B->T[i][j][g]; } } } }
}

void operation_conv(ternix Q, ternix *RX, ternix Hx, ternix Hy, ternix Hz,
                    ternix Ur, ternix Us, ternix Ut, struct paramstype *params)
/* Given Q, produce UR, US, and UT by faked transformation. HX, HY, and HZ
   are temporary space. RX is the list of transformation ternices. */
{

  /* Generate three random constants. */

  dtype a = ((dtype) rand() / (RAND_MAX));
  dtype b = ((dtype) rand() / (RAND_MAX));
  dtype c = ((dtype) rand() / (RAND_MAX));

  int k, j, i;

  /* First, make HX, HY, and HZ, which are used in the next step. */

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {
        Hx->T[i][j][k] = a * Q->T[i][j][k];
        Hy->T[i][j][k] = b * Q->T[i][j][k];
        Hz->T[i][j][k] = c * Q->T[i][j][k];
      }
    }
  }

  /* Then, produce our outputs using HX, HY, HZ, and RX. */

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {

        /* Generate UR. */
        Ur->T[i][j][k] = ( RX[0]->T[i][j][k] * Hx->T[i][j][k] +
                           RX[1]->T[i][j][k] * Hy->T[i][j][k] +
                           RX[2]->T[i][j][k] * Hz->T[i][j][k] );

        /* Generate UR. */
        Us->T[i][j][k] = ( RX[3]->T[i][j][k] * Hx->T[i][j][k] +
                           RX[4]->T[i][j][k] * Hy->T[i][j][k] +
                           RX[5]->T[i][j][k] * Hz->T[i][j][k] );

        /* Generate UT. */
        Ut->T[i][j][k] = ( RX[6]->T[i][j][k] * Hx->T[i][j][k] +
                           RX[7]->T[i][j][k] * Hy->T[i][j][k] +
                           RX[8]->T[i][j][k] * Hz->T[i][j][k] );

      }
    }
  }
}

void operation_sum(ternix X, ternix Y, ternix Z, ternix R, struct paramstype *params)
/* Add three ternices together and put the result in R. */
{
  int k, j, i;

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {
        R->T[i][j][k] = X->T[i][j][k] + Y->T[i][j][k] + Z->T[i][j][k];
      }
    }
  }
}

void operation_rk(ternix Q, ternix R, struct paramstype *params)
/* Perform a faked Runge Kutta stage (no previous stage information used). */
{
  int k, j, i;

  for (k = 0; k < params->ELEMENT_SIZE; k++) {
    for (j = 0; j < params->ELEMENT_SIZE; j++) {
      for (i = 0; i < params->ELEMENT_SIZE; i++) {
        Q->T[i][j][k] = ( R->T[i][j][k] * 0.5 +
                          R->T[i][j][k] * 0.25 +
                          Q->T[i][j][k] * 0.5 );
      }
    }
  }
}


