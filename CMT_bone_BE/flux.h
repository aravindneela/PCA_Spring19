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


#ifndef FLUX_H_
#define FLUX_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "dstructs.h"


/* ---------------------------- Face Functions ----------------------------- */

/* Return a collection of faces from a set of elements, where the faces
   for each physical parameter have been clumped together in anticipation
   of a transfer operation. Possible values for arguments:

     axis:  {0, 1, 2}  |  X, Y, or Z
     sign:  {-1, 1}    |  Minus or Plus Face

   The resulting vector output is of size: PHYSICAL_PARAMS * FACE_SIZE
   multiplied by the number of elements on the face of interest. */
vector new_extracted_faces(element *elements, int axis, int sign, struct paramstype *params);

/* Same as above, but intended for the recv side, so not initialized. */
vector new_empty_faces(int axis, struct paramstype *params);


/* ------------------------ Faked CMT-Nek Operations ----------------------- */

/* Perform the R axis derivative operation, with kernel A and result C. */
void operation_dr(matrix A, ternix B, ternix C, struct paramstype *params);

/* Perform the S axis derivative operation, with kernel A and result C. */
void operation_ds(matrix A, ternix B, ternix C, struct paramstype *params);

/* Perform the T axis derivative operation, with kernel A and result C. */
void operation_dt(matrix A, ternix B, ternix C, struct paramstype *params);

/* Given Q, produce UR, US, and UT by faked transformation. HX, HY, and HZ
   are temporary space. RX is the list of transformation ternices. */
void operation_conv(ternix Q, ternix *RX, ternix Hx, ternix Hy, ternix Hz,
                    ternix Ur, ternix Us, ternix Ut, struct paramstype *params);

/* Add three ternices together and put the result in R. */
void operation_sum(ternix X, ternix Y, ternix Z, ternix R, struct paramstype *params);

/* Perform a faked Runge Kutta stage (no previous stage information used). */
void operation_rk(ternix Q, ternix R, struct paramstype *params);

#endif
