#include "sor.h"
#include <math.h>
#include <stdio.h>
#include "helper.h"

void sor(double omg, double dx, double dy, int imax, int jmax, double **P,
         double **RS, double *res, int **Flag, int no_fluid_cells) {
  int i, j;
  double rloc;
  double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

  /* SOR iteration */
  for (i = 1; i <= imax; i++) {
    for (j = 1; j <= jmax; j++) {
      if (Flag[i][j] & 1) {  // only for fluid cells
        P[i][j] = (1.0 - omg) * P[i][j] +
                  coeff * ((P[i + 1][j] + P[i - 1][j]) / (dx * dx) +
                           (P[i][j + 1] + P[i][j - 1]) / (dy * dy) - RS[i][j]);
      }
    }
  }

  /* compute the residual */
  rloc = 0;
  for (i = 1; i <= imax; i++) {
    for (j = 1; j <= jmax; j++) {
      if (Flag[i][j] & 1) {  // only for fluid cells
        rloc += ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx) +
                 (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy) -
                 RS[i][j]) *
                ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx) +
                 (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy) -
                 RS[i][j]);
      }
    }
  }
  rloc = rloc / no_fluid_cells;
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;

  /* set boundary values at boundary stripe */
  for (i = 1; i <= imax; i++) {
    P[i][0] = P[i][1];
    P[i][jmax + 1] = P[i][jmax];
  }
  for (j = 1; j <= jmax; j++) {
    P[0][j] = P[1][j];
    P[imax + 1][j] = P[imax][j];
  }

  // set boundary values at interior cells
  for (i = 1; i <= imax; i++) {
    for (j = 1; j <= jmax; j++) {
      if (!(Flag[i][j] & 1))
        switch ((Flag[i][j] >> SHIFT_EWSN) &
                EWSN_MASK) {  // extract if the our cell is fluid / bound and,
                              // if bound, where is fluid around
            // edges
          case N:
            P[i][j] = P[i][j + 1];
            break;
          case E:
            P[i][j] = P[i + 1][j];
            break;
          case W:
            P[i][j] = P[i - 1][j];
            break;
          case S:
            P[i][j] = P[i][j - 1];
            break;
            // corners
          case SE:
            P[i][j] = (P[i + 1][j] + P[i][j - 1]) / 2;
            break;
          case NW:
            P[i][j] = (P[i - 1][j] + P[i][j + 1]) / 2;
            break;
          case NE:
            P[i][j] = (P[i + 1][j] + P[i][j + 1]) / 2;
            break;
          case SW:
            P[i][j] = (P[i - 1][j] + P[i][j - 1]) / 2;
            break;
        }
    }
  }
}
