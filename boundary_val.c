#include "boundary_val.h"
#include <string.h>
#include "helper.h"

#define PI 3.14159265358979323846

void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flag,
                    char *problem, double **T) {
  // General problem: boundary values on the boundary stripe
  for (int j = 1; j <= jmax; j++) {
    // W bound & fluid to E (i = 0)
    if (Flag[1][j] & 1) {  // put values only if we have a fluid to E
      switch (Flag[0][j] & FLAGS_BC) {
        case NO_SLIP:
          U[0][j] = 0;         // no-slip (15)
          V[0][j] = -V[1][j];  // no-slip (16)
          break;
        case FREE_SLIP:
          U[0][j] = 0;        // free-slip (1.1)
          V[0][j] = V[1][j];  // free-slip (1.2)
          break;
        case OUTFLOW:
          U[0][j] = U[1][j];  // outflow (1.3)
          V[0][j] = V[1][j];  // outflow (1.3)
          break;
      }
      T[0][j] = T[1][j];  // adiabatic (2.10)
    }

    // E bound & fluid to W (i = imax)
    if (Flag[imax][j] & 1) {  // put values only if we have a fluid to W
      switch (Flag[imax + 1][j] & FLAGS_BC) {
        case NO_SLIP:
          U[imax][j] = 0;                // no-slip (15)
          V[imax + 1][j] = -V[imax][j];  // no-slip  (16)
          break;
        case FREE_SLIP:
          U[imax][j] = 0;               // free-slip (1.1)
          V[imax + 1][j] = V[imax][j];  // free-slip (1.2)
          break;
        case OUTFLOW:
          U[imax][j] = U[imax - 1][j];  // outflow (1.3)
          V[imax + 1][j] = V[imax][j];  // outflow (1.3)
          break;
      }
      T[imax + 1][j] = T[imax][j];  // adiabatic (2.11)
    }
  }

  for (int i = 1; i <= imax; i++) {
    // S bound & fluid to N (j = 0)
    if (Flag[i][1] & 1) {  // put values only if we have a fluid to N
      switch (Flag[i][0] & FLAGS_BC) {
        case NO_SLIP:
          V[i][0] = 0;         // no-slip (15)
          U[i][0] = -U[i][1];  // no-slip (16)
          break;
        case FREE_SLIP:
          V[i][0] = 0;        // free-slip (1.1)
          U[i][0] = U[i][1];  // free-slip (1.2)
          break;
        case OUTFLOW:
          V[i][0] = V[i][1];  // outflow (1.3)
          U[i][0] = U[i][1];  // outflow (1.3)
          break;
      }
      T[i][0] = T[i][1];  // adiabatic (2.12)
    }

    // N bound & fluid to S (j = jmax)
    if (Flag[i][jmax] & 1) {  // put values only if we have a fluid to S
      switch (Flag[i][jmax + 1] & FLAGS_BC) {
        case NO_SLIP:
          V[i][jmax] = 0;                // no-slip (15)
          U[i][jmax + 1] = -U[i][jmax];  // no-slip (16)
          break;
        case FREE_SLIP:
          V[i][jmax] = 0;               // free-slip (1.1)
          U[i][jmax + 1] = U[i][jmax];  // free-slip (1.2)
          break;
        case OUTFLOW:
          V[i][jmax] = V[i][jmax - 1];  // outflow (1.3)
          U[i][jmax + 1] = U[i][jmax];  // outflow (1.3)
          break;
      }
      T[i][jmax + 1] = T[i][jmax];  // adiabatic (2.13)
    }
  }
}

void specialboundaryvalues(int imax, int jmax, double **U, double **V,
                           int **Flag, char *problem, double **T, double T_h,
                           double T_c, double t, int special_input) {
  // Specify the inflow conditions at West edge - only for the problems
  // "channel-bfs" and "karman_vortex"
  if (!strcmp(problem, "channel-bfs")) {
    // West edge has last 1/2 no-slip BC and first 1/2 Inflow BC; we set only
    // the input conditions
    for (int j = jmax; j >= 1 && (Flag[1][j] & 1); --j) {
      U[0][j] = 2 - U[1][j];
      V[0][j] = 0;
    }
  } else if (!strcmp(problem, "karman_vortex")) {
    // West edge has Inflow BC
    for (int j = 1; j <= jmax; ++j) {
      U[0][j] = 2 - U[1][j];
      V[0][j] = 0;
    }
  } else if (!strcmp(problem, "cavity100")) {
    double sin_arg;
    if (special_input)  // sine at top edge
      sin_arg = sin(2 * PI * t / 10);
    else
      sin_arg = 1;
    for (int i = 1; i <= imax; i++) {
      // formula (15)
      V[i][jmax] = 0;
      // U[i][jmax + 1] = 2 - U[i][jmax];  // moving wall
      U[i][jmax + 1] = sin_arg * 2 - U[i][jmax];  // time-dependent moving wall
    }
  } else if (!strcmp(problem, "fluid_trap") ||  // fluid trap
             !strcmp(problem,
                     "natural_conv")) {  // Natural convection problems
    // only no-slip boundary conditions on all edges
    // temperature Dirichlet for left & right
    for (int j = 1; j <= jmax; j++) {
      // Dirichlet for temperature on left boundary (2.6)
      T[0][j] = 2 * T_h - T[1][j];
      // Dirichlet for temperature on right boundary (2.7)
      T[imax + 1][j] = 2 * T_c - T[imax][j];
    }
  } else if (!strcmp(problem, "Rayleigh_Benard")) {
    for (int i = 1; i <= imax; i++) {
      // Dirichlet for hot temperature on bottom boundary (2.8)
      T[i][0] = 2 * T_h - T[i][1];
      // Dirichlet for cold temperature on top boundary (2.9)
      T[i][jmax + 1] = 2 * T_c - T[i][jmax];
    }
  }

  // Obstacles / Interior cells
  // edges or corners; (+) check if our cell is fluid or boundary
  if (!strcmp(problem, "karman_vortex") || !strcmp(problem, "channel-bfs") ||
      !strcmp(problem, "fluid_trap")) {
    for (int i = 1; i <= imax; i++) {
      for (int j = 1; j <= jmax; j++) {
        if (!(Flag[i][j] & 1))
          switch ((Flag[i][j] >> SHIFT_EWSN) &
                  EWSN_MASK) {  // extract if the our cell is fluid / bound
                                // and, if bound, where is fluid around
              // edges
              // N and E first - for "Flow over a Step" we don't have W and S
            case N:
              V[i][j] = 0;
              U[i - 1][j] = -U[i - 1][j + 1];
              U[i][j] = -U[i][j + 1];
              T[i][j] = T[i][j + 1];  // adiabatic condition
              break;
            case E:
              U[i][j] = 0;
              V[i][j] = -V[i + 1][j];
              V[i][j - 1] = -V[i + 1][j - 1];
              T[i][j] = T[i + 1][j];  // adiabatic condition
              break;
            case W:
              U[i - 1][j] = 0;
              V[i][j - 1] = -V[i - 1][j - 1];
              V[i][j] = -V[i - 1][j];
              T[i][j] = T[i - 1][j];  // adiabatic condition
              break;
            case S:
              V[i][j - 1] = 0;
              U[i][j] = -U[i][j - 1];
              U[i - 1][j] = -U[i - 1][j - 1];
              T[i][j] = T[i][j + 1];  // adiabatic condition
              break;
              // corners
            case SE:
              U[i][j] = 0;
              V[i][j - 1] = 0;
              V[i][j] = -V[i + 1][j];
              U[i - 1][j] = -U[i - 1][j - 1];
              T[i][j] = (T[i + 1][j] + T[i][j - 1]) / 2;  // adiabatic condition
              break;
            case NW:
              U[i - 1][j] = 0;
              V[i][j] = 0;
              U[i][j] = -U[i][j + 1];
              V[i][j - 1] = -V[i - 1][j - 1];
              T[i][j] = (T[i - 1][j] + T[i][j + 1]) / 2;  // adiabatic condition
              break;
            case NE:
              U[i][j] = 0;
              V[i][j] = 0;
              U[i - 1][j] = -U[i - 1][j + 1];
              V[i][j - 1] = -V[i + 1][j - 1];
              T[i][j] = (T[i + 1][j] + T[i][j + 1]) / 2;  // adiabatic condition
              break;
            case SW:
              U[i - 1][j] = 0;
              V[i][j - 1] = 0;
              U[i][j] = -U[i][j - 1];
              V[i][j] = -V[i - 1][j];
              T[i][j] = (T[i - 1][j] + T[i][j - 1]) / 2;  // adiabatic condition
              break;
          }
      }
    }
  }
}
