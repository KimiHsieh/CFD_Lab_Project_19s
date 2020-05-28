#include "uvp.h"
#include <math.h>
#include "helper.h"

double pow2(double a) { return a * a; }

double max_abs_matrix(double **M, int nrl, int nrh, int ncl, int nch) {
  double max_abs = 0.;
  for (int i = nrl; i <= nrh; ++i)
    for (int j = ncl; j <= nch; ++j)
      if (fabs(M[i][j]) > max_abs) max_abs = fabs(M[i][j]);
  return max_abs;
}

void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
                  int imax, int jmax, double **U, double **V, double Pr) {
  // formula (14)
  if (tau > 0) {
    double max_abs_U = max_abs_matrix(U, 0, imax, 0, jmax + 1);
    double max_abs_V = max_abs_matrix(V, 0, imax + 1, 0, jmax);
    double aux = Re / 2 / (1. / pow2(dx) + 1. / pow2(dy));
    *dt = tau * fmin(fmin(dx / max_abs_U, dy / max_abs_V), aux);
    if (Pr < 1.)  // condition if temperature problem, with Pr < 1.0
      *dt = fmin(*dt, Pr * aux);
  }
}

void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
                  double dx, double dy, int imax, int jmax, double **U,
                  double **V, double **F, double **G, int **Flag, double **T,
                  double beta) {
  // formulae (10),(11) for F, G and formula (18) at the boundary
  // alpha is gamma from the formulae
  // ===== compute F
  for (int i = 1; i <= imax - 1; ++i)
    for (int j = 1; j <= jmax; ++j) {
      if (Flag[i][j] & 1) {                    // I am fluid cell
        if ((Flag[i][j] >> SHIFT_EWSN) & E) {  // I have fluid at E => compute F

          // formula (5)
          double d2udx2 = (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / pow2(dx);
          double d2udy2 = (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / pow2(dy);
          double du2dx = 1. / dx *
                             (pow2((U[i][j] + U[i + 1][j]) / 2) -
                              pow2((U[i - 1][j] + U[i][j]) / 2)) +
                         alpha * (1. / dx) *
                             (fabs(U[i][j] + U[i + 1][j]) / 2 *
                                  (U[i][j] - U[i + 1][j]) / 2 -
                              fabs(U[i - 1][j] + U[i][j]) / 2 *
                                  (U[i - 1][j] - U[i][j]) / 2);
          double duvdy =
              1. / dy *
                  ((V[i][j] + V[i + 1][j]) / 2 * (U[i][j] + U[i][j + 1]) / 2 -
                   (V[i][j - 1] + V[i + 1][j - 1]) / 2 *
                       (U[i][j - 1] + U[i][j]) / 2) +
              alpha * (1. / dy) *
                  (fabs(V[i][j] + V[i + 1][j]) / 2 * (U[i][j] - U[i][j + 1]) /
                       2 -
                   fabs(V[i][j - 1] + V[i + 1][j - 1]) / 2 *
                       (U[i][j - 1] - U[i][j]) / 2);
          // formula (10) and (2.14)
          F[i][j] = (1. / Re * (d2udx2 + d2udy2) - du2dx - duvdy + GX) -
                    (beta / 2 * (T[i][j] + T[i + 1][j]) * GX);
        }
      }
    }
  // ===== end of compute F

  // ===== compute G
  for (int i = 1; i <= imax; ++i)
    for (int j = 1; j <= jmax - 1; ++j) {
      if (Flag[i][j] & 1) {                    // I am fluid cell
        if ((Flag[i][j] >> SHIFT_EWSN) & N) {  // I have fluid at N => compute G
                                               // formula (6)
          double d2vdx2 = (V[i + 1][j] - 2 * V[i][j] + V[i - 1][j]) / pow2(dx);
          double d2vdy2 = (V[i][j + 1] - 2 * V[i][j] + V[i][j - 1]) / pow2(dy);
          double dv2dy = 1. / dy *
                             (pow2((V[i][j] + V[i][j + 1]) / 2) -
                              pow2((V[i][j - 1] + V[i][j]) / 2)) +
                         alpha * (1. / dy) *
                             (fabs(V[i][j] + V[i][j + 1]) / 2 *
                                  (V[i][j] - V[i][j + 1]) / 2 -
                              fabs(V[i][j - 1] + V[i][j]) / 2 *
                                  (V[i][j - 1] - V[i][j]) / 2);
          double duvdx =
              1. / dx *
                  ((U[i][j] + U[i][j + 1]) / 2 * (V[i][j] + V[i + 1][j]) / 2 -
                   (U[i - 1][j] + U[i - 1][j + 1]) / 2 *
                       (V[i - 1][j] + V[i][j]) / 2) +
              alpha * (1. / dx) *
                  (fabs(U[i][j] + U[i][j + 1]) / 2 * (V[i][j] - V[i + 1][j]) /
                       2 -
                   fabs(U[i - 1][j] + U[i - 1][j + 1]) / 2 *
                       (V[i - 1][j] - V[i][j]) / 2);
          // formula (11) and (2.15)
          G[i][j] = (1. / Re * (d2vdx2 + d2vdy2) - duvdx - dv2dy + GY) -
                    (beta / 2 * (T[i][j] + T[i][j + 1]) * GY);
        }
      }
    }
  // ===== end of compute G
}

void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
                  double **F, double **G, double **RS, int **Flag, double **U,
                  double **V) {
  // RHS of formula (12)
  for (int i = 1; i <= imax; ++i)
    for (int j = 1; j <= jmax; ++j)
      if (Flag[i][j] & 1)  // only for fluid cells
        RS[i][j] =
            1 / dt *
                ((U[i][j] - U[i - 1][j]) / dx + (V[i][j] - V[i][j - 1]) / dy) +
            (F[i][j] - F[i - 1][j]) / dx + (G[i][j] - G[i][j - 1]) / dy;
}

void calculate_uv_stagek(double dx, double dy, int imax, int jmax,
                         double ***Udot, double ***Vdot, double **F, double **G,
                         double **P, int **Flag, int k) {
  // formula (8)
  for (int i = 1; i <= imax - 1; ++i)
    for (int j = 1; j <= jmax; ++j)
      if ((Flag[i][j] & 1) &&                // fluid cell
          ((Flag[i][j] >> SHIFT_EWSN) & E))  // only edge to East
        Udot[k][i][j] = F[i][j] - 1. / dx * (P[i + 1][j] - P[i][j]);
  // formula (9)
  for (int i = 1; i <= imax; ++i)
    for (int j = 1; j <= jmax - 1; ++j)
      if ((Flag[i][j] & 1) &&                // fluid cell
          ((Flag[i][j] >> SHIFT_EWSN) & N))  // only edge to North
        Vdot[k][i][j] = G[i][j] - 1. / dy * (P[i][j + 1] - P[i][j]);
}

void calculate_uv_tilde(double dt, int imax, int jmax, double **U, double **V,
                        double **U_tilde, double **V_tilde, double ***Udot,
                        double ***Vdot, int method, int k) {
  switch (method) {
    case 0:  // Explicit Euler (only 1 stage)
      switch (k) {
        case 0:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j];
            }
          break;
        default:
          printf("Error! Stage %d not defined in method %d!", k, method);
          break;
      }
      break;

    case 1:  // RK2 (Heun) - 2 stages
      switch (k) {
        case 0:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j];
            }
          break;
        case 1:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * Udot[0][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * Vdot[0][i][j];
            }
          break;
        default:
          printf("Error! Stage %d not defined in method %d!", k, method);
          break;
      }
      break;

    case 2:  // RK4 - 4 stages
      switch (k) {
        case 0:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j];
            }
          break;

        case 1:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 0.5 * Udot[0][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 0.5 * Vdot[0][i][j];
            }
          break;

        case 2:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 0.5 * Udot[1][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 0.5 * Vdot[1][i][j];
            }
          break;
        case 3:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * Udot[2][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * Vdot[2][i][j];
            }
          break;

        default:
          printf("Error! Stage %d not defined in method %d!", k, method);
          break;
      }
      break;

    case 3:  // ode23 - 3/4 stages
      switch (k) {
        case 0:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j];
            }
          break;

        case 1:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 0.5 * Udot[0][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 0.5 * Vdot[0][i][j];
            }
          break;

        case 2:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 0.75 * Udot[1][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 0.75 * Vdot[1][i][j];
            }
          break;

        case 3:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] =
                  U[i][j] + dt *
                                (2. * Udot[0][i][j] + 3. * Udot[1][i][j] +
                                 4. * Udot[2][i][j]) /
                                9.;
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] =
                  V[i][j] + dt *
                                (2. * Vdot[0][i][j] + 3. * Vdot[1][i][j] +
                                 4. * Vdot[2][i][j]) /
                                9.;
            }
          break;

        default:
          printf("Error! Stage %d not defined in method %d!", k, method);
          break;
      }
      break;

    case 4:  // ode45 - 4/6 stages
      switch (k) {
        case 0:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j];
            }
          break;

        case 1:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 0.25 * Udot[0][i][j];
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 0.25 * Vdot[0][i][j];
            }
          break;

        case 2:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 3. * Udot[0][i][j] / 32. +
                              dt * 9. * Udot[1][i][j] / 32.;
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 3. * Vdot[0][i][j] / 32. +
                              dt * 9. * Vdot[1][i][j] / 32.;
            }
          break;

        case 3:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] =
                  U[i][j] + dt *
                                (1932. * Udot[0][i][j] - 7200. * Udot[1][i][j] +
                                 7296. * Udot[2][i][j]) /
                                2197.;
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] =
                  V[i][j] + dt *
                                (1932. * Vdot[0][i][j] - 7200. * Vdot[1][i][j] +
                                 7296. * Vdot[2][i][j]) /
                                2197.;
            }
          break;

        case 4:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] + dt * 439. * Udot[0][i][j] / 216. -
                              dt * 8. * Udot[1][i][j] +
                              dt * 3680. * Udot[2][i][j] / 513. -
                              dt * 845. * Udot[3][i][j] / 4104.;
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] + dt * 439. * Vdot[0][i][j] / 216. -
                              dt * 8. * Vdot[1][i][j] +
                              dt * 3680. * Vdot[2][i][j] / 513. -
                              dt * 845. * Vdot[3][i][j] / 4104.;
            }
          break;

        case 5:
          for (int i = 0; i <= imax; ++i)
            for (int j = 0; j <= jmax + 1; ++j) {
              U_tilde[i][j] = U[i][j] - dt * 8. * Udot[0][i][j] / 27. +
                              dt * 2. * Udot[1][i][j] -
                              dt * 3544. * Udot[2][i][j] / 2565. +
                              dt * 1859. * Udot[3][i][j] / 4104. -
                              dt * 11. * Udot[4][i][j] / 40.;
            }

          for (int i = 0; i <= imax + 1; ++i)
            for (int j = 0; j <= jmax; ++j) {
              V_tilde[i][j] = V[i][j] - dt * 8. * Vdot[0][i][j] / 27. +
                              dt * 2. * Vdot[1][i][j] -
                              dt * 3544. * Vdot[2][i][j] / 2565. +
                              dt * 1859. * Vdot[3][i][j] / 4104. -
                              dt * 11. * Vdot[4][i][j] / 40.;
            }
          break;

        default:
          printf("Error! Stage %d not defined in method %d!", k, method);
          break;
      }
      break;

    default:
      printf("Error! Stage %d not defined in method %d!", k, method);
      break;
  }
}

void calculate_uv(double dt, int imax, int jmax, double **U, double **V,
                  double ***Udot, double ***Vdot, int **Flag, int method) {
  switch (method) {
    case 0:  // Explicit Euler (only 1 stage)

      // formula (8)
      for (int i = 1; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax; ++j)
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & E))  // only edge to East
            U[i][j] = U[i][j] + dt * Udot[0][i][j];
      // formula (9)
      for (int i = 1; i <= imax; ++i)
        for (int j = 1; j <= jmax - 1; ++j)
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & N))  // only edge to North
            V[i][j] = V[i][j] + dt * Vdot[0][i][j];
      break;

    case 1:  // RK2 (Heun) - 2 stages
      for (int i = 1; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax; ++j)
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & E))  // only edge to East
            U[i][j] = U[i][j] + dt * 0.5 * (Udot[0][i][j] + Udot[1][i][j]);
      // formula (9)
      for (int i = 1; i <= imax; ++i)
        for (int j = 1; j <= jmax - 1; ++j)
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & N))  // only edge to North
            V[i][j] = V[i][j] + dt * 0.5 * (Vdot[0][i][j] + Vdot[1][i][j]);

      break;

    case 2:  // RK4 - 4 stages
      for (int i = 1; i <= imax - 1; ++i)
        for (int j = 1; j <= jmax; ++j) {
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & E))  // only edge to East
            U[i][j] = U[i][j] + dt *
                                    (Udot[0][i][j] + 2. * Udot[1][i][j] +
                                     2. * Udot[2][i][j] + Udot[3][i][j]) /
                                    6.;
        }

      for (int i = 1; i <= imax; ++i)
        for (int j = 1; j <= jmax - 1; ++j) {
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & N))  // only edge to North
            V[i][j] = V[i][j] + dt *
                                    (Vdot[0][i][j] + 2. * Vdot[1][i][j] +
                                     2. * Vdot[2][i][j] + Vdot[3][i][j]) /
                                    6.;
        }
      break;

    case 3:  // ode 23
      for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax; ++j) {
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & E))  // only edge to East
            U[i][j] = U[i][j] + dt *
                                    (2. * Udot[0][i][j] + 3. * Udot[1][i][j] +
                                     4. * Udot[2][i][j]) /
                                    9.;
        }
      }

      for (int i = 1; i <= imax; ++i)
        for (int j = 1; j <= jmax - 1; ++j) {
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & N))  // only edge to North
            V[i][j] = V[i][j] + dt *
                                    (2. * Vdot[0][i][j] + 3. * Vdot[1][i][j] +
                                     4. * Vdot[2][i][j]) /
                                    9.;
        }
      break;

    case 4:  // ode 45
      for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax; ++j) {
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & E))  // only edge to East
            U[i][j] = U[i][j] + dt * 25. * Udot[0][i][j] / 216. +
                      dt * 1408. * Udot[2][i][j] / 2565. +
                      dt * 2197. * Udot[3][i][j] / 4104. -
                      dt * Udot[4][i][j] / 5.;
        }
      }

      for (int i = 1; i <= imax; ++i)
        for (int j = 1; j <= jmax - 1; ++j) {
          if ((Flag[i][j] & 1) &&                // fluid cell
              ((Flag[i][j] >> SHIFT_EWSN) & N))  // only edge to North
            V[i][j] = V[i][j] + dt * 25. * Vdot[0][i][j] / 216. +
                      dt * 1408. * Vdot[2][i][j] / 2565. +
                      dt * 2197. * Vdot[3][i][j] / 4104. -
                      dt * Vdot[4][i][j] / 5.;
        }
      break;

    default:
      printf("Error! Method %d not defined!", method);
      break;
  }
}

void calculate_adaptive_dt(double *dt, double dt_ref, double ***Udot,
                           double ***Vdot, int imax, int jmax, double toler,
                           int *indicator, double **error_U, double **error_V,
                           int method) {
  int p;               // order of the method
  double alpha = 0.9;  // safety factor
  switch (method) {
    case 3:   // ode 23
      p = 3;  // method of third order
      for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax; ++j) {
          error_U[i][j] = (-5. * Udot[0][i][j] + 6. * Udot[1][i][j] +
                           8. * Udot[2][i][j] - 9. * Udot[3][i][j]) /
                          72.;
        }
      }
      for (int i = 1; i <= imax; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
          error_V[i][j] = (-5. * Vdot[0][i][j] + 6. * Vdot[1][i][j] +
                           8. * Vdot[2][i][j] - 9. * Vdot[3][i][j]) /
                          72.;
        }
      }
      break;

    case 4:   // ode 45
      p = 5;  // method of fifth order
      for (int i = 1; i <= imax - 1; ++i) {
        for (int j = 1; j <= jmax; ++j) {
          error_U[i][j] = Udot[0][i][j] / 360. - 128. * Udot[2][i][j] / 4275. -
                          2197. * Udot[3][i][j] / 75240. + Udot[4][i][j] / 50. +
                          2. * Udot[5][i][j] / 55.;
        }
      }
      for (int i = 1; i <= imax; ++i) {
        for (int j = 1; j <= jmax - 1; ++j) {
          error_V[i][j] = Vdot[0][i][j] / 360. - 128. * Vdot[2][i][j] / 4275. -
                          2197. * Vdot[3][i][j] / 75240. + Vdot[4][i][j] / 50. +
                          2. * Vdot[5][i][j] / 55.;
        }
      }
      break;

    default:
      printf(
          "Error! Method %d doesn't have adaptive time-step version "
          "implemented!\n",
          method);
      break;
  }
  double max_abs_error_U = max_abs_matrix(error_U, 1, imax - 1, 1, jmax);
  double max_abs_error_V = max_abs_matrix(error_V, 1, imax, 1, jmax - 1);
  double max_abs_error = fmax(max_abs_error_U, max_abs_error_V);
  *indicator = 1;  // default - no need to adapt time
  double dt_new = (*dt) * alpha * pow((toler / max_abs_error), 1. / (p + 1));
  *dt = fmin(dt_new, dt_ref);  // using dt_ref to satisfy the CFL condition
  if (max_abs_error > toler) {
    printf(
        "Rejected time-step: old dt = %.7f; new dt = %.7f; error = %.5f; toler "
        "= "
        "%.5f; dt_ref = "
        "%.5f\n",
        *dt, dt_new, max_abs_error, toler, dt_ref);
    *indicator = 0;  // repeat the steps
  } else {
    printf(
        "Accepted time-step: old dt = %.7f; new dt = %.7f; error = %.5f; toler "
        "= "
        "%.5f; dt_ref = "
        "%.5f\n",
        *dt, dt_new, max_abs_error, toler, dt_ref);
  }
}

void calculate_temp(double Re, double Pr, double alpha, double dt, double dx,
                    double dy, int imax, int jmax, double **U, double **V,
                    double **T, double **T_temp, int **Flag) {
  // formula (2.4)
  for (int i = 1; i <= imax; ++i) {
    for (int j = 1; j <= jmax; ++j) {
      if (Flag[i][j] & 1) {  // only for fluid cells
        double duTdx = 1. / dx *
                           (U[i][j] * (T[i][j] + T[i + 1][j]) / 2 -
                            (U[i - 1][j] * (T[i - 1][j] + T[i][j]) / 2)) +
                       alpha * (1. / dx) *
                           (fabs(U[i][j]) * (T[i][j] - T[i + 1][j]) / 2 -
                            fabs(U[i - 1][j]) * (T[i - 1][j] - T[i][j]) / 2);

        double dvTdy = 1. / dy *
                           (V[i][j] * (T[i][j] + T[i][j + 1]) / 2 -
                            (V[i][j - 1] * (T[i][j - 1] + T[i][j]) / 2)) +
                       alpha * (1. / dy) *
                           (fabs(V[i][j]) * (T[i][j] - T[i][j + 1]) / 2 -
                            fabs(V[i][j - 1]) * (T[i][j - 1] - T[i][j]) / 2);

        double d2Tdx2 = (T[i + 1][j] - 2 * T[i][j] + T[i - 1][j]) / pow2(dx);
        double d2Tdy2 = (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) / pow2(dy);

        T_temp[i][j] =
            T[i][j] + dt * (1. / (Re * Pr) * (d2Tdx2 + d2Tdy2) - duTdx - dvTdy);
      }
    }
  }
  for (int i = 1; i <= imax; ++i) {
    for (int j = 1; j <= jmax; ++j) {
      if (Flag[i][j] & 1) {  // only for fluid cells
        T[i][j] = T_temp[i][j];
      }
    }
  }
}
