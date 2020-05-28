#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "boundary_val.h"
#include "helper.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "visual.h"

void copy_matrices_dot(double ***Udot, double ***Vdot, int dest, int source,
                       int imax, int jmax) {
  for (int i = 1; i <= imax - 1; ++i)
    for (int j = 1; j <= jmax; ++j) Udot[dest][i][j] = Udot[source][i][j];

  for (int i = 1; i <= imax; ++i)
    for (int j = 1; j <= jmax - 1; ++j) Vdot[dest][i][j] = Vdot[source][i][j];
}

int main(int argn, char **args) {
  // time measurement structure
  clock_t clocks;
  clocks = clock();

  int imax, jmax, itermax, it, n, index_dt_value, num_fluid_cells, method,
      special_input;
  int **Flag, **Geom;
  double Re, UI, VI, PI, TI, GX, GY, t_end, xlength, ylength, dt, dt_ref, dx,
      dy, alpha, omg, tau, eps, dt_value, T_c, T_h, Pr, beta, toler;
  double res, t;
  double **U, **V, **P, **F, **G, **RS, **T;
  const char *problemName;
  char *szProblem, *szFileName, *dirName;
  char *problem, *geometry;
  // name of input file
  if (argn == 2)
    problemName = args[1];
  else {
    printf(
        "\nError! Please use the format:"
        "\t ./sim [problem_name], \n "
        "where [problem_name] is one of our defined problems:\n"
        "\t- cavity100 \t\t\t stationary problem: Driven cavity (analyzed in "
        "WS1)\n"
        "\t- karman_vortex \t\t non-stationary problem: The Karman Vortex "
        "Street (analyzed in WS2)\n"
        "\t- channel-bfs \t\t\t stationary problem: Flow over a Step (analyzed "
        "in WS2)\n"
        "or your own configuration file put in the 'configs' folder. \n"
        "Moreover, take care that the geometry files should be put in the "
        "'geometry' folder, as .pgm files.\n\n");
    exit(0);
  }

  // Allocate memory for the strings read from config file
  szFileName = malloc(strlen(problemName) + 15);
  sprintf(szFileName, "./configs/%s.dat", problemName);
  // create directory for .vtk files
  mkdir("paraviewFiles", 0777);
  dirName = malloc(strlen(problemName) + 17);
  sprintf(dirName, "./paraviewFiles/%s", problemName);
  mkdir(dirName, 0777);

  /* Input datafile*/
  geometry = (char *)malloc(50 * sizeof(char));
  problem = (char *)malloc(30 * sizeof(char));
  read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
                  &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,
                  &itermax, &eps, &dt_value, problem, geometry, &Pr, &TI, &T_h,
                  &T_c, &beta, &method, &toler, &special_input);

  // Method initialization
  int methods[2] = {3, 4};  // 5 time-stepping methods implemented
  int adaptive_method;
  search_ivalue(methods, 2, method, &adaptive_method);
  int k_stages;
  switch (method) {
    case 0:
      k_stages = 1;  // Explicit Euler
      break;
    case 1:
      k_stages = 2;  // RK2 - Heun method
      break;
    case 2:
      k_stages = 4;  // RK4
      break;
    case 3:
      k_stages = 4;  // ode23 (3 in, 4 outer)
      break;
    case 4:
      k_stages = 6;  // ode45 (4 in, 6 outer)
      break;
    default:
      printf("Error! Method %d not defined!\n", method);
      exit(0);
      break;
  }

  printf(
      "Welcome! Our time integration methods are:"
      "\n \t Method 0: Explicit Euler"
      "\n \t Method 1: Runge-Kutta 2 (Heun method)"
      "\n \t Method 2: Runge-Kutta 4"
      "\n \t Method 3: Bogacki-Shampine (adaptive time-step) (ode 23)"
      "\n \t Method 4: Fehlberg (adaptive time-step) (ode 45)\n"
      "Now you run Method %d\n\n",
      method);

  // names of .vtk files
  szProblem = malloc(strlen(dirName) + strlen(problem) + 15);
  sprintf(szProblem, "./%s/%s_%d", dirName, problem, method);

  // initialize U, V, P, T
  U = matrix(0, imax, 0, jmax + 1);
  V = matrix(0, imax + 1, 0, jmax);
  P = matrix(0, imax + 1, 0, jmax + 1);
  T = matrix(0, imax + 1, 0, jmax + 1);
  double **T_temp = matrix(1, imax, 1, jmax);
  F = matrix(0, imax, 1, jmax);
  G = matrix(1, imax, 0, jmax);
  RS = matrix(1, imax, 1, jmax);
  // Flag matrix
  Flag = imatrix(0, imax + 1, 0, jmax + 1);
  Geom = imatrix(1, imax, 1, jmax);
  init_flag(problem, geometry, imax, jmax, Flag, &num_fluid_cells, Geom);

  // allocate and initialize supplementary auxiliary matrices
  double **error_U = matrix(1, imax - 1, 1, jmax);
  double **error_V = matrix(1, imax, 1, jmax - 1);

  double **U_tilde = matrix(0, imax, 0, jmax + 1);
  double **V_tilde = matrix(0, imax + 1, 0, jmax);
  double ***Udot =
      (double ***)malloc(sizeof(double **) * k_stages);  // [U1, U2, U3, U4];
  double ***Vdot =
      (double ***)malloc(sizeof(double **) * k_stages);  // [U1, U2, U3, U4];
  for (int k = 0; k < k_stages; ++k) {
    Udot[k] = matrix(0, imax, 0, jmax + 1);  // U[k] - U1, U2, U3, U4, ...
    Vdot[k] = matrix(0, imax + 1, 0, jmax);  // V[k] - V1, V2, V3, V4, ...
    for (int i = 0; i <= imax; ++i)
      for (int j = 0; j <= jmax + 1; ++j) Udot[k][i][j] = 0;
    for (int i = 0; i <= imax + 1; ++i)
      for (int j = 0; j <= jmax; ++j) Vdot[k][i][j] = 0;
  }

  int total_sor_it = 0;
  int indicator = 1;  // 0 - we have the adaptive dt; 1 - we still search for it

  // main loop
  init_uvp(UI, VI, PI, TI, imax, jmax, U, V, P, T, F, G);
  index_dt_value = 0;  // counts number of writings in .vtk
  for (t = 0, n = 0; t < t_end; n++) {
    calculate_dt(Re, tau, &dt_ref, dx, dy, imax, jmax, U, V, Pr);
    if ((!adaptive_method) ||
        (n == 0)) {  // for adapting-time methods we use the optimal dt computed
                     // at step n-1
      dt = dt_ref;
    }
    do {
      boundaryvalues(imax, jmax, U, V, Flag, problem, T);
      specialboundaryvalues(imax, jmax, U, V, Flag, problem, T, T_h, T_c, t,
                            special_input);
      // printing vtk Files with timestep dt_value, at step n
      if (t >= index_dt_value * dt_value) {
        write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V,
                      P, T, Geom);
        index_dt_value += 1;
        printf("Writing vtkFile at t = %.3f, n = %d, using dt = %.8f\n", t, n,
               dt);
      }
      for (int k = 0; k < k_stages; ++k) {
        if ((method == 3) && (k == 0) && (n > 1)) {
          // Udot[0] at step (n) is Udot[3] from step (n-1) for ode23
          copy_matrices_dot(Udot, Vdot, 0, k_stages - 1, imax, jmax);
        } else {
          calculate_uv_tilde(dt, imax, jmax, U, V, U_tilde, V_tilde, Udot, Vdot,
                             method, k);
          boundaryvalues(imax, jmax, U_tilde, V_tilde, Flag, problem, T);
          specialboundaryvalues(imax, jmax, U_tilde, V_tilde, Flag, problem, T,
                                T_h, T_c, t, special_input);
          calculate_temp(Re, Pr, alpha, dt, dx, dy, imax, jmax, U_tilde,
                         V_tilde, T, T_temp, Flag);
          calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U_tilde,
                       V_tilde, F, G, Flag, T, beta);
          calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, Flag, U_tilde,
                       V_tilde);
          res = eps + 1;
          for (it = 0; it < itermax && res > eps; it++)
            sor(omg, dx, dy, imax, jmax, P, RS, &res, Flag, num_fluid_cells);
          total_sor_it += it;
          if (res > eps) {
            printf(
                "Warning: SOR iterations did not converge!\n\t n = %d: res = "
                "%.5f, "
                " eps = % .5f, it = % d, t = % .3f, omg = % .3f\n",
                n, res, eps, it, t, omg);
          }
          calculate_uv_stagek(dx, dy, imax, jmax, Udot, Vdot, F, G, P, Flag, k);
        }
      }
      calculate_uv(dt, imax, jmax, U, V, Udot, Vdot, Flag, method);
      if (adaptive_method) {  // calculate_optimal_dt
        calculate_adaptive_dt(&dt, dt_ref, Udot, Vdot, imax, jmax, toler,
                              &indicator, error_U, error_V, method);
      }
    } while (!indicator);
    t += dt;
    if (!(n % 10000)) printf("n = %d, t = %.9f\n", n, t);
  }
  write_vtkFile(szProblem, n + 1, xlength, ylength, imax, jmax, dx, dy, U, V, P,
                T, Geom);
  printf(
      "\nSuccessfully ended! The visualization files are in the folder %s.\n"
      "The name of the steady-state file is %s_%d.%d.vtk\n\n",
      dirName, problemName, method, n + 1);
  float avg_sor_steps = (float)total_sor_it;
  switch (method) {
    case 1:  // RK2
      avg_sor_steps /= 2;
      break;
    case 2:  // RK4
      avg_sor_steps /= 4;
      break;

    case 3:  // ode23
      avg_sor_steps /= 4;
      break;

    case 4:  // ode45
      avg_sor_steps /= 6;
      break;
  }

  printf(
      "Method %d:"
      "\n\t Total number of SOR iterations: %d"
      "\n\t Average number of SOR iterations per integration step: %.2f\n\n",
      method, total_sor_it, avg_sor_steps);

  // deallocate matrices
  free_matrix(U, 0, imax, 0, jmax + 1);
  free_matrix(V, 0, imax + 1, 0, jmax);
  free_matrix(P, 0, imax + 1, 0, jmax + 1);
  free_matrix(T, 0, imax + 1, 0, jmax + 1);
  free_matrix(T_temp, 1, imax, 1, jmax);
  free_matrix(F, 0, imax, 1, jmax);
  free_matrix(G, 1, imax, 0, jmax);
  free_matrix(RS, 1, imax, 1, jmax);
  free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);
  free_imatrix(Geom, 1, imax, 1, jmax);

  // deallocate supplementary auxiliary matrices
  free_matrix(error_U, 1, imax - 1, 1, jmax);
  free_matrix(error_V, 1, imax, 1, jmax - 1);

  free_matrix(U_tilde, 0, imax, 0, jmax + 1);
  free_matrix(V_tilde, 0, imax + 1, 0, jmax);

  for (int k = 0; k < k_stages; ++k) {
    free_matrix(Udot[k], 0, imax, 0, jmax + 1);  // U[k] - U1, U2, U3, U4, ...
    free_matrix(Vdot[k], 0, imax + 1, 0, jmax);  // V[k] - V1, V2, V3, V4, ...
  }
  free(Udot);
  free(Vdot);

  // deallocate strings
  // file names
  free(szFileName);
  free(dirName);
  free(szProblem);
  // paths
  free(problem);
  free(geometry);

  // time measurements
  clocks = clock() - clocks;
  double time_taken = clocks / CLOCKS_PER_SEC;  // in seconds
  printf("The simulation took %.1f seconds \n", time_taken);

  return -1;
}
