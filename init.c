#include "init.h"
#include "helper.h"

// 0 - Fluid; 1 - NO_SLIP; 2 - FREE_SLIP; 3 - OUTFLOW; 4 - INFLOW; 5 - COUPLING
// SET_FLAG(Flag_aux[i][j])
// Example:
//   Flag_aux[i][j] = 4 => (1 << ((4+1) % (4+1))) = 1 << 0 = 1 (FLUID)
#define SET_FLAG(x) (1 << ((x + 1) % (MASK_FLUID + 1)))  // set flag for BC / F

#define SET_E (SHIFT_EWSN + 3)  // shift 1 to left with 8
#define SET_W (SHIFT_EWSN + 2)  // shift 1 to left with 7
#define SET_S (SHIFT_EWSN + 1)  // shift 1 to left with 6
#define SET_N (SHIFT_EWSN)      // shift 1 to left with 5

int read_parameters(
    const char *szFileName,  // name of the file
    double *Re,              // Reynolds number
    double *UI,              // velocity x-direction
    double *VI,              // velocity y-direction
    double *PI,              // pressure
    double *GX,              // gravitation x-direction
    double *GY,              // gravitation y-direction
    double *t_end,           // end time
    double *xlength,         // length of the domain x-dir.
    double *ylength,         // length of the domain y-dir.
    double *dt,              // time step
    double *dx,              // length of a cell x-dir.
    double *dy,              // length of a cell y-dir.
    int *imax,               // number of cells x-direction
    int *jmax,               // number of cells y-direction
    double *alpha,           // uppwind differencing factor
    double *omg,             // relaxation factor
    double *tau,             // safety factor for time step
    int *itermax,      // max. number of iterations for pressure per time step
    double *eps,       // accuracy bound for pressure
    double *dt_value,  // time for output
    char *problem,     // scenario to be calculated (e.g., driven-cavity)
    char *geometry,    // file name/path to geometry.pgm
    double *Pr,        // Prandtl number
    double *TI,        // initial temperature
    double *T_h,       // hot temperature value
    double *T_c,       // cold temperature value
    double *beta,      // thermal expansion
    int *method,       // integration method
    double *toler,     // tolerance for time stepping
    int *special_input) {  // use or not sine inputs
  READ_DOUBLE(szFileName, *xlength);
  READ_DOUBLE(szFileName, *ylength);

  READ_DOUBLE(szFileName, *Re);
  READ_DOUBLE(szFileName, *t_end);
  READ_DOUBLE(szFileName, *dt);

  READ_INT(szFileName, *imax);
  READ_INT(szFileName, *jmax);

  READ_DOUBLE(szFileName, *omg);
  READ_DOUBLE(szFileName, *eps);
  READ_DOUBLE(szFileName, *tau);
  READ_DOUBLE(szFileName, *alpha);

  READ_INT(szFileName, *itermax);
  READ_DOUBLE(szFileName, *dt_value);

  READ_DOUBLE(szFileName, *UI);
  READ_DOUBLE(szFileName, *VI);
  READ_DOUBLE(szFileName, *GX);
  READ_DOUBLE(szFileName, *GY);
  READ_DOUBLE(szFileName, *PI);

  READ_STRING(szFileName, problem);
  READ_STRING(szFileName, geometry);

  READ_DOUBLE(szFileName, *Pr);
  READ_DOUBLE(szFileName, *beta);
  READ_DOUBLE(szFileName, *TI);
  READ_DOUBLE(szFileName, *T_h);
  READ_DOUBLE(szFileName, *T_c);

  READ_INT(szFileName, *method);
  READ_DOUBLE(szFileName, *toler);
  READ_INT(szFileName, *special_input);

  *dx = *xlength / (double)(*imax);
  *dy = *ylength / (double)(*jmax);

  return 1;
}

void init_uvp(double UI, double VI, double PI, double TI, int imax, int jmax,
              double **U, double **V, double **P, double **T, double **F,
              double **G) {
  init_matrix(U, 0, imax, 0, jmax + 1, UI);
  init_matrix(V, 0, imax + 1, 0, jmax, VI);
  init_matrix(P, 0, imax + 1, 0, jmax + 1, PI);
  init_matrix(T, 0, imax + 1, 0, jmax + 1, TI);
  init_matrix(F, 0, imax, 1, jmax, 0);
  init_matrix(G, 1, imax, 0, jmax, 0);
}

void init_flag(char *problem, char *geometry, int imax, int jmax, int **Flag,
               int *no_fluid_cells, int **Geom) {
  int **Flag_aux = read_pgm(geometry);
  int i, j;
  *no_fluid_cells = 0;
  // ==================== NORTH CELLS - (0:imax+1, jmax+1)
  j = jmax + 1;
  // NW - (0, jmax+1)
  i = 0;
  Flag[i][j] = ((MASK_FLUID == Flag_aux[i + 1][j]) << SET_E) |  // E
               ((MASK_FLUID == Flag_aux[i][j - 1]) << SET_S) |  // S
               SET_FLAG(Flag_aux[i][j]);                        // BC / F
  assert(!((Flag[i][j] >> SHIFT_EWSN) & SE));                   // !SE

  // N - (1:imax, jmax+1)
  for (i = 1; i <= imax; i++) {
    Flag[i][j] = ((MASK_FLUID == Flag_aux[i + 1][j]) << SET_E) |  // E
                 ((MASK_FLUID == Flag_aux[i - 1][j]) << SET_W) |  // W
                 ((MASK_FLUID == Flag_aux[i][j - 1]) << SET_S) |  // S
                 SET_FLAG(Flag_aux[i][j]);                        // BC / F
    assert((!(Flag[i][j] >> SHIFT_EWSN)) |     // no fluid around
           ((Flag[i][j] >> SHIFT_EWSN) & W) |  // fluid to W
           ((Flag[i][j] >> SHIFT_EWSN) & S) |  // fluid to S
           ((Flag[i][j] >> SHIFT_EWSN) & E));  // fluid to E
  }
  // NE - (imax+1, jmax+1)
  Flag[i][j] = ((MASK_FLUID == Flag_aux[i - 1][j]) << SET_W) |  // W
               ((MASK_FLUID == Flag_aux[i][j - 1]) << SET_S) |  // S
               SET_FLAG(Flag_aux[i][j]);                        // BC / F
  assert(!((Flag[i][j] >> SHIFT_EWSN) & SW));                   // EWSN != 0110

  for (j = jmax; j >= 1; j--) {
    // ============= WEST CELLS - (0,1:jmax
    i = 0;
    Flag[i][j] = ((MASK_FLUID == Flag_aux[i + 1][j]) << SET_E) |  // E
                 ((MASK_FLUID == Flag_aux[i][j - 1]) << SET_S) |  // S
                 ((MASK_FLUID == Flag_aux[i][j + 1]) << SET_N) |  // N
                 (SET_FLAG(Flag_aux[i][j]));
    assert((!(Flag[i][j] >> SHIFT_EWSN)) |     // EWSN = 0000
           ((Flag[i][j] >> SHIFT_EWSN) & E) |  // EWSN = 1000
           ((Flag[i][j] >> SHIFT_EWSN) & S) |  // EWSN = 0010
           ((Flag[i][j] >> SHIFT_EWSN) & N));  // EWSN = 0001

    // ==================== Interior cells - (1:imax, 1:jmax)
    for (i = 1; i <= imax; i++) {
      Flag[i][j] = ((MASK_FLUID == Flag_aux[i + 1][j]) << SET_E) |  // E
                   ((MASK_FLUID == Flag_aux[i - 1][j]) << SET_W) |  // W
                   ((MASK_FLUID == Flag_aux[i][j - 1]) << SET_S) |  // S
                   ((MASK_FLUID == Flag_aux[i][j + 1]) << SET_N) |  // N
                   SET_FLAG(Flag_aux[i][j]);                        // BC / F
      if (Flag[i][j] & 1) {
        (*no_fluid_cells)++;
        Geom[i][j] = 1;
      } else {
        Geom[i][j] = 0;
      }
      // fluid or (EW!=11 && SN != 11)
      assert((Flag[i][j] & 1) |                            // fluid or
             (((Flag[i][j] >> (SHIFT_EWSN + 2)) ^ EW) &&   // EW != 11 &&
              (((Flag[i][j] >> SHIFT_EWSN) & SN) ^ SN)));  // SN != 11
    }
    // ================= EAST CELLS - (imax+1,1:jmax)
    Flag[i][j] = ((MASK_FLUID == Flag_aux[i - 1][j]) << SET_W) |  // W
                 ((MASK_FLUID == Flag_aux[i][j - 1]) << SET_S) |  // S
                 ((MASK_FLUID == Flag_aux[i][j + 1]) << SET_N) |  // N
                 SET_FLAG(Flag_aux[i][j]);
    assert((!(Flag[i][j] >> SHIFT_EWSN)) |     // EWSN = 0000
           ((Flag[i][j] >> SHIFT_EWSN) & W) |  // EWSN = 0100
           ((Flag[i][j] >> SHIFT_EWSN) & S) |  // EWSN = 0010
           ((Flag[i][j] >> SHIFT_EWSN) & N));  // EWSN = 0001
  }
  // ==================== SOUTH CELLS  (0:imax+1, 0)
  // SW (0,0)
  i = 0;
  Flag[i][j] = ((MASK_FLUID == Flag_aux[i + 1][j]) << SET_E) |  // E
               ((MASK_FLUID == Flag_aux[i][j + 1]) << SET_N) |  // N
               SET_FLAG(Flag_aux[i][j]);                        // BC / F
  assert(!((Flag[i][j] >> SHIFT_EWSN) & NE));                   // EWSN != 1001

  // S (1:imax, 0)
  for (i = 1; i <= imax; i++) {
    Flag[i][j] = ((MASK_FLUID == Flag_aux[i + 1][j]) << SET_E) |  // E
                 ((MASK_FLUID == Flag_aux[i - 1][j]) << SET_W) |  // W
                 ((MASK_FLUID == Flag_aux[i][j + 1]) << SET_N) |  // N
                 SET_FLAG(Flag_aux[i][j]);                        // BC / F
    assert((!(Flag[i][j] >> SHIFT_EWSN)) |                        // EWSN = 0000
           ((Flag[i][j] >> SHIFT_EWSN) & E) |                     // EWSN = 1000
           ((Flag[i][j] >> SHIFT_EWSN) & W) |                     // EWSN = 0100
           ((Flag[i][j] >> SHIFT_EWSN) & N));                     // EWSN = 0001
  }
  // SE - (imax+1, 0)
  Flag[i][j] = ((MASK_FLUID == Flag_aux[i - 1][j]) << SET_W) |  // W
               ((MASK_FLUID == Flag_aux[i][j + 1]) << SET_N) |  // N
               SET_FLAG(Flag_aux[i][j]);                        // BC / F
  assert(!((Flag[i][j] >> SHIFT_EWSN) & NW));                   // EWSN = 0101

  // free Flag_aux, i.e. the 'pic' matrix from read_pgm()
  free_imatrix(Flag_aux, 0, imax + 1, 0, jmax + 1);
}
