#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flag,
                    char *problem, double **T);

void specialboundaryvalues(int imax, int jmax, double **U, double **V,
                           int **Flag, char *problem, double **T, double T_h,
                           double T_c, double t, int special_input);

#endif
