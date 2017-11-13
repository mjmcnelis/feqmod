
#include <stdlib.h>

#ifndef GAUSS_INTEGRATION_H

#define GAUSS_INTEGRATION_H


double Gauss1D(double thermal_1D_integrand(double pbar, double mbar, double alphaB, int baryon, int sign), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar, double alphaB, int baryon, int sign);

double GaussMod3D(double mod_thermal_3D_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign), double * xphi_root, double * xphi_weight, double * costheta_root, double * costheta_weight, double * pbar_root, double * pbar_weight, int xphi_pts, int costheta_pts, int pbar_pts, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);


#endif
