
#include <stdlib.h>

#ifndef THERMAL_INTEGRANDS_H

#define THERMAL_INTEGRANDS_H


double nBeq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double Eeq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double Peq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);


double bn2_J10_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double bn2_J11_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double bn_J20_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double J30_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double J32_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);

double neq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double bn_J10_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double J20_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);
double J21_integrand(double pbar, double mbar, double alphaB, int baryon, int sign);


double modE_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTxx_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTyy_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTzz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTxy_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTxz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTyz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTtx_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTty_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modTtz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);

double modVx_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modVy_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);
double modVz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);

double modN_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign);


#endif