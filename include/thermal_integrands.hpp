
#include <stdlib.h>

#ifndef THERMAL_INTEGRANDS_H

#define THERMAL_INTEGRANDS_H


double nB_int(double pbar, double mbar, double T, double muB, double b, double a);
double E_int(double pbar, double mbar, double T, double muB, double b, double a);
double P_int(double pbar, double mbar, double T, double muB, double b, double a);


double Z10_int(double pbar, double mbar, double T, double muB, double b, double a);
double Z11_int(double pbar, double mbar, double T, double muB, double b, double a);
double N20_int(double pbar, double mbar, double T, double muB, double b, double a);
double J30_int(double pbar, double mbar, double T, double muB, double b, double a);
double J32_int(double pbar, double mbar, double T, double muB, double b, double a);

double neq_int(double pbar, double mbar, double T, double muB, double b, double a);
double N10_int(double pbar, double mbar, double T, double muB, double b, double a);
double J20_int(double pbar, double mbar, double T, double muB, double b, double a);
double J21_int(double pbar, double mbar, double T, double muB, double b, double a);


double modE_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTxx_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTyy_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTzz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTxy_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTxz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTyz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTtx_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTty_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modTtz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);

double modVx_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modVy_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);
double modVz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);

double modn_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a);


#endif