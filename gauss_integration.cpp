
#include <stdlib.h>
#include "gauss_integration.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                    1D and 3D GAUSS INTEGRATION                     ::
//                                                                    ::
//     Compute 1D thermal integrals over radial momentum              ::
//     bar using Gauss Laguerre quadrature. Compute 3D test thermal   ::
//	   and modded thermal integrals. The two angular integrals        ::
//     are computed with Clenshaw Curtis quadrature and the radial    ::
//     momentum bar with Gauss Laguerre quadrature.                   ::
//                                                                    ::
//                  Gauss1D                 GaussMod3D                ::
//                                                                    ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



// for thermal equilibrium functions (1D integral over pbar coordinate from 0 to infinity)
double Gauss1D(double thermal_1D_integrand(double pbar, double mbar, double alphaB, int baryon, int sign), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar, double alphaB, int baryon, int sign)
{
	double sum = 0.0;
	for(int k = 0; k < pbar_pts; k++) sum += pbar_weight[k] * thermal_1D_integrand(pbar_root[k], mbar, alphaB, baryon, sign);
	return sum;
}


// for modded thermal functions (3D integral now in spherical p_prime coordinates; i.e. {xphi, costheta, pbar} -> {xphi_prime, costheta_prime, pbar_prime})
double GaussMod3D(double mod_thermal_3D_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign), double * xphi_root, double * xphi_weight, double * costheta_root, double * costheta_weight, double * pbar_root, double * pbar_weight, int xphi_pts, int costheta_pts, int pbar_pts, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double sum = 0.0;
	for(int i = 0; i < xphi_pts; i++)
	{
		for(int j = 0; j < costheta_pts; j++)
		{
			for(int k = 0; k < pbar_pts; k++)
			{
				sum += xphi_weight[i] * costheta_weight[j] * pbar_weight[k] * mod_thermal_3D_integrand(xphi_root[i], costheta_root[j], pbar_root[k], A, V_q, V_alpha, n, mbar, alphaB, baryon, sign);
			}
		}
	}
	return sum;
}
