
#include <stdlib.h>
#include <math.h>
#include "thermal_integrands.hpp"


// equilibrium net baryon density (make separate nB type because may confuse the signs)
double nBeq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gauss laguerre (a = 1)
	return bn * pbar * exp(pbar) / (exp(Ebar-bn*alphaB)+a);
}

// equilibrium energy density
double Eeq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gauss laguerre (a = 2)
	return Ebar * exp(pbar) / (exp(Ebar-bn*alphaB)+a);
}

// equilibrium pressure
double Peq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gauss laguerre (a = 2)
	return pbar*pbar/Ebar * exp(pbar) / (exp(Ebar-bn*alphaB)+a);
}


////////////////////////////////////////////////////////////////////////


// for dalpha, dT, betaPi, betapi, betaV
double bn2_J10_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 1)
	return bn*bn * pbar * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}

double bn2_J11_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 1)
	return bn * bn * (pbar*pbar*pbar)/(Ebar*Ebar) * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}

double bn_J20_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 2)
	return bn * Ebar * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}

double J30_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 3)
	return Ebar*Ebar/pbar * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}

double J32_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 3)
	return (pbar*pbar*pbar / (Ebar*Ebar)) * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}


////////////////////////////////////////////////////////////////////////


// for linearized and modified particle densities
double neq_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	// gauss laguerre (a = 1)
	return pbar * exp(pbar) / (exp(Ebar-bn*alphaB)+a);
}

double bn_J10_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 1)
	return bn * pbar * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}


double J20_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 2)
	return Ebar * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}

double J21_integrand(double pbar, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;
	double qstat = exp(Ebar-bn*alphaB)+a;
	// gauss laguerre (a = 2)
	return pbar*pbar/Ebar * exp(pbar+Ebar-bn*alphaB)/(qstat*qstat);
}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                MODDED 3D THERMAL FUNCTION INTEGRANDS               ::
//                                                                    ::
//      Hydrodynamic moments of the modded thermal distribution       ::
//      function. Integrate with GaussMod3D. The p-coordinates        ::
//      are converted to the p_prime coordinates, which is easier     ::
//      for the most general momentum rescaling matrix A              ::
//                                                                    ::
//						 (10 components of T^munu)                    ::
//                                                                    ::
//              mod_E             mod_Txx          mod_Tyy            ::
//              mod_Tzz           mod_Txy          mod_Txz            ::
//              mod_Tyz           mod_Ttx          mod_Tty            ::
//              mod_Ttz                                               ::
// 																	  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double modE_integrand(double xphi, double costheta, double pbar, double ** A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

    // here {xphi, costheta, pbar} stand for the primed coordinates
	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta * costheta); // sintheta in right domain

	// convert p to p-prime coordinates with matrix transformation

	double p_unit[3] = {0.0,0.0,0.0}; // p_unit = p / pbar_prime = (A * p_prime) / pbar_prime

	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta}; // p_prime / pbar_prime in spherical coordinates

	// p components in terms of p_prime
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i]; // E(p) / pbar_prime
	energy_unit = sqrt(energy_unit);

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);

}


double modTxx_integrand(double xphi, double costheta, double pbar, double ** A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}

	//p_unit[0] = (pbar*unit_vec[0]+bn*V_alpha[0] - V_q[0]*sqrt((mbar*mbar + pbar*pbar*(1.0 - unit_vec[0]*unit_vec[0]))*(1.0 - V_q[0]*V_q[0]) + (pbar*unit_vec[0]+bn*V_alpha[0])*(pbar*unit_vec[0]+bn*V_alpha[0])))/(pbar*(1.0 - V_q[0]*V_q[0]));

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	//double detJ = 1.0 / (1.0 + px_unit * V_q[0] / energy_unit);

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;  // this could affect renormalization of particle density
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * px_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTyy_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_a[i] - V_p[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * py_unit * py_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTzz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * pz_unit * pz_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}



double modTxy_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double py_unit = p_unit[1];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * py_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTxz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double pz_unit = p_unit[2];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * pz_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTyz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];
	double pz_unit = p_unit[2];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * py_unit * pz_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTtx_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTty_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * py_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modTtz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * pz_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}



///////////////////////////////////////////////////////////////////////////////////


// for baryon diffusion


double modVx_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i]; // if include bulk would be bn * T/Tp * V_alpha[i]
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 1)
	return detJ * bn * pbar * px_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modVy_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 1)
	return detJ * bn * pbar * py_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}


double modVz_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;

	// gauss laguerre (a = 1)
	return detJ * bn * pbar * pz_unit / energy_unit * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}




// modified particle density

double modN_integrand(double xphi, double costheta, double pbar, double **A, double * V_q, double * V_alpha, int n, double mbar, double alphaB, int baryon, int sign)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double a = (double)sign;
	double bn = (double)baryon;

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k];
		}
		//p_unit[i] += (bn*V_alpha[i] - V_q[i]*Ebar)/pbar;
		p_unit[i] -= V_q[i]*Ebar/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double dalphaX = 0.0;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] / energy_unit;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*V_alpha[i] * pbar/ Ebar;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*V_alpha[i]/energy_unit;
	

	double detJ = 1.0 - pbar * unit_vec[0] * V_q[0] / Ebar;  // this is very makeshift

	// gauss laguerre (a = 1)
	return detJ * pbar * exp(pbar) / (exp(Ebar-bn*(alphaB+dalphaX)) + a);
}



