
#include <stdlib.h>
#include <math.h>
#include "thermal_integrands.hpp"


// equilibrium net baryon density (make separate nB type because may confuse the signs)
double nB_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	// gauss laguerre (a = 1)
	return b * pbar * exp(pbar) / (exp(Ebar-b*muB/T)+a);
}

// equilibrium energy density
double E_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	// gauss laguerre (a = 2)
	return Ebar * exp(pbar) / (exp(Ebar-b*muB/T)+a);
}

// equilibrium pressure
double P_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	// gauss laguerre (a = 2)
	return pbar*pbar/Ebar * exp(pbar) / (exp(Ebar-b*muB/T)+a);
}


////////////////////////////////////////////////////////////////////////


// for dmuB, dT, betaPi, betapi, betaV
double Z10_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 1)
	return b*b * pbar * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}

double Z11_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 1)
	return b * b * pbar*pbar*pbar / (Ebar*Ebar) * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}

double N20_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 2)
	return b * Ebar * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}

double J30_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 3)
	return Ebar*Ebar / pbar * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}

double J32_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 3)
	return pbar*pbar*pbar / (Ebar*Ebar) * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}


////////////////////////////////////////////////////////////////////////


// for linearized particle densities
double neq_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	// gauss laguerre (a = 1)
	return pbar * exp(pbar) / (exp(Ebar-b*muB/T)+a);
}

double N10_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 1)
	return b * pbar * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}


double J20_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 2)
	return Ebar * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
}

double J21_int(double pbar, double mbar, double T, double muB, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);
	double qstat = exp(Ebar-b*muB/T)+a;
	// gauss laguerre (a = 2)
	return pbar*pbar/Ebar * exp(pbar+Ebar-b*muB/T)/(qstat*qstat);
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


double modE_int(double xphi, double costheta, double pbar, double ** A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);


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
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i]; // E(p) / pbar_prime
	energy_unit = sqrt(energy_unit);

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	// gauss laguerre (a = 2)
	return detJ * pbar * energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);

}


double modTxx_int(double xphi, double costheta, double pbar, double ** A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	// I could take out the Va part:
	//p_unit[0] = (pbar*unit_vec[0]+bn*Va[0] + Vq[0]*sqrt((mbar*mbar + pbar*pbar*(1.0 - unit_vec[0]*unit_vec[0]))*(1.0 - Vq[0]*Vq[0]) + (pbar*unit_vec[0]+bn*Va[0])*(pbar*unit_vec[0]+bn*Va[0])))/(pbar*(1.0 - Vq[0]*Vq[0]));

	//p_unit[0] = (pbar*unit_vec[0] - Vq[0]*sqrt((mbar*mbar + pbar*pbar*(1.0 - unit_vec[0]*unit_vec[0]))*(1.0 - Vq[0]*Vq[0]) + (pbar*unit_vec[0])*(pbar*unit_vec[0])))/(pbar*(1.0 - Vq[0]*Vq[0]));

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	//for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	//double detJ = 1.0 / (1.0 + px_unit * Vq[0] / energy_unit);

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*Va[i] / energy_unit;  // this could affect renormalization of particle density
	//for(int i = 0; i < n; i++) dalphaX += p_unit[i]*Va[i] * pbar/ Ebar;
	//for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]/energy_unit;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * px_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTyy_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * py_unit * py_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTzz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * pz_unit * pz_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}



double modTxy_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double py_unit = p_unit[1];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * py_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTxz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];
	double pz_unit = p_unit[2];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * pz_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTyz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];
	double pz_unit = p_unit[2];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * py_unit * pz_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTtx_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * px_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTty_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * py_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modTtz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 2)
	return detJ * pbar * pz_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}



///////////////////////////////////////////////////////////////////////////////////


// for baryon diffusion


double modVx_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double px_unit = p_unit[0];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 1)
	return detJ * b * pbar * px_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modVy_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double py_unit = p_unit[1];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 1)
	return detJ * b * pbar * py_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}


double modVz_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};
	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;

	}
	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double pz_unit = p_unit[2];

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift

	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 1)
	return detJ * b * pbar * pz_unit / energy_unit * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}




// modified particle density
double modn_int(double xphi, double costheta, double pbar, double **A, double * Vq, double * Va, int n, double mbar, double T, double Tp, double muBp, double b, double a)
{
	double Ebar = sqrt(pbar*pbar + mbar*mbar);

	double cosphi = cos(M_PI*(1.0+xphi));
	double sinphi = sin(M_PI*(1.0+xphi));
	double sintheta = sqrt(1.0 - costheta*costheta);

	double p_unit[3] = {0.0,0.0,0.0};
	double unit_vec[3] = {cosphi*sintheta, sinphi*sintheta, costheta};

	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < n; k++)
		{
			p_unit[i] += A[i][k]*unit_vec[k] - 0.0*Vq[i]*Vq[k]*unit_vec[k];
		}
		//p_unit[i] += (b*(T/Tp)*Va[i] - Vq[i]*Ebar)/pbar;
		//p_unit[i] += b*(T/Tp)*Va[i]/pbar;
	}

	double energy_unit = (mbar*mbar) / (pbar*pbar);
	for(int i = 0; i < n; i++) energy_unit += p_unit[i] * p_unit[i];
	energy_unit = sqrt(energy_unit);

	double detJ = 1.0 - 0.0*Vq[0]*Vq[0] - pbar*unit_vec[0]*Vq[0]/Ebar;  // this is very makeshift


	detJ = 1.0;

	double heat = 0.0;
	//for(int i = 0; i < n; i++) heat += Vq[i]*unit_vec[i]*pbar;
	for(int i = 0; i < n; i++) heat += Vq[i]*p_unit[i]*pbar;

	double dalphaX = 0.0;
	for(int i = 0; i < n; i++) dalphaX += unit_vec[i]*Va[i]*pbar/Ebar;

	// gauss laguerre (a = 1)
	return detJ * pbar * exp(pbar) / (exp(Ebar+heat-b*(muBp/Tp+dalphaX)) + a);
}



