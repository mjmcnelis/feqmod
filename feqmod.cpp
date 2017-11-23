#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "gauss_integration.hpp"
#include "thermal_integrands.hpp"
#include "freearray.hpp"
//#include <gsl/gsl_sf.h>

#define GEV_TO_INVERSE_FM 5.067731

// temporary
const int a = 21;
const int gla_pts = 32;
double root_gla[a][gla_pts];
double weight_gla[a][gla_pts];

int load_gauss_laguerre_data()
{
  FILE *fp;

  stringstream laguerre_roots_weights;
  laguerre_roots_weights << "gla_roots_weights_" << gla_pts << "_points.txt";

  if((fp = fopen(laguerre_roots_weights.str().c_str(), "r")) == NULL)
  {
     return 1;
  }
  for(int i = 0; i < a; i++)
  {
   for(int j = 0; j < gla_pts; j++)
   {
      if(fscanf(fp, "%i %lf %lf", &i, &root_gla[i][j], &weight_gla[i][j])!= 3)
      	{
        	printf("error loading roots/weights of Gauss-Laguerre Quadradture at %d %d\n", i, j);
    		return 1;
    	}
   }
  }
  fclose(fp);
  return 0;
}




int main()
{
	// load hadron resonance gas particles' mass, degeneracy, baryon, sign
	FILE *HRG;
	stringstream resonances;
	resonances << "pdg.dat";
	HRG = fopen(resonances.str().c_str(),"r");

	// pdg.dat contains (anti)mesons and baryons, but not antibaryons
	// so add antibaryons manually
	int N_mesons, N_baryons, N_antibaryons;
	// read 1st line: number of mesons
	fscanf(HRG, "%d", &N_mesons);
	// read 2nd line: number of baryons
	fscanf(HRG, "%d", &N_baryons);

	N_antibaryons = N_baryons;

	// total number of resonances
	int N_resonances = N_mesons + N_baryons + N_antibaryons;

	int particle_id;
	char name[20];
	double mass[N_resonances]; // file units = [GeV]
	double width;
	int degeneracy[N_resonances];
	int baryon[N_resonances], strange, charm, bottom, isospin;
	double charge;
	int decays;

	int m = 0; // antibaryon marker

	// load data of mesons+baryons
	for(int k = 0; k < N_mesons+N_baryons; k++)
	{
		fscanf(HRG, "%d %s %lf %lf %d %d %d %d %d %d %lf %d", &particle_id, name, &mass[k], &width, &degeneracy[k], &baryon[k], &strange, &charm, &bottom, &isospin, &charge, &decays);

		if(baryon[k] == 1)
		{
			// manually add data of antibaryons at end of array
			mass[m+N_mesons+N_baryons] = mass[k];
			degeneracy[m+N_mesons+N_baryons] = degeneracy[k];
			baryon[m+N_mesons+N_baryons] = -1;
			m++;
		}
	}

	// sign array for bose/fermi distributions
	int sign[N_resonances];
	for(int k = 0; k < N_resonances; k++)
	{
		// degeneracy = 2*spin + 1
		if(degeneracy[k] % 2 == 0)
			sign[k] = 1;  // fermions
		else if(degeneracy[k] % 2 == 1)
			sign[k] = -1; // bosons
		// convert resonance masses to fm^-1
		mass[k] *= GEV_TO_INVERSE_FM;
	}

	fclose(HRG);






	//  Set up the pbar roots/weights for:   Eeq, PTeq, PLeq, piTxx, pixy, pixz, piyz, Qx, Qy, Qz (aA = 2)
	//                                       J32 (aR = 3)
	const int aT = 2; // associated Laguerre polynomials a = 2 for Eeq, PTeq, PLeq and other T^munu components
	const int aJ = 3; // associated Laguerre polynomials a = 3 for J32
	const int aN = 1; // associated Laguerre polynomials a = 1 for Neq
	double * pbar_rootT = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightT = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_rootJ = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightJ = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_rootN = (double *)malloc(gla_pts * sizeof(double));
	double * pbar_weightN = (double *)malloc(gla_pts * sizeof(double));


	// Load gauss laguerre roots-weights
	printf("Start loading gauss data...");
	int num_error;
	if((num_error = load_gauss_laguerre_data()) != 0)
	{
		fprintf(stderr, "Error loading gauss data (%d)!\n", num_error);
		return 1;
	}
	printf("done\n\n");



	// Set momentum bar roots-weights
	for(int i = 0; i < gla_pts; i++)
	{
		pbar_rootT[i] = root_gla[aT][i];
		pbar_weightT[i] = weight_gla[aT][i];
		pbar_rootJ[i] = root_gla[aJ][i];
		pbar_weightJ[i] = weight_gla[aJ][i];
		pbar_rootN[i] = root_gla[aN][i];
		pbar_weightN[i] = weight_gla[aN][i];
	}



	// Chebyshev root-weight generator (temporary)

	// phi = M_PI * (1 + xphi)        (variable substitution)

	const int Ncheby = 32; // set to an even value

	double * cheby_root = (double *)malloc((Ncheby+1) * sizeof(double));
	double * cheby_weight = (double *)malloc((Ncheby+1) * sizeof(double));

	double cheby_interval = M_PI / (double)Ncheby; // spacings of the cosine argument

	cheby_weight[0] = 1.0 / ((double)Ncheby * (double)Ncheby - 1.0);
	cheby_weight[Ncheby] = cheby_weight[0];

	for(int i = 0; i < Ncheby+1; i++) cheby_root[i] = cos((double)i * cheby_interval);

	double subw = 0.0;

	if(Ncheby % 2 == 0)
	{
		for(int i = 1; i < Ncheby; i++)
		{
			for(int j = 1; j < Ncheby/2; j++)
			{
				subw += 2.0 / (1.0 - 4.0 * (double)j * (double)j) * cos(2.0 * (double)i * (double)j * cheby_interval);
			}
			cheby_weight[i] = 2.0 / (double)Ncheby * (1.0 + subw + cos((double)i * M_PI) / (1.0 - (double)Ncheby * (double)Ncheby));
			subw = 0.0;
		}
	}
	else printf("Ncheby is not even!\n");



	// set angular roots-weights
	const int angle_pts = Ncheby+1; // number of angular evaluation points
	double * xphi_root = (double *)malloc(angle_pts * sizeof(double));
	double * xphi_weight = (double *)malloc(angle_pts * sizeof(double));
	double * costheta_root = (double *)malloc(angle_pts * sizeof(double));
	double * costheta_weight = (double *)malloc(angle_pts * sizeof(double));

	for(int i = 0; i < angle_pts; i++)
	{
		xphi_root[i] = cheby_root[i];
		xphi_weight[i] = cheby_weight[i];
		costheta_root[i] = cheby_root[i];
		costheta_weight[i] = cheby_weight[i];
	}



	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//                                                       ::
	//                   TESTING GENERAL CASE                ::
	//                                                       ::
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	// thermodynamic variables
	double T = 0.155 * GEV_TO_INVERSE_FM;         // input temperature in fm^-1
	double alphaB = 1.0;						  // baryon chemical potential over temperature



	// equilibrium quantities
	double factor_nBeq = pow(T,3) / (2.0*M_PI*M_PI);   // net baryon density
	double factor_Eeq = pow(T,4) / (2.0*M_PI*M_PI);    // energy density
	double factor_Peq = pow(T,4) / (6.0*M_PI*M_PI);    // pressure

	double nBeq = 0.0;
	double Eeq = 0.0;
	double Peq = 0.0;

	for(int k = 0; k < N_resonances; k++)
	{
		if(baryon[k] != 0)
			nBeq += factor_nBeq * (double)degeneracy[k] * Gauss1D(nBeq_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		Eeq += factor_Eeq * (double)degeneracy[k] * Gauss1D(Eeq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		Peq += factor_Peq * (double)degeneracy[k] * Gauss1D(Peq_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
	}



	// Initialize bulk pressure and shear stress tensor
	double Pi = 0.0 * Peq;
	double pixx = 0.0 * Peq;
	double piyy = 0.0 * Peq;
	double pixy = 0.0 * Peq;
	double pixz = 0.0 * Peq;
	double piyz = 0.0 * Peq;
	double pizz = - (pixx + piyy);


	// Initialize stress energy tensor
	double Txx = Peq + pixx + Pi;
	double Tyy = Peq + piyy + Pi;
	double Tzz = Peq + pizz + Pi;
	double Txy = pixy;
	double Txz = pixz;
	double Tyz = piyz;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttz = 0.0;


	// Initialize baryon diffusion current
	double Vx = 0.0 * nBeq;
	double Vy = 0.0 * nBeq;
	double Vz = 0.0 * nBeq;





	// functions for computing dalphaB, dT, betaPi, and betapi
	double factor_bn2_J10 = pow(T,3) / (2.0*M_PI*M_PI);   // b_n^2 * J10_n
	double factor_bn_J20 = pow(T,4) / (2.0*M_PI*M_PI);    // b_n * J20_n
	double factor_J30 = pow(T,5) / (2.0*M_PI*M_PI);       // J30
	double factor_J32 = pow(T,5) / (30.0*M_PI*M_PI);      // J32


	double bn2_J10 = 0.0;
	double bn_J20 = 0.0;
	double J30 = 0.0;
	double J32 = 0.0;


	for(int k = 0; k < N_resonances; k++)
	{
		if(baryon[k] != 0)
		{
			bn2_J10 += factor_bn2_J10 * (double)degeneracy[k] * Gauss1D(bn2_J10_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
			bn_J20 += factor_bn_J20 * (double)degeneracy[k] * Gauss1D(bn_J20_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		}
		J30 += factor_J30 * (double)degeneracy[k] * Gauss1D(J30_integrand, pbar_rootJ, pbar_weightJ, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		J32 += factor_J32 * (double)degeneracy[k] * Gauss1D(J32_integrand, pbar_rootJ, pbar_weightJ, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
	}


	double G = (J30*nBeq - (Eeq+Peq)*bn_J20)/(bn_J20*bn_J20 - J30*bn2_J10);
	double F = ((Eeq+Peq)*bn2_J10 - bn_J20*nBeq)/(bn_J20*bn_J20 - J30*bn2_J10);
	double betaPi = 5.0*J32/(3.0*T) + G*nBeq*T + F*(Eeq+Peq)*T;

	double dT = (Pi*F*T*T)/betaPi;
	double dalphaB = (Pi*G)/(betaPi);

	double Tp = T + dT;
	double alphaBp = alphaB + dalphaB;  // modified temperature and baryon chemical potential
	double Piterm = Pi / (3.0*betaPi);

	cout << "dT = " << dT << endl;
	cout << "dalphaB = " << dalphaB << endl;
	cout << "betaPi = " << betaPi << endl;
	cout << "dnB = " << dalphaB*bn2_J10 + dT*bn_J20/(T*T) + Pi*nBeq/betaPi << endl;
	cout << "dE = " << dalphaB*bn_J20 + dT*J30/(T*T) + Pi*(Eeq+Peq)/betaPi << endl;
	cout << "dPi = " << dalphaB*nBeq*T + dT*(Eeq+Peq)/(T) + 5.0*Pi*J32/(3.0*T*betaPi) - Pi << endl;
	cout << "Piterm = " << Piterm << endl;



	// term for baryon diffusion
	double factor_bn2_J11 = pow(T,3) / (6.0*M_PI*M_PI);   // b_n^2 * J11_n
	double bn2_J11 = 0.0;


	for(int k = 0; k < N_resonances; k++)
	{
		if(baryon[k] != 0)
			bn2_J11 += factor_bn2_J11 * (double)degeneracy[k] * Gauss1D(bn2_J11_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
	}

	double betaV = bn2_J11 - nBeq*nBeq*T/(Eeq+Peq);


	// momentum rescaling matrix A (p = A * p_prime)
	const int n = 3;
  	double **A = (double **) malloc(n * sizeof(double *));
  	for(int i = 0; i < n; i++) A[i] = (double *) malloc(n* sizeof(double));

    A[0][0] = 1.0 + (pixx*T)/(2.0*J32) + Piterm;  A[0][1] = (pixy*T)/(2.0*J32);				    A[0][2] = (pixz*T)/(2.0*J32);

    A[1][0] = (pixy*T)/(2.0*J32);				  A[1][1] = 1.0 + (piyy*T)/(2.0*J32) + Piterm;  A[1][2] = (piyz*T)/(2.0*J32);

  	A[2][0] = (pixz*T)/(2.0*J32);				  A[2][1] = (piyz*T)/(2.0*J32);				    A[2][2] = 1.0 + (pizz*T)/(2.0*J32) + Piterm;

	// detA
	double detA = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);


	// chemical diffusion current
	double V_alpha[n] = {0.0,0.0,0.0};

	// heat diffusion current
	double V_q[n] = {0.0,0.0,0.0};

	if(alphaB != 0.0)
	{
		V_alpha[0] = Vx / betaV;
		V_alpha[1] = Vy / betaV;
		V_alpha[2] = Vz / betaV;

		V_q[0] = Vx * nBeq*T / (betaV*(Eeq+Peq));
		V_q[1] = Vy * nBeq*T / (betaV*(Eeq+Peq));
		V_q[2] = Vz * nBeq*T / (betaV*(Eeq+Peq));
	}





	// modified net baryon current
	double fact_N = detA * pow(Tp,3) / (2.0*M_PI*M_PI);
	double fact_V = detA * pow(Tp,3) / (8.0*M_PI*M_PI);
	double modnB = 0.0;
	double modVx = 0.0;
	double modVy = 0.0;
	double modVz = 0.0;

	// modified energy momentum tensor
	double fact_T = detA * pow(Tp,4) / (8.0*M_PI*M_PI);
	double modE = 0.0;
	double modTxx = 0.0;
	double modTyy = 0.0;
	double modTzz = 0.0;
	double modTxy = 0.0;
	double modTxz = 0.0;
	double modTyz = 0.0;
	double modTtx = 0.0;
	double modTty = 0.0;
	double modTtz = 0.0;

	double fact_modN = detA * pow(Tp,3) / (8.0*M_PI*M_PI);



	double renormalize[N_resonances]; // renormalization for individual species

	double neq, bn_J10, J20, J21;
	double n_linear, n_prime;

	for(int k = 0; k < N_resonances; k++)
	{
		// need to construct these for the linearized particle densities of mesons and baryons
		neq = factor_nBeq * (double)degeneracy[k] * Gauss1D(neq_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		if(baryon[k] != 0)
			bn_J10 = factor_nBeq * (double)degeneracy[k] * Gauss1D(bn_J10_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		else
			bn_J10 = 0.0;
		J20 = factor_Eeq * (double)degeneracy[k] * Gauss1D(J20_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);
		J21 = factor_Peq * (double)degeneracy[k] * Gauss1D(J21_integrand, pbar_rootT, pbar_weightT, gla_pts, mass[k]/T, alphaB, baryon[k], sign[k]);

		// linearized particle density
		n_linear = neq + bn_J10*dalphaB + J20*dT/(T*T) + J21*Pi/(betaPi*T);
		//n_linear = neq;

		// modified particle density
		//n_prime = (double)degeneracy[k] * fact_N * Gauss1D(neq_integrand, pbar_rootN, pbar_weightN, gla_pts, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		n_prime = (double)degeneracy[k] * fact_modN * GaussMod3D(modN_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		renormalize[k] = n_linear / n_prime;
	}





	// Compute modification outputs
	for(int k = 0; k < N_resonances; k++)
	{
		// I expect modnB to have no errors
		if(baryon[k] != 0)
		{
			modnB += renormalize[k] * (double)baryon[k] * (double)degeneracy[k] * fact_modN * GaussMod3D(modN_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

			modVx += renormalize[k] * (double)degeneracy[k] * fact_V * GaussMod3D(modVx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

			modVy += renormalize[k] * (double)degeneracy[k] * fact_V * GaussMod3D(modVy_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

			modVz += renormalize[k] * (double)degeneracy[k] * fact_V * GaussMod3D(modVz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);
		}


		modE += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modE_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTxx += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTxx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTyy += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTyy_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTzz += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTzz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTxy += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTxy_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTxz += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTxz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTyz += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTyz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTtx += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTtx_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTty += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTty_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);

		modTtz += renormalize[k] * (double)degeneracy[k] * fact_T * GaussMod3D(modTtz_integrand, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, V_q, V_alpha, n, mass[k]/Tp, alphaBp, baryon[k], sign[k]);
	}


	double piTxx = 0.5 * (Txx - Tyy);
	double modpiTxx = 0.5 * (modTxx - modTyy);

	double bulk_in = (Txx + Tyy + Tzz)/3.0 - Peq;
	double modbulk = (modTxx + modTyy + modTzz)/3.0 - Peq;


	//modVx += bn2_J11 * V_alpha[0];
	//modTtx += nBeq * T * V_alpha[0];  // add truncated V_alpha correction


	// print input/output results for comparison

	printf("\n");

	cout << setprecision(5) << "nB       " << nBeq << "         " << setprecision(5) << (modnB / nBeq - 1.0) * 100 << " \% error" << "\n" << "modnB    " << setprecision(5) << modnB << endl;

	printf("\n");

	cout << setprecision(5) << "Vx       " << Vx << "         " << setprecision(5) << (modVx / Vx - 1.0) * 100 << " \% error" << "\n" << "modVx    " << setprecision(5) << modVx << endl;

	printf("\n");

	cout << setprecision(5) << "Vy       " << Vy << "         " << setprecision(5) << (modVy / Vy - 1.0) * 100 << " \% error" << "\n" << "modVy    " << setprecision(5) << modVy << endl;

	printf("\n");

	cout << setprecision(5) << "Vz       " << Vz << "         " << setprecision(5) << (modVz / Vz - 1.0) * 100 << " \% error" << "\n" << "modVz    " << setprecision(5) << modVz << endl;

	printf("\n");

	cout << setprecision(5) << "Ttt       " << Eeq << "         " << setprecision(5) << (modE / Eeq - 1.0) * 100 << " \% error" << "\n" << "modTtt    " << setprecision(5) << modE << endl;

	printf("\n");

	cout << setprecision(5) << "Txx       " << Txx << "        " << setprecision(5) << (modTxx / Txx - 1.0) * 100 << " \% error" << "\n" << "modTxx    " << setprecision(5) << modTxx << endl;

	printf("\n");

	cout << setprecision(5) << "Tyy       " << Tyy << "        " << setprecision(5) << (modTyy / Tyy - 1.0) * 100 << " \% error" << "\n" << "modTyy    " << setprecision(5) << modTyy << endl;

	printf("\n");

	cout << setprecision(5) << "Tzz       " << Tzz << "        " << setprecision(5) << (modTzz / Tzz - 1.0) * 100 << " \% error" << "\n" << "modTzz    " << setprecision(5) << modTzz << endl;

	printf("\n");

	cout << setprecision(5) << "Txy       " << Txy << "        " << setprecision(5) << (modTxy / Txy - 1.0) * 100 << " \% error" << "\n" << "modTxy    " << setprecision(5) << modTxy << endl;

	printf("\n");

	cout << setprecision(5) << "Txz       " << Txz << "        " << setprecision(5) << (modTxz / Txz - 1.0) * 100 << " \% error" << "\n" << "modTxz    " << setprecision(5) << modTxz << endl;

	printf("\n");

	cout << setprecision(5) << "Tyz       " << Tyz << "        " << setprecision(5) << (modTyz / Tyz - 1.0) * 100 << " \% error" << "\n" << "modTyz    " << setprecision(5) << modTyz << endl;

	printf("\n");

	cout << setprecision(5) << "Ttx       " << Ttx << "\n" << "modTtx    " << setprecision(5) << modTtx << endl;

	printf("\n");

	cout << setprecision(5) << "Tty       " << Tty << "\n" << "modTty    " << setprecision(5) << modTty << endl;

	printf("\n");

	cout << setprecision(5) << "Ttz       " << Ttz << "\n" << "modTtz    " << setprecision(5) << modTtz << endl;

	printf("\n");
	printf("\n");
	printf("Plots:\n\n");


	cout << setprecision(5) << "dE/Eeq       " << 0.0 << "\n" << "dEmod/Eeq    " << setprecision(5) << (modE / Eeq - 1.0) << endl;

	printf("\n");

	cout << setprecision(5) << "Txx/Peq      " << Txx / Peq  << "\n" << "Txxmod/Peq   " << setprecision(5) << modTxx / Peq << "\n" << "dTxx/Peq     " << setprecision(5) << (modTxx - Txx) / Peq << endl;

	printf("\n");

	cout << setprecision(5) << "Tyy/Peq      " << Tyy / Peq  << "\n" << "Tyymod/Peq   " << setprecision(5) << modTyy / Peq << "\n" << "dTyy/Peq     " << setprecision(5) << (modTyy - Tyy) / Peq << endl;

	printf("\n");

	cout << setprecision(5) << "Tzz/Peq      " << Tzz / Peq  << "\n" << "Tzzmod/Peq   " << setprecision(5) << modTzz / Peq << "\n" << "dTzz/Peq     " << setprecision(5) << (modTzz - Tzz) / Peq<< endl;

	printf("\n");

	cout << setprecision(5) << "Txy/Peq      " << pixy / Peq << "\n" << "Txymod/Peq   " << setprecision(5) << modTxy / Peq << "\n" << "dTxy/Peq     " << setprecision(5) << (modTxy - Txy) / Peq << endl;

	printf("\n");

	cout << setprecision(5) << "Txz/Peq      " << pixz / Peq << "\n" << "Txzmod/Peq   " << setprecision(5) << modTxz / Peq << "\n" << "dTxz/Peq     " << setprecision(5) << (modTxz - Txz) / Peq << endl;

	printf("\n");

	cout << setprecision(5) << "Tyz/Peq      " << piyz / Peq << "\n" << "Tyzmod/Peq   " << setprecision(5) << modTyz / Peq << "\n" << "dTyz/Peq     " << setprecision(5) << (modTyz - Tyz) / Peq << endl;

	printf("\n\n");

	// for(int i = 0; i < N_resonances; i++)
	// {
	// 	cout << mass[i]/GEV_TO_INVERSE_FM << "\t" << baryon[i] << "\t" << renormalize[i] << endl;
	// }

	printf("Freeing memory...");

	free(pbar_rootT);
	free(pbar_weightT);
	free(pbar_rootJ);
	free(pbar_weightJ);
	free(pbar_rootN);
	free(pbar_weightN);
	free(cheby_root);
	free(cheby_weight);
	free(xphi_root);
	free(xphi_weight);
	free(costheta_root);
	free(costheta_weight);

	free_2D(A,n);

	printf("done\n\n");

	return 0;
}






