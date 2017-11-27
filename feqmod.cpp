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
	// load hadron resonakce gas particles' mass, degeneracy, baryon, sign
	FILE *HRG;
	stringstream resonances;
	resonances << "pdg.dat";
	HRG = fopen(resonances.str().c_str(),"r");

	// pdg.dat contains (akti)mesons akd baryons, but not aktibaryons
	// so add aktibaryons makually
	int N_mesons, N_baryons, N_aktibaryons;
	// read 1st line: number of mesons
	fscanf(HRG, "%d", &N_mesons);
	// read 2nd line: number of baryons
	fscanf(HRG, "%d", &N_baryons);

	N_aktibaryons = N_baryons;

	// total number of resonances
	int N_resonances = N_mesons + N_baryons + N_aktibaryons;

	int particle_id;
	char name[20];
	double mass[N_resonances]; // file units = [GeV]
	double width;
	int degeneracy[N_resonances];
	int baryon[N_resonances], strakge, charm, bottom, isospin;
	double charge;
	int decays;

	int m = 0; // aktibaryon marker

	// load data of mesons+baryons
	for(int k = 0; k < N_mesons+N_baryons; k++)
	{
		fscanf(HRG, "%d %s %lf %lf %d %d %d %d %d %d %lf %d", &particle_id, name, &mass[k], &width, &degeneracy[k], &baryon[k], &strakge, &charm, &bottom, &isospin, &charge, &decays);

		if(baryon[k] == 1)
		{
			// makually add data of aktibaryons at end of array
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
		// convert resonakce masses to fm^-1
		mass[k] *= GEV_TO_INVERSE_FM;
	}

	fclose(HRG);






	//  Set up the pbar roots/weights for:   Eeq, PTeq, PLeq, piTxx, pixy, pixz, piyz, Qx, Qy, Qz (aA = 2)
	//                                       J32 (aR = 3)
	const int aT = 2; // associated Laguerre polynomials a = 2 for Eeq, PTeq, PLeq akd other T^munu components
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

	const int Ncheby = 32; // set to ak even value

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



	// set akgular roots-weights
	const int angle_pts = Ncheby+1; // number of akgular evaluation points
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


	// chemical freezeout T-muB data fit (Beccatini 2006 paper)
	// T = T0 - b0 * muB^2
	double T0 = 0.1675 * GEV_TO_INVERSE_FM;       // RHIC + SPS data
	double b0 = 0.1583 / GEV_TO_INVERSE_FM;

	// thermodynamic variables
	double T = 0.120 * GEV_TO_INVERSE_FM;         // temperature in fm^-1
	double muB = sqrt((T0 - T)/b0);               // baryon chemical potential in fm^-1

	double aB = muB/T;						      // baryon chemical potential over temperature

	cout << "muB = " << muB / GEV_TO_INVERSE_FM << " GeV" << endl;

	cout << "muB = " << muB << " fm^-1" << endl;
	cout << "T = " << T << " fm^-1" << endl;


	// equilibrium quaktities
	double factnB = pow(T,3) / (2.0*M_PI*M_PI);   // net baryon density
	double factE = pow(T,4) / (2.0*M_PI*M_PI);    // energy density
	double factP = pow(T,4) / (6.0*M_PI*M_PI);    // pressure

	double nB = 0.0;
	double E = 0.0;
	double P = 0.0;

	double dof, mbar, bk, ak;

	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];	mbar = mass[k]/T;
		bk = (double)baryon[k];			ak = (double)sign[k];

		if(baryon[k] != 0)
			nB += factnB * dof * Gauss1D(nB_int, pbar_rootN, pbar_weightN, gla_pts, mbar, T, muB, bk, ak);
		E += factE * dof * Gauss1D(E_int, pbar_rootT, pbar_weightT, gla_pts, mbar, T, muB, bk, ak);
		P += factP * dof * Gauss1D(P_int, pbar_rootT, pbar_weightT, gla_pts, mbar, T, muB, bk, ak);
	}

	cout << "Energy density = " << E << " fm^-4" << endl;
	cout << "Pressure = " << P << " fm^-4" << endl;
	cout << "nB = " << nB << " fm^-3" << endl;

	// Initialize viscous input
	double Pi = - 0.0 * P;
	double pixx = 0.5 * P;
	double piyy = - 0.5 * P;
	double pixy = 0.0 * P;
	double pixz = 0.0 * P;
	double piyz = 0.0 * P;
	double pizz = - (pixx + piyy);
	double Vx = 0.0 * nB;
	double Vy = 0.0 * nB;
	double Vz = 0.0 * nB;


	// Initialize stress energy tensor
	double Txx = P + pixx + Pi;
	double Tyy = P + piyy + Pi;
	double Tzz = P + pizz + Pi;
	double Txy = pixy;
	double Txz = pixz;
	double Tyz = piyz;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttz = 0.0;



	// Compute dmuB, dT, betaPi, betapi, betaV
	double factZ10 = pow(T,3) / (2.0*M_PI*M_PI);
	double factN20 = pow(T,4) / (2.0*M_PI*M_PI);
	double factJ30 = pow(T,5) / (2.0*M_PI*M_PI);
	double factJ32 = pow(T,5) / (30.0*M_PI*M_PI);
	double factZ11 = pow(T,3) / (6.0*M_PI*M_PI);

	double Z10 = 0.0;
	double N20 = 0.0;
	double J30 = 0.0;
	double J32 = 0.0;
	double Z11 = 0.0;

	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];	mbar = mass[k]/T;
		bk = (double)baryon[k];			ak = (double)sign[k];

		if(baryon[k] != 0)
		{
			Z10 += factZ10 * dof * Gauss1D(Z10_int, pbar_rootN, pbar_weightN, gla_pts, mbar, T, muB, bk, ak);
			N20 += factN20 * dof * Gauss1D(N20_int, pbar_rootT, pbar_weightT, gla_pts, mbar, T, muB, bk, ak);
			Z11 += factZ11 * dof * Gauss1D(Z11_int, pbar_rootN, pbar_weightN, gla_pts, mbar, T, muB, bk, ak);
		}
		J30 += factJ30 * dof * Gauss1D(J30_int, pbar_rootJ, pbar_weightJ, gla_pts, mbar, T, muB, bk, ak);
		J32 += factJ32 * dof * Gauss1D(J32_int, pbar_rootJ, pbar_weightJ, gla_pts, mbar, T, muB, bk, ak);
	}

	double G = T * ((E+P+muB*nB)*N20 - J30*nB - muB*(E+P)*Z10) / (J30*Z10 - N20*N20);
	double F = T * T * (N20*nB - (E+P)*Z10) / (J30*Z10 - N20*N20);

	double BPi = G*nB + F*(E+P-muB*nB)/T + 5.0*J32/(3.0*T);
	double BV = Z11 - nB*nB*T/(E+P);
	double Bpi = J32/T;


	double dT = Pi*F/BPi;
	double dmuB = Pi*G/BPi;

	double Tp = T + dT;
	double muBp = muB + dmuB;  // modified temperature akd baryon chemical potential


	cout << "dT = " << dT << endl;
	cout << "dmuB = " << dmuB << endl;
	cout << "betaPi = " << BPi << endl;
	cout << "dnB = " << dmuB*Z10/T + dT*N20/(T*T) - muB*Z10*dT/(T*T) + Pi*nB/BPi << endl;
	cout << "dE = " << dmuB*N20/T + dT*J30/(T*T) - muB*N20*dT/(T*T) + Pi*(E+P)/BPi << endl;
	cout << "dPi = " << dmuB*nB + dT*(E+P-muB*nB)/T + 5.0*Pi*J32/(3.0*T*BPi) - Pi << endl;
	cout << "Piterm = " << Pi/(3.0*BPi) << endl;

	// double G = (J30*nB - (E+P)*N20)/(N20*N20 - J30*Z10);
	// double F = T*T*((E+P)*Z10 - N20*nB)/(N20*N20 - J30*Z10);
	//double BPi = 5.0*J32/(3.0*T) + G*nB*T + F*(E+P)/T;

	// double dT = Pi*F/BPi;
	// double daB = Pi*G/BPi;
	// double Tp = T + dT;
	// double aBp = aB + daB;  // modified temperature akd baryon chemical potential



	// cout << "dT = " << dT << endl;
	// cout << "dalphaB = " << dalphaB << endl;
	// cout << "betaPi = " << betaPi << endl;
	// cout << "dnB = " << dalphaB*bk2_J10 + dT*bk_J20/(T*T) + Pi*nBeq/betaPi << endl;
	// cout << "dE = " << dalphaB*bk_J20 + dT*J30/(T*T) + Pi*(Eeq+Peq)/betaPi << endl;
	// cout << "dPi = " << dalphaB*nBeq*T + dT*(Eeq+Peq)/(T) + 5.0*Pi*J32/(3.0*T*betaPi) - Pi << endl;
	// cout << "Piterm = " << Piterm << endl;


	// momentum rescaling matrix A
	const int n = 3;
  	double **A = (double **) malloc(n * sizeof(double *));
  	for(int i = 0; i < n; i++) A[i] = (double *) malloc(n* sizeof(double));

    A[0][0] = 1.0+pixx/(2.0*Bpi)+Pi/(3.0*BPi); A[0][1] = pixy/(2.0*Bpi);				  A[0][2] = pixz/(2.0*Bpi);

    A[1][0] = pixy/(2.0*Bpi);				   A[1][1] = 1.0+piyy/(2.0*Bpi)+Pi/(3.0*BPi); A[1][2] = piyz/(2.0*Bpi);

  	A[2][0] = pixz/(2.0*Bpi);				   A[2][1] = piyz/(2.0*Bpi);				  A[2][2] = 1.0+pizz/(2.0*Bpi)+Pi/(3.0*BPi);

	// detA
	double detA = A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

	cout << "detA = " << detA << endl;


	// diffusion currents
	double Va[n] = {0.0,0.0,0.0};
	double Vq[n] = {0.0,0.0,0.0};

	if(muB != 0.0)
	{
		Va[0] = Vx / BV;
		Va[1] = Vy / BV;
		Va[2] = Vz / BV;

		Vq[0] = Vx * nB*T / (BV*(E+P));
		Vq[1] = Vy * nB*T / (BV*(E+P));
		Vq[2] = Vz * nB*T / (BV*(E+P));
	}

	// modified net baryon current
	double factN = detA * pow(Tp,3) / (8.0*M_PI*M_PI);
	double factV = detA * pow(Tp,3) / (8.0*M_PI*M_PI);
	double modnB = 0.0;
	double modVx = 0.0;
	double modVy = 0.0;
	double modVz = 0.0;

	// modified energy momentum tensor
	double factT = detA * pow(Tp,4) / (8.0*M_PI*M_PI);
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


	double renormalize[N_resonances];

	double neq, N10, J20, J21;

	double nlin, modn;

	double mbarp;

	// renormalization for individual species
	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];
		mbar = mass[k]/T;
		mbarp = mass[k]/Tp;
		bk = (double)baryon[k];
		ak = (double)sign[k];

		neq = factnB * dof * Gauss1D(neq_int, pbar_rootN, pbar_weightN, gla_pts, mbar, T, muB, bk, ak);
		if(baryon[k] != 0)
		{
			N10 = factnB * dof * Gauss1D(N10_int, pbar_rootN, pbar_weightN, gla_pts, mbar, T, muB, bk, ak);
		}
		else
		{
			N10 = 0.0;
		}
		J20 = factE * dof * Gauss1D(J20_int, pbar_rootT, pbar_weightT, gla_pts, mbar, T, muB, bk, ak);
		J21 = factP * dof * Gauss1D(J21_int, pbar_rootT, pbar_weightT, gla_pts, mbar, T, muB, bk, ak);


		// linearized akd modified particle density
		nlin = neq + Pi/(BPi*T)*(G*N10 + F*(J20-muB*N10)/T + J21);

		// neq + Pi/(BPi*T)*(F*J20/T + J21)

		//modn = factN * dof * GaussMod3D(modn_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);


		modn = 4.0 * dof * factN * Gauss1D(neq_int, pbar_rootN, pbar_weightN, gla_pts, mbarp, Tp, muBp, bk, ak);

		renormalize[k] = nlin / modn;

		if(nlin < 0.0)
		{
			printf("\nNegative linear density\n");
		}
		if(modn < 0.0)
		{
			printf("\nNegative mod density\n");
		}
	}



	double Z = 0.0;  // renormalization

	// Compute modification hydro outputs
	for(int k = 0; k < N_resonances; k++)
	{
		dof = (double)degeneracy[k];	mbar = mass[k]/T;
		mbarp = mass[k]/Tp;				Z = renormalize[k];
		bk = (double)baryon[k];			ak = (double)sign[k];

		if(baryon[k] != 0)
		{
			modnB += Z * bk * dof * factN * GaussMod3D(modn_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

			modVx += Z * dof * factV * GaussMod3D(modVx_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

			modVy += Z * dof * factV * GaussMod3D(modVy_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

			modVz += Z * dof * factV * GaussMod3D(modVz_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootN, pbar_weightN, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);
		}

		modE += Z * dof * factT * GaussMod3D(modE_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTxx += Z * dof * factT * GaussMod3D(modTxx_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTyy += Z * dof * factT * GaussMod3D(modTyy_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTzz += Z * dof * factT * GaussMod3D(modTzz_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTxy += Z * dof * factT * GaussMod3D(modTxy_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTxz += Z * dof * factT * GaussMod3D(modTxz_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTyz += Z * dof * factT * GaussMod3D(modTyz_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTtx += Z * dof * factT * GaussMod3D(modTtx_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTty += Z * dof * factT * GaussMod3D(modTty_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);

		modTtz += Z * dof * factT * GaussMod3D(modTtz_int, xphi_root, xphi_weight, costheta_root, costheta_weight, pbar_rootT, pbar_weightT, angle_pts, angle_pts, gla_pts, A, Vq, Va, n, mbarp, T, Tp, muBp, bk, ak);
	}


	double piTxx = 0.5 * (Txx - Tyy);
	double modpiTxx = 0.5 * (modTxx - modTyy);

	double bulk_in = (Txx + Tyy + Tzz)/3.0 - P;
	double modbulk = (modTxx + modTyy + modTzz)/3.0 - P;


	//modVx += bk2_J11 * V_alpha[0];
	//modTtx += nBeq * T * V_alpha[0];  // add truncated V_alpha correction


	// print input/output results for comparison

	printf("\n");

	cout << setprecision(5) << "nB       " << nB << "         " << setprecision(5) << (modnB / nB - 1.0) * 100 << " \% error" << "\n" << "modnB    " << setprecision(5) << modnB << endl;

	printf("\n");

	cout << setprecision(5) << "Vx       " << Vx << "         " << setprecision(5) << (modVx / Vx - 1.0) * 100 << " \% error" << "\n" << "modVx    " << setprecision(5) << modVx << endl;

	printf("\n");

	cout << setprecision(5) << "Vy       " << Vy << "         " << setprecision(5) << (modVy / Vy - 1.0) * 100 << " \% error" << "\n" << "modVy    " << setprecision(5) << modVy << endl;

	printf("\n");

	cout << setprecision(5) << "Vz       " << Vz << "         " << setprecision(5) << (modVz / Vz - 1.0) * 100 << " \% error" << "\n" << "modVz    " << setprecision(5) << modVz << endl;

	printf("\n");

	cout << setprecision(5) << "Ttt       " << E << "         " << setprecision(5) << (modE / E - 1.0) * 100 << " \% error" << "\n" << "modTtt    " << setprecision(5) << modE << endl;

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


	cout << setprecision(5) << "dE/E       " << 0.0 << "\n" << "dEmod/Eeq    " << setprecision(5) << (modE / E - 1.0) << endl;

	printf("\n");

	cout << setprecision(5) << "RPi/P      " << bulk_in / P  << "\n" << "RmodPi/Peq   " << setprecision(5) << modbulk / P << endl;

	printf("\n");

	cout << setprecision(5) << "Rpi/P      " << sqrt(2.0) * piTxx / P  << "\n" << "modRpi/Peq   " << setprecision(5) << sqrt(2.0) * modpiTxx / P << endl;

	cout << setprecision(5) << "dTxx/Txx      " << (modTxx - Txx) / Txx  << "\n" << endl;

	cout << setprecision(5) << "dTyy/Tyy      " << (modTyy - Tyy) / Tyy  << "\n" << endl;

	printf("\n");


	printf("\n\n");


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






