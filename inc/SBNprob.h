#ifndef SBNprob_H_
#define SBNprob_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <map>
#include <globes/globes.h>   /* GLoBES library */


#include <prob.h>


namespace sbn{

	struct SBNprob{
		public:	
			SBNprob(int);
			SBNprob(int, std::vector<double>,std::vector<double>, std::vector<double>);

			int init();

			double t12,t13,t23,t14,t24,t34;
			double Dm21, Dm31,Dm41;
			double d13,d24,d34;

			double Vcc,Vnc;
			bool useMatterEffect;
			bool useNCMatterEffect;
			bool useAntiNeutrino;

			int dimension;
			double eV2GeV;	
			double km2invGeV;
			double rho;
			double degree;
			double Nneutrino;
			double conversion_parameter;

			double gaussian(double x, double mean, double sigma);

			complex_matrix hamiltonian;
			complex_matrix hamil_kin;
			complex_matrix potential;
			complex_matrix UtVU;
			complex_matrix U;
			complex_matrix Uconj;	
	
			int setMatterEffect(bool);
			int setAntiNeutrinoMode(bool);
			int setNCMatterEffect(bool);
			int setParameters( std::vector<double>,std::vector<double>, std::vector<double>);

			double probabilityVacuumExact(int a, int b ,double E, double L);
			
			double probabilityMatterExact(int a, int b ,double E, double L);

			double probabilityGlobes(int a, int b, int panti, double E, double L );

			double probabilityMatterExactSmear(int, int ,double, double, double p, double n);

			int plotProbabilityMatter(int a, int b, double Emin, double Emax, double L, double percen, double n, std::ofstream *filestream);

	};



}

#endif
