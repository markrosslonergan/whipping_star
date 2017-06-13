#ifndef SBNprob_H_
#define SBNprob_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <array>
#include <map>


#include <prob.h>


namespace sbn{

	struct SBNprob{
		public:	
			SBNprob(int);
			double t12,t13,t23,t14,t24,t34;
			double Dm21, Dm31,Dm41;
			double d13,d24,d34;

			double Vcc,Vnc;

			int dimension;
			double eV2GeV;	
			double km2invGeV;
			double rho;
			double degree;
			double Nneutrino;
			double conversion_parameter;

			complex_matrix hamiltonian;
			complex_matrix hamil_kin;
			complex_matrix potential;
			complex_matrix UtVU;
			complex_matrix U;
			complex_matrix Uconj;	




			double probabilityMatterExact(int, int ,double, double);

	};



}

#endif
