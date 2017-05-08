#ifndef SBNCHI_H_
#define SBNCHI_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNconfig.h"

#include "TMatrixT.h"
#include "TRandom3.h"
#include "TFile.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <ctime>
#include "params.h"

namespace sbn{


class SBNchi : public SBNconfig{
	
	int matrix_size ;
	int matrix_size_c ;
	int bigMsize;
	int contMsize;

	public:
	SBNspec bkgSpec;

	double lastChi;


	SBNchi(SBNspec); 
	SBNchi(SBNspec,std::string); 
	
	int load_bkg(SBNspec);


	double CalcChi(SBNspec sigSpec);
	double CalcChi(std::vector<double> );

	std::vector<std::vector<double >> vMcI;
	std::vector<std::vector<double >> vMc;


	TMatrixT<double> sys_fill_direct(std::string, std::string);
	TMatrixT<double> sys_fill_direct();

	void fake_fill(TMatrixT <double>&  M);
	void stats_fill(TMatrixT <double>&  M, std::vector<double> diag);


	std::vector<std::vector<double >> to_vector(TMatrixT <double> McI);

	void collapse_layer1(TMatrixT <double> & M, TMatrixT <double> & Mc);
	void collapse_layer2(TMatrixT <double> & M, TMatrixT <double> & Mc);
	void collapse_layer3(TMatrixT <double> & M, TMatrixT <double> & Mc);


	// Todo:
	std::vector<double >  calc_signal_events(struct neutrinoModel &nuModel);


	//int init_minim();
	//double minim_CalcChi(const double * x);
	//double minimize(neutrinoModel newModel, double ipot, double ipotbar);

};


};
#endif
