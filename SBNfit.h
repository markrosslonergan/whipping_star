#ifndef SBNFIT_H_
#define SBNFIT_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

#include "SBNconfig.h"
#include "SBNspec.h"
#include "SBNchi.h"

#include <TH1D.h>
#include <string>
#include <TF1.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Math/GSLMinimizer.h"
#include "Math/GSLSimAnMinimizer.h"

class SBNfit : public SBNchi {
	
	protected:

	std::vector<double> fInitialValues;
	std::vector<double> fUpperValues;
	std::vector<double> fLowerValues;
	std::vector<double> fStepSizes;
	std::vector<int>  fFixed;
	std::vector<std::string> fNames;

	SBNspec fOsc;

	std::string fMinMode;
	std::string fMinAlgo;

	public:
	SBNspec sigOsc;

	int Npar;
	std::vector<std::pair<std::string, int>> vecScales;
	double * fX;

	double BFchi;
	const double * BFparam;
	int BFncalls;

	//using SBNchi::SBNchi;
	SBNfit(SBNspec bk, SBNspec sk,int npa);
	int load_signal(SBNspec);

	int initialize_norm(std::vector< std::pair<std::string, int>> );
	virtual double minim_calc_chi(const double * X);
	int init_minim();
	double minimize();	

	//ROOT::Math::Minimizer* min ;     

	int setInitialValues(std::vector<double>);
	int setUpperValues(std::vector<double>);
	int setLowerValues(std::vector<double>);
	int setStepSizes(std::vector<double>);
	int setFixed(std::vector<int>);
	int setNames(std::vector<std::string>);

};

#endif
