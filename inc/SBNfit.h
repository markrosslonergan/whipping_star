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
namespace sbn{


class SBNfit : public SBNchi {
	
	protected:

	std::vector<double> f_initial_values;
	std::vector<double> f_upper_values;
	std::vector<double> f_lower_values;
	std::vector<double> f_step_sizes;
	std::vector<int>    f_is_fixed;
	std::vector<std::string> f_param_names;

	std::string f_minimizer_mode;
	std::string f_minimizer_algo;


	SBNspec fOsc;

	public:
	SBNspec sigOsc;

	int num_params;
	std::vector<std::pair<std::string, int>> vec_scales;
	double * fX;

	double bf_chi;
	const double * bf_params;
	
	int num_func_calls;

	//using SBNchi::SBNchi;
	SBNfit(SBNspec bk, SBNspec sk,int npa);
	int load_signal(SBNspec);

	int initialize_norm(std::vector< std::pair<std::string, int>> );
	virtual double MinimizerCalcChi(const double * X);
	double Minimize();	

	//ROOT::Math::Minimizer* min ;     

	int setMethod(std::string, std::string);

	int setInitialValues(std::vector<double>);
	int setInitialValues(double);
	
	int setUpperValues(std::vector<double>);
	int setUpperValues(double);

	int setLowerValues(std::vector<double>);
	int setLowerValues(double);
	
	int setStepSizes(std::vector<double>);
	int setStepSizes(double);
	
	int setFixed(std::vector<int>);
	int setFixed(int);
	
	int setNames(std::vector<std::string>);

};


};
#endif
