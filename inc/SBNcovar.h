#ifndef SBNCOVAR_H_
#define SBNCOVAR_H_

#include <cmath>
#include <vector>
#include <iostream>
#include "SBNspec.h"
#include "SBNconfig.h"

#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TNtuple.h"

#include "TROOT.h"
#include "TRint.h"
#include "TSystem.h"
#include "TStyle.h"

#include <map>
#include <ctime>
#include "params.h"

namespace sbn{


class SBNcovar : public SBNconfig{
	bool is_small_negative_eigenvalue;
	

	public:
	double tolerence_positivesemi;

	SBNspec spec_CV;	
	SBNspec spec_CV2;	
	SBNcovar(std::string rootfile, std::string xmlname);

	// a vector of num_multisim vectors, with a vector of subchannel*bin histograms in each	
	std::vector<SBNspec> multi_hists;
	
	int formCovarianceMatrix();
	TMatrixD full_covariance;
	TMatrixD frac_covariance;
	TMatrixD full_correlation;

	std::vector<std::vector<int> > vars_i;
	std::vector<std::vector<double> > vars_d;

	virtual bool eventSelection(int file);
	virtual int fillHistograms(int file, int uni, double wei);

};


};
#endif
