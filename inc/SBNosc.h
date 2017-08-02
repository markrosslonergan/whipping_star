#ifndef SBNOSC_H_
#define SBNOSC_H_

#include "SBNspec.h"
#include "prob.h"
#include "params.h"
#include <string>
#include <utility>

/**********************************************
 *	This is a less general 3+N osc spectrum
 * *******************************************/
namespace sbn{
// Note for this to work, you need to have precomputed spectra in data/precomp
// As this the classes used for our SBN paper, we assume that the precomputed frequencies 
// consist of 100 samples (in both Sin^2 and Sin)  from 0.01 eV^2 to 100 eV^2 
// They are labeled SBN_FREW_MASS_.root where  FREQ is either SIN or SINSQ and MASS is the log10 of the sterile ev^2, e.g -0.04, or 1.20 
// to 2 sig figures


class SBNosc : public SBNspec{
	protected:
	public:
	//this is the structure contains 3+N oscillation parameters (find in prob.h)
	struct neutrinoModel workingModel;	
	
	// which_mode to oscillate in  (APP, DIS, etc..) 
	int which_mode;
	double mStepSize;	//has to be 0.04 for now

	SBNosc(const char *, std::string); //constructor
	SBNosc(const char *, std::string, neutrinoModel); //constructor

	//find all the frequencies! Needs to know if a frequency corresponds to 41 51 54..etc.. so thats the int
	std::vector< std::pair <double, int>> mass_splittings;	

	//Oscillate the contained std::vector<TH1D> hists 
	int OscillateThis();	
	// Or just oscillate a copy and return the ompressed vector
	std::vector<double> Oscillate();
	std::vector<double> Oscillate(double);

	int load_model(neutrinoModel);	
	int calcMassSplittings();	

	//Oscillation mode 
	int setMode(int);
	void setAppMode();
	void setDisMode();
	void setBothMode();
	void setWierdMode();
	void setDisEMode();


};

};
#endif
