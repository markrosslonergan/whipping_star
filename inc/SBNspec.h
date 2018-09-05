#ifndef SBNSPEC_H_
#define SBNSPEC_H_
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include "SBNconfig.h"
#include <TH1D.h>
#include <TF1.h>
#include <string>
#include <TLine.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
//#include <TROOT.h>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <TStyle.h>
#include <ctime>
#include <TFile.h>
#include "params.h"
#include <TRandom3.h>
namespace sbn{
//This is the basic class that holds all spectral information in whatever reco or true variable you have decided you want in the xml files.
// Inherits from SBNconfig as thats how its all configured/kept equal! :


template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

class SBNspec : public SBNconfig{
	
	public:
	// this vector of hists contains all spectra used.
	// The order of filling is the same as the order defined in xml file!
	std::vector<TH1D > hist;
	std::map<std::string, int> map_hist;

	//This is the full concatanated vector (in xml order)	
	std::vector<double > fullVec;
	//This is the compessed vector, collapsing all subchannels down to a single channel
	std::vector<double > compVec;

	SBNspec(const char *, std::string); //Load in root file and config file
	SBNspec(const char *, std::string, bool); //Load in root file and config file, optional verbosity
	SBNspec(std::string); //Load in config file, create EMPTY hists 
	SBNspec(std::string, int); //Load in config file, create EMPTY hists, with optional numbering (e.g for multisims!) 
	SBNspec() {};

	int randomScale(); //mainly a debugging function, just randomly scales each hist by 0-2
	int poissonScale(); //Scales every histogram by a poissonian random number
	int poissonScale(TRandom3*); //Scales every histogram by a poissonian random number
	int setAsGaussian(double mean, double sigma, int n);
	int setAsFlat(double val);

	//Scales all vectors in hist by double
	int ScaleAll(double);
	//Scales all vectors whose xml names contains the string name. so you can scale all nu mode at one with Scale("nu");
	int Scale(std::string name, double val);
	//Same as above, but applies a bin-center dependant fucntion to the selected histograms
	int Scale(std::string name, TF1 *);
	//Same as above but normalises to value rather than scales
	int NormAll(double);
	int Norm(std::string name, double val);

	//need to store a history of the scales for oscillation purposes.  FIX THIS
	std::string scale_hist_name;
	double scale_hist_val;
	bool has_been_scaled;

	int Clear();

	//Recaculates the fullVec and compVec's
	int calcFullVector();
	int compressVector();

	double getTotalEvents();

	int getGlobalBinNumber(double invar, int which_hist);
	int getLocalBinNumber(double invar, int which_hist);
	

	//Just some debugging/checking
	int printFullVec();
	int printCompVec();
	//writeOut saves all to an externam rootfile, each individual subchannel and a stacked channel plot.
	int writeOut(std::string);
	
int writeSpec(std::string filename);


	int compareSBNspecs(SBNspec * compsec, std::string filename);
	//Addes two SBNspec together. must have same xml!
	int Add(SBNspec*);

};



}

#endif
