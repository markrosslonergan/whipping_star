#ifndef SBNCONFIG_H_
#define SBNCONFIG_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>

#include "TH1D.h"
#include "TFile.h"

#include "tinyxml.h"

namespace sbn{
// All declarations are within the namespace scope.

//Order is important, "xml order" means that we loop over all modes,detectors channels and subchannels, But finished each mode,det and channels before moving on.
// if we have modes m1 and m2, detectors d1 and d2, channels c1 and c2 and subchannels c1s1 c1s2 and c2s1 then the order of spectra would be..
// m1_d1_c1_c1s1
// m1_d1_c1_c1s2
// m1_d1_c2_c2s1
// m1_d2_c1_c1s1
// m1_d2_c1_c1s2
// m1_d2_c2_c2s1

// m2_d1_c1_c1s1
// m2_d1_c1_c1s2
// m2_d1_c2_c2s1
// m2_d2_c1_c1s1
// m2_d2_c1_c1s2
// m2_d2_c2_c2s1

class SBNconfig {

	protected:
	
	int num_detectors;
	int num_channels;
	int num_modes;
	
	//vectors of length num_channels
	std::vector<int> num_subchannels; 
	std::vector<int> num_bins;

	int num_bins_detector_block;
	int num_bins_mode_block;
	int num_bins_total;

	int num_bins_detector_block_compressed;
	int num_bins_mode_block_compressed;
	int num_bins_total_compressed;

	std::string correlation_matrix_rootfile;
	std::string correlation_matrix_name;

	public:
	std::string xmlname;	
	std::string multisim_name;	
	std::string data_path;


	//the xml names are the way we track which channels and subchannels we want to use later
	std::vector<std::string> mode_names; 			
	std::vector<std::string> detector_names; 		
	std::vector<std::string> channel_names; 		
	std::vector<std::vector<std::string >> subchannel_names; 


	// vector Bools for turning on and off certain modes/detectors..etc..
	std::vector<bool> mode_bool; 
	std::vector<bool> detector_bool; 
	std::vector<bool> channel_bool; 
	std::vector<std::vector<bool >> subchannel_bool; 

	//self explanatory
	std::vector<std::vector<double> > bin_edges;
	std::vector<std::vector<double> > bin_widths;

	//Given a string e.g "nu_ICARUS_elike_intrinisc" this map returns the index of the corresponding covariance matrix. Not really used.
	std::map <std::string, std::vector<int> > mapIndex;

	// Fullnames is kinda important, it contains all the concatanated names of all individual histograms that have been configured with the "use=1" attribute
	// The order is IMPORTANT its the same as defined in xml
	std::vector<std::string> fullnames;
	// If you have a large covariance matrix/spectrum (say containing nu and nubar mode) but only want to run with nu-mode (so you set use=0 in nubarmode) the vector used_bins contains all the bins that are actually in use. 
	std::vector<int> used_bins; 

	SBNconfig(std::string);
	SBNconfig(){};

	//For generating a covariance matrix from scratch, this contains the number of multisims (weights in weight vector) and their names.
	// For some reason I have decided that the first multisim, weight[0] must be the central value, =1
	int num_multisim;
	std::vector<std::string> branch_names;
};

}


#endif
