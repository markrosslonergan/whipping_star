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

#include "tinyxml/tinyxml.h"

namespace sbn{
// All declarations are within the namespace scope.

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
	
	std::vector<std::string> mode_names; 			
	std::vector<std::string> detector_names; 		
	std::vector<std::string> channel_names; 		
	std::vector<std::vector<std::string >> subchannel_names; 

	// vector Bools for turning on and off
	std::vector<bool> mode_bool; 
	std::vector<bool> detector_bool; 
	std::vector<bool> channel_bool; 
	std::vector<std::vector<bool >> subchannel_bool; 


	std::vector<std::vector<double> > bin_edges;
	std::vector<std::vector<double> > bin_widths;


	std::map <std::string, std::vector<int> > mapIndex;

	std::vector<std::string> fullnames;
	std::vector<int> used_bins; 

	SBNconfig(std::string);
	SBNconfig(){};
//	SBNConfig(const char *);
//	int printnum_detectors();

	int num_multisim;
	std::vector<std::string> branch_names;
};

}


#endif
