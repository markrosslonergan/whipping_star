#ifndef SBNCONFIG_H_
#define SBNCONFIG_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <map>
//#include <libconfig.h++>

#include "tinyxml/tinyxml.h"


class SBNconfig {

	protected:
	
	int Ndet;
	int Nchan;
	int Nmode;
	std::vector<int> Chan;
	std::vector<int> Bins;

	int Tdet;
	int Tmode;
	int Tall;

	int TdetComp;
	int TmodeComp;
	int TallComp;


	std::string CorrMatRoot;
	std::string CorrMatName;
	public:
	std::string xmlname;	
	
	std::vector<std::string> mname; 		//modename
	std::vector<std::string> dname; 		//detectorname
	std::vector<std::string> cname; 		//channel_name
	std::vector<std::vector<std::string >> scname; //subchannel_name

	// vector Bools for turning on and off
	std::vector<bool> mBool; //modename
	std::vector<bool> dBool; //detectorname
	std::vector<bool> cBool; //channel_name
	std::vector<std::vector<bool >> scBool; //subchannel_name


	std::vector<std::vector<double> > binEdges;
	std::vector<std::vector<double> > binWidths;


	std::map <std::string, std::vector<int> > mapIndex;

	std::vector<std::string> fullnames;
	std::vector<int> useBins; 

	SBNconfig(std::string);
//	SBNConfig(const char *);
//	int printNdet();

};




#endif
