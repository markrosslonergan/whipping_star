#ifndef SBNCONFIG_H_
#define SBNCONFIG_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
//#include <libconfig.h++>

class SBNconfig {

	protected:
	int Ndet;
	int Nchan;
	int Nmode;
	int Chan[100];
	int Bins[100];

	public:
	
	std::vector<std::string> mname; //modename
	std::vector<std::string> dname; //detectorname
	std::vector<std::string> cname; //channel_name
	std::vector<std::vector<std::string >> scname; //subchannel_name
//	std::map<int, std::string> getPos;	

	std::vector<std::string> fullnames;

	SBNconfig();
//	SBNConfig(const char *);
//	int printNdet();

};




#endif
