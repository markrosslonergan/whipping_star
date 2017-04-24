#include "SBNconfig.h"

SBNconfig::SBNconfig(){
	Ndet = 3;
	Nchan = 2;
	Nmode = 2;
	
	Chan[0] = 7;
	Chan[1] = 2;
	Bins[0] = 11;
	Bins[1] = 19;

	mname.push_back("nu");
	mname.push_back("nubar");

	dname.push_back("SBND");
	dname.push_back("uBooNE");
	dname.push_back("ICARUS");

	cname.push_back("elike");
	cname.push_back("mlike");

	scname.resize(Nchan);
	scname[0].push_back("fulloscnu"); 
	scname[0].push_back("fulloscnuebar"); 
	scname[0].push_back("intrinsic");
	scname[0].push_back("mismuon");
	scname[0].push_back("misphoton");
	scname[0].push_back("dirt");
	scname[0].push_back("cosmic");
	
	scname[1].push_back("intrinsic");
	scname[1].push_back("misncpion");


	std::string tempn;
	for(auto imode: mname){
		for(auto idet: dname){
			for(int ic = 0; ic < Nchan; ic++){
				for(int sc = 0; sc < Chan[ic]; sc++){
				tempn = imode +"_" +idet+"_"+cname[ic]+"_"+scname[ic][sc];
				std::cout<<tempn<<std::endl;
				fullnames.push_back(tempn);
				}	
			}
		}
	}

}
