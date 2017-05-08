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

class SBNosc : public SBNspec{
	protected:
	public:

	struct neutrinoModel workingModel;	
	int which_mode;
	double mStepSize;	

	SBNosc(const char *, std::string); //constructor
	SBNosc(const char *, std::string, neutrinoModel); //constructor

	std::vector< std::pair <double, int>> mass_splittings;	


	int OscillateThis();	
	std::vector<double> Oscillate();

	int load_model(neutrinoModel);	
	int calcMassSplittings();	

	//Oscillation mode 
	int setMode(int);
	void setAppMode();
	void setDisMode();
	void setBothMode();
	void setWierdMode();

};


};
#endif
