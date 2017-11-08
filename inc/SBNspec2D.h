#ifndef SBNSPEC2D_H_
#define SBNSPEC2D_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom3.h>


#include "params.h"
#include "SBNconfig.h"
#include "SBNspec.h"


namespace sbn{
//This is the basic class that holds all spectral information in whatever reco or true variable you have decided you want in the xml files.
// Inherits from SBNconfig as thats how its all configured/kept equal! :
class SBNspec2D : public SBNconfig{
	
	public:

	//Currently
	



};



}

#endif
