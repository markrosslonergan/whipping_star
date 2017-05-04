#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"

#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNspec.h"
#include "SBNosc.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


/*************************************************************
 *************************************************************
 *		Define a working instance class as I'm sick
 *		of re-nitialising and doing chi^2 every time
 ************************************************************
 ************************************************************/



/*************************************************************
 *************************************************************
 *		BEGIN Main::sbnfit.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{


	/*	int bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti;
		int contMsize = (N_e_bins+N_m_bins)*N_dets*N_anti;
		
		TMatrixT <double> Msys(bigMsize,bigMsize);
		TMatrixT <double> MsysC(contMsize,contMsize);

		for(int i =0; i< Msys.GetNcols(); i++){
		for(int j =0; j< Msys.GetNcols(); j++){
		Msys(i,j) = 1;
		}}

		contract_signal2_anti(Msys,MsysC);
		
	for(int i =0; i< MsysC.GetNcols(); i++){
		for(int j =0; j< MsysC.GetNcols(); j++){
			std::cout<<i<<" "<<j<<" "<<MsysC(i,j)<<std::endl;
		}}
exit(EXIT_FAILURE);
*/









//gSystem->Load("libTree");


// Just a few flags to control program flow.
bool fit_flag = false;
bool verbose_flag = false;
bool test_flag=false;
bool bkg_flag= false;
bool sample_flag = false;
bool cov_flag = false;
bool sens_flag=false;
bool comb_flag = false;

bool both_flag = true;
bool dis_flag = false;
bool app_flag = false;
bool wierd_flag = false;
int which_channel = BOTH_ONLY;

bool unit_flag = false;
bool fraction_flag = false;
bool anti_flag = false;
int anti_mode = 0;
bool inject_flag = false;


double bfmn[3] = {0.398107,1.0,0};
double bfue[3] = {0.13,0.14,0};
double bfum[3] = {0.15,0.13,0};
double bfphi[3] = {0.0,0.0,0.0};



bool margin = false;


int plotmode = 1;

bool pot_flag = false;

int mode_flag = 0;

int parallel_split=0;

int  sens_num=1;
double dm41 = -1.0;

double in_dm = 0;
double in_ue4 = 0;
double in_um4=0;


double pot_num =0;
double pot_num_bar =0;

double inPhi45 = 0;

bool stat_only = false;
int dis_which = 1;
int num_ster = 0;
int index; 
int iarg = 0;
opterr=1;
/*************************************************************
 *************************************************************
 *		Command Line Argument Reading
 ************************************************************
 ************************************************************/

const struct option longopts[] = 
{
	{"fit", 		no_argument, 		0, 'F'},
	{"bkg",			no_argument,		0, 'B'},
	{"mode",		required_argument,	0, 'M'},
	{"test", 		no_argument, 		0, 'T'},
	{"mass",		required_argument,	0, 'm'},
	{"anti",		no_argument,		0, 'A'},
	{"ue4", 		required_argument,	0, 'e'},
	{"help",		no_argument, 		0, 'h'},
	{"um4"	,		required_argument,	0, 'u'},
	{"inject",	required_argument,		0, 'I'},
	{"verbose",		no_argument, 		0, 'v'},
	{"sensitivity",		required_argument,	0, 'S'},
	{"stat-only",		no_argument,		0, 'l'},
	{"sample",		no_argument,		0, 's'},
	{"cov",			no_argument, 		0, 'c'},
	{"dis",			no_argument,	 	0, 'd'},
	{"wierd",		no_argument,		0, 'w'},
	{"signal",		required_argument,	0, 'g'},
	{"both",		no_argument,		0, 'b'},
	{"unitary",		no_argument,		0, 'n'},
	{"num",			required_argument,	0, 'N'},
	{"fraction",		required_argument,	0, 'f'},
	{"pot",			required_argument,	0, 'p'},
	{"plotmode",		required_argument,	0, 'k'},
	{"margin",		no_argument,		0, 'r'},
	{"app",			no_argument,		0, 'a'},
	{"phi",			required_argument,	0, 'P'},
	{0,			no_argument, 		0,  0},
};


while(iarg != -1)
{
	iarg = getopt_long(argc,argv, "I:daP:lf:nuM:e:m:svp:hS:crFTABN:bwk:g:", longopts, &index);

	switch(iarg)
	{
		case 'r':
			margin=true;
			break;
		case 'g':
			{
			
			double tm[3] ={0,0,0};
			double te[3] ={0,0,0};
			double tu[3] ={0,0,0};
			double tp[3] ={0,0,0};

			std::string s = optarg;
			std::string delimiter = ":";
			int cnt=0;
			size_t pos = 0;
			std::string token;
			while ((pos = s.find(delimiter)) != std::string::npos) {
			    token = s.substr(0, pos);
			    //std::cout << token << std::endl;
			    s.erase(0, pos + delimiter.length());

			    if(cnt<3){
				double mm=pow(10,round(log10(fabs( atof(token.c_str())))/0.04)*0.04);	
				tm[cnt] =mm; 
			    }
 			    if(cnt>=3&&cnt<6){
				te[cnt-3] = atof(token.c_str());	
			    } 
 			    if(cnt>=6&&cnt<9){
				tu[cnt-6] = atof(token.c_str());	
			    } 
			    if(cnt>=9&&cnt<12){
				tp[cnt-9] = atof(token.c_str());	
			    } 
			   
			    cnt++;
			}
			//std::cout << s << std::endl;//last one remains in s
			tp[2]=atof(s.c_str());
					

			}
			break;
		case 'P':
			inPhi45 = strtof(optarg,NULL);
			break;
		case 'I':
			{			
			inject_flag = true;
			std::string s = optarg;
			std::string delimiter = ":";
			//std::cout<<s.find(delimiter)<<std::endl;
			std::string spot = s.substr(0, s.find(delimiter));	
			std::string spotbar = s.substr(s.find(delimiter)+1);

			pot_num= atof(spot.c_str());
			pot_num_bar = atof(spotbar.c_str());
			}
		

			break;
		case 'F':
			fit_flag = true;
			//mS = strtof(optarg,NULL);
			break;
		case 'B':
			bkg_flag = true;
			break;
		case 'M':
			mode_flag = strtof(optarg,NULL);
			break;
		case 'A':
			anti_flag = true;
			anti_mode = 1;
			break;
		case 'n':
			unit_flag = true;
			break;
		case 'T':
			test_flag = true;
			break;	
		case 'v':
			verbose_flag = true;
			break;
		case 'p':
			{
			pot_flag = true;
			std::string s = optarg;
			std::string delimiter = ":";
			//std::cout<<s.find(delimiter)<<std::endl;
			std::string spot = s.substr(0, s.find(delimiter));	
			std::string spotbar = s.substr(s.find(delimiter)+1);

			pot_num= atof(spot.c_str());
			pot_num_bar = atof(spotbar.c_str());
			}
			break;
		case 'k':
			plotmode = strtof(optarg,NULL);
			break;
		case 'N':
			num_ster = strtof(optarg,NULL);
			break;
		case 's':
			sample_flag = true;
			break;
		case 'l':
			stat_only = true;
			break;
		case 'd':
			dis_flag = true;
			//	dis_which =strtof(optarg,NULL); //obsolete!
			which_channel = DIS_ONLY;
			break;
		case 'a':
			app_flag = true;
			which_channel = APP_ONLY;
			break;
		case 'w':
			wierd_flag = true;
			which_channel = WIERD_ONLY;
			break;
		case 'b':
			both_flag = true;
			which_channel = BOTH_ONLY;
			break;
		case 'm':
			in_dm  = strtof(optarg,NULL);
			break;
		case 'f':
			fraction_flag = true;
			parallel_split  = strtof(optarg,NULL);
			break;
		case 'u':
			in_um4  = strtof(optarg,NULL);
			break;
		case 'e':
			in_ue4  = strtof(optarg,NULL);
			break;
		case 'c':
			cov_flag = true;
			break;
		case 'S':
			sens_flag = true;
			sens_num  = strtof(optarg,NULL);
			break;
		case '?':
		case 'h':
			std::cout<<"Allowed arguments:"<<std::endl;
			std::cout<<"\t-F\t--fit\t\tRun a single SBN fit. Used in conjuction with --mass, --ue4, --um4"<<std::endl;
			std::cout<<"\t-m\t--mass\t\tSet 3p1 dm41 value"<<std::endl;
			std::cout<<"\t-u\t--um4\t\tSet 3p1 Um4 value"<<std::endl;
			std::cout<<"\t-e\t--ue4\t\tSet 3p1 Ue4 value"<<std::endl;
			std::cout<<"\t-T\t--test\t\trun SBN test code"<<std::endl;
			std::cout<<"\t-S\t--sensitivity\t\trun a full sensitivity fit. Required argument, number of steriles. Run with -a -d"<<std::endl;
			std::cout<<"\t\t\t --sensitivity 1 --app --dis 1 --mode 1,  makes 3p1_both_em.dat"<<std::endl;
			std::cout<<"\t-f\t--fraction\t\tRun fraction analysis, takes 1 argument for 3+N"<<std::endl;
			std::cout<<"\t-d\t--app\t\tRun app only sensitivity"<<std::endl;
			std::cout<<"\t-a\t--dis\t\tRun dis only sensitivity, takes 1 arg 0 for e-dis 1 for mudis"<<std::endl;
			std::cout<<"\t-l\t--stat-only\t\t Run with no systematics, ony stats"<<std::endl;
			std::cout<<"\t-B\t--bkg\t\trun bkg generating test code"<<std::endl;
			std::cout<<"\t-h\t--help\t\tDisplays this help message"<<std::endl;
			std::cout<<"\t-s\t--sample\t\t generate all the sin and sin^2 samples"<<std::endl;
			std::cout<<"\t-d\t\t\tRequired Argument. Creates a sin and sin^2 frequency ntuples for a dmsq."<<std::endl;
			std::cout<<"\t-v\t\t\tVerbose run, mostly debugging"<<std::endl;	
			std::cout<<"\t-n\t\t\tUnitary run"<<std::endl;
			std::cout<<"\t-A\t--anti\t\t Anti neutrino running mode"<<std::endl;
			return 0;
	}

}


SBNspec bkgSpec("precomp/SBN_bkg_all", "sbn.xml");
bkgSpec.calcFullVector();
bkgSpec.compressVector();

SBNchi chi(bkgSpec);

/*SBNspec sig("precomp/SBN_bkg_all", "sbn.xml");

for(double x =0.5 ;x <= 1.5; x+=0.025){
	SBNspec loopSpec = sig;
	
	loopSpec.Scale("elike_misphoton", x);
	loopSpec.calcFullVector();
	loopSpec.compressVector();
	std::cout<<"scaling: "<<x<<" "<<"Chi^2: "<< chi.calc_chi(loopSpec)<<std::endl;
}

std::cout<<"Begining SBNosc study"<<std::endl;
*/

	bkgSpec.ScaleAll(0.5);

	for(double m = -2.00; m <=2.04; m=m+0.04){
	 for(double sins2 = log10(0.25) ; sins2 > -5; sins2 = sins2 - 0.2){
	 	double uei = 0.1;
		double umi = sqrt(pow(10,sins2))/(2*uei);
		
				SBNosc oscSig("precomp/SBN_bkg_all","sbn.xml");
				oscSig.setAppMode();
				oscSig.ScaleAll(0.5);		

				double imn[3] = {sqrt(pow(10,m)),0,0};
				double iue[3] = {umi,0,0};
				double ium[3] = {uei,0,0};
				double iph[3] = {0,0,0};
				neutrinoModel signalModel(imn,iue,ium,iph);
				signalModel.numsterile = 1;
				
				oscSig.load_model(signalModel);
				oscSig.oscillate();

				std::cout<<m<<" "<<sins2<<" "<<chi.calc_chi(oscSig)<<std::endl;

	 }
	}



}
