#ifndef SBNDET_H_
#define SBNDET_H_

#include <vector>
#include <string>
#include <TRandom3.h>
#include <cmath>
#include <iostream>
#include "params.h"


namespace sbn{

/*************************************************************
 *************************************************************
 *	TODO:
 *	    (4) overload detectors so I can just pass identifiers DONE
 ************************************************************
 ************************************************************/
	double massive_smear_energy(double En, double Percen, TRandom3 * rangen,double mass);
	
	double bethe(double beta);

	double CSDA_integrand(double Emu);
	double pion_containment(double posX, double posY, double posZ, TRandom3 * r);




class SBNdet {
	double dh, dw, dl;


	public:
	char const * name;
	bool mumode;
	double height, width, length;
	double f_height, f_width, f_length;
        double volume;
	double f_volume;
	double mass;
	double f_mass;	
	double baseline;
	char const * fname; // location of root ntuple containing data file		
	char const * foscname;	//location of full oscillated
	char const * fbarname; // location of root ntuple containing data file		
	char const * fbaroscname;	//location of full oscillated 
	double potmodifier;
	int identifier;
	double proposal_modifier;

	SBNdet (double h, double w, double l, double fh, double fw, double fl,double base);
	SBNdet (int identifier, bool ismue = false);

	double osc_length(TRandom3 * rangen);

	bool is_active(double * pos);
	bool is_fiducial(double * pos);
	
	int random_pos(TRandom3 * rangen, double * vec);

	double track_length_escape(double * inside, double * outside);
	bool is_fully_contained(double *vertex,double * endpoint);
	

	double smear_energy(double En, double Percen, TRandom3 *rangen);
	double smear_angle(double the, double ang, TRandom3 *rangen);

	double sanford_wang(double,double);
	double get_pion(TRandom3*);
	double get_baseline(TRandom3*, double);

	double muon_track_length(double El);

	double pion_track_length(double El);
	double photon_conversion_length(double ep, TRandom3 * r);

	int get_endpoint(double *vertex,double track_L,double * pl,double *  endpoint);
};

};

#endif
