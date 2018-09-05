#include "SBNdet.h"
using namespace sbn;

SBNdet::SBNdet(double h, double w, double l, double fh, double fw, double fl, double base){
	height = h;
	length = l;
	width = w;
	volume = h*l*w;

	f_height = fh;
	f_length = fl;
	f_width = fw;
	f_volume = fh*fl*fw;	

	baseline=base;

	dh = (h-fh)/2.0;
	dl = (l-fl)/2.0;
	dw = (w-fw)/2.0;

	proposal_modifier =1;

	fname ="../rootfiles/ntuple.ICARUS.root";		
	potmodifier = 1.0;
	identifier = DET_ICARUS;

}

SBNdet::SBNdet(int ident, bool ismu ){


	switch (ident)
	{
		case(DET_SBND):
			height = 400;
			length = 500;
			width = 2*200;
			volume = height*length*width;
			mass =112;

			if(ismu)
			{
				f_height = 370;
				f_length = 405;
				f_width =2*183.5 ;
				f_mass = 77.0;

				proposal_modifier = 0.8*5212690.0/(6.6*887966.0*MET2IMP) ;

			}
			else
			{

				f_height =350;
				f_length = 420;
				f_width = 2*173.5 ;
				f_mass = 71.4;

				proposal_modifier = 36798.0/(6.6*6931*MET2IMP) ;
			}


			f_volume = f_height*f_length*f_width;	

			baseline=110;

			dh = (height-f_height)/2.0;
			dl = (length-f_length)/2.0;
			dw = (width-f_width)/2.0;

			name = "SBND";
			fname ="../rootfiles/ntuple.SBND.root";	
			foscname ="../rootfiles/ntuple.SBND_fullosc.root";	
			fbarname ="../rootfiles/NUBAR_MODE/ntuple.SBND.root";	
			fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.SBND_fullosc.root";	

			potmodifier = 1.0;
			identifier = DET_SBND;

			break;
		case(DET_UBOONE):
			height = 233.0;
			length = 1037.0;
			width = 256.0;
			mass = 86.6;
			volume = height*length*width;

			if(ismu)
			{

				f_height = 203;
				f_length = 942;
				f_width = 226;
				f_mass = 60.5;
				proposal_modifier =1.1*0.907*mass/89.0*173302/(6.6*26182);
				//proposal_modifier =  5212690.0/(6.6*887966.0*MET2IMP) ;
			}
			else
			{
				f_height = 183.0;
				f_length = 957.0;
				f_width = 206.0;
				f_mass = 47.9;

				//proposal_modifier =36798.0/(6.6*6931*MET2IMP) ;
				proposal_modifier =0.9*  0.90*86.6/89.0*1469/(6.6*258.7);
				//proposal_modifier   =1469/1417.9;
			}

			f_volume = f_height*f_length*f_width;	

			baseline=470.0;


			dh = (height-f_height)/2.0;
			dl = (length-f_length)/2.0;
			dw = (width-f_width)/2.0;

			name = "uBooNE";
			fname ="../rootfiles/ntuple.uBooNE.root";	
			foscname ="../rootfiles/ntuple.uBooNE_fullosc.root";	
			fbarname ="../rootfiles/NUBAR_MODE/ntuple.uBooNE.root";	
			fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.uBooNE_fullosc.root";
			potmodifier = 2.0;
			identifier = DET_UBOONE;
			break;
		case(DET_ICARUS):

			height = 316.0;
			length = 1795;
			width = 4*150;
			volume = height*length*width;
			mass = 476;

			if(ismu)
			{

				f_height = 286;
				f_length = 1700;
				f_width = 4*133;
				f_mass = 363;
				//proposal_modifier =0.9*mass/476*173302*MET2IMP/(6.6*28182);
				proposal_modifier =1.1 *0.907*mass/476*173302/(6.6*26182);
				//proposal_modifier =  5212690.0/(6.6*887966.0*MET2IMP) ;
				//proposal_modifier = 173302*87*MET2IMP/(6.6*28182*89);



			}
			else
			{

				f_height = 266;
				f_length = 1715;
				f_width = 4*123.5;
				f_mass = 315;
				proposal_modifier =0.9*  0.90*86.6/89.0*1469/(6.6*258.7);
				//proposal_modifier = 1469*87*MET2IMP/(6.6*238.7*89);
			}




			f_volume = f_height*f_length*f_width;	

			baseline=600;

			dh = (height-f_height)/2.0;
			dl = (length-f_length)/2.0;
			dw = (width-f_width)/2.0;

			name = "ICARUS";
			fname ="../rootfiles/ntuple.ICARUS.root";	
			foscname ="../rootfiles/ntuple.ICARUS_fullosc.root";	
			fbarname ="../rootfiles/NUBAR_MODE/ntuple.ICARUS.root";	
			fbaroscname ="../rootfiles/NUBAR_MODE/ntuple.ICARUS_fullosc.root";	

			potmodifier = 1.0;
			identifier = DET_ICARUS;
			break;

		default:
			std::cout<<"#ERROR: SBNdet::SBNdet(int ident) basefun.c, ident not one of 1,2,3. "<<std::endl;
			exit(EXIT_FAILURE);
	}

}

bool SBNdet::is_active(double * pos){
	bool ans = false;

	if(0.0 <= pos[0] && pos[0]<= height && 0.0 <= pos[1] && pos[1]<= width && 0.0 <= pos[2] && pos[2]<= length  )
	{
		ans = true;
	}


	return ans;
}

bool SBNdet::is_fiducial(double * pos){
	bool ans = false;

	if(dh <= pos[0] && pos[0]<= height-dh && dw <= pos[1] && pos[1]<= width-dw && dl <= pos[2] && pos[2]<= length-dl  )
	{
		ans = true;
	}


	return ans;
}

int SBNdet::random_pos(TRandom3 * rangen, double * vec){

	vec[0]=dh + rangen->Uniform(height-dh);
	vec[1]=dw + rangen->Uniform(width-dw);
	vec[2]=dl + rangen->Uniform(length-dl);
return 0;
}


bool SBNdet::is_fully_contained(double *vertex,double * endpoint){
	bool ans = false;
	if( SBNdet::is_active(vertex)&& SBNdet::is_active(endpoint)){

		ans = true;
	}


	return ans;
}


double SBNdet::track_length_escape(double * in, double *out){
	// Returns in CM

	if(is_active(out) || !is_active(in) )
	{
		std::cout<<"#ERROR: SBNdet::track_length_escape. in array is outside OR out array is inside!"<<std::endl;
		return 0;
	}


	double a[3] = {in[0],in[1],in[2]};
	double m[3] = {out[0]-in[0],out[1]-in[1],out[2]-in[2]};
	double minlen= 1e12;
	int npos = 0;

	double t[6] = {-a[0]/m[0], (-a[0]+height)/m[0],-a[1]/m[1],(-a[1]+width)/m[1],-a[2]/m[2],(-a[2]+length)/m[2]};
	std::vector<double > tpos;

	for(int i=0;i<6;i++)
	{
		if(t[i]>=0)
		{
			npos++;
			tpos.push_back(t[i]);
		}
	}

	for(int j=0; j<npos; j++){

		double pt[3] = {a[0]+m[0]*tpos[j],a[1]+m[1]*tpos[j],a[2]+m[2]*tpos[j]};	
		double len = sqrt(pow(pt[0]-in[0],2)+pow(pt[1]-in[1],2)+pow(pt[2]-in[2],2));

		if(len<=minlen)
		{
			minlen=len;
		}

	}


	return minlen;


}


double SBNdet::osc_length(TRandom3 * rangen){

	return baseline - rangen->Uniform(50.0);


}

double SBNdet::sanford_wang(double p, double theta){
	double c9=1;
	double pb=7;
	double c1=220.7;
	double c2=1.080;
	double c3=1.0;
	double c4=1.978;
	double c5=1.32;
	double c6=5.572;
	double c7=0.0868;
	double c8=9.86;

	return c1*pow(p,c2)*(1.0-p/(pb-c9))*exp( -c3*pow(p,c4)/pow(pb,c5) -c6*theta*(pb-c7*pow(cos(theta),c8)));
	
}

double SBNdet::get_pion(TRandom3 * rangen ){

	//boy this needs to be fixed
	
	double tp = rangen->Uniform(5.0);
	double wei = rangen->Uniform(120);

	while(wei > sanford_wang(tp,150*0.001)){
	tp = rangen->Uniform(5.0);
	wei = rangen->Uniform(120);
	}	

	return tp;
	
}

double SBNdet::get_baseline(TRandom3*rangen, double ppion){
	double mpion=0.13957;
	double pionlifetime = 2.6033*pow(10,-8.0);
	double invsecond2meter = 3.336*pow(10,-11.0)*100;

	double tl = rangen->Uniform(50.0);
	double wei = rangen->Uniform(1.0);

	while(wei > exp(-mpion*tl*invsecond2meter/(pionlifetime*ppion)) ){
		 tl = rangen->Uniform(50.0);
		 wei = rangen->Uniform(1.0);
	}


	return tl;
}


double SBNdet::smear_energy(double En, double Percen, TRandom3 * rangen){
	double ans = 0;
	while(ans <=0){
		ans = rangen->Gaus(En,Percen*En/sqrt(En));
	}
	return ans;
} 

double SBNdet::smear_angle(double th, double an, TRandom3 * rangen){
	double ans = 0;
	ans = rangen->Gaus(th,an);
	return ans;
}
double massive_smear_energy(double En, double Percen, TRandom3 * rangen, double mass){
	double ans = 0;
	while(ans <=mass){
		ans = rangen->Gaus(En,Percen*En/sqrt(En));
	}
	return ans;
} 


double SBNdet::muon_track_length(double El){
	double	ans = 0.0;
	double list[390] ={0.000925764,0.022103,0.0648654,0.122844,0.192052,0.269917,0.354701,0.445186,0.540483,0.639929,0.743014,0.849338,0.958581,1.07048,1.18483,1.30145,1.42017,1.54089,1.66347,1.78783,1.91388,2.04155,2.17076,2.30146,2.43359,2.5671,2.70195,2.83811,2.97552,3.11416,3.25399,3.39499,3.53713,3.68038,3.82471,3.97012,4.11657,4.26405,4.41253,4.56201,4.71246,4.86386,5.01621,5.16949,5.32368,5.47878,5.63476,5.79162,5.94935,6.10794,6.26737,6.42763,6.58873,6.75063,6.91335,7.07687,7.24118,7.40627,7.57213,7.73877,7.90617,8.07432,8.24321,8.41285,8.58323,8.75433,8.92615,9.09869,9.27194,9.4459,9.62056,9.79591,9.97195,10.1487,10.3261,10.5042,10.6829,10.8623,11.0424,11.2232,11.4046,11.5866,11.7693,11.9526,12.1366,12.3212,12.5065,12.6923,12.8788,13.0659,13.2537,13.442,13.631,13.8205,14.0107,14.2015,14.3929,14.5848,14.7774,14.9705,15.1643,15.3586,15.5535,15.749,15.9451,16.1417,16.339,16.5368,16.7351,16.934,17.1335,17.3336,17.5342,17.7354,17.9371,18.1394,18.3422,18.5456,18.7495,18.954,19.159,19.3646,19.5707,19.7773,19.9845,20.1922,20.4004,20.6092,20.8185,21.0283,21.2387,21.4495,21.6609,21.8729,22.0853,22.2983,22.5117,22.7257,22.9402,23.1553,23.3708,23.5868,23.8034,24.0205,24.238,24.4561,24.6747,24.8937,25.1133,25.3334,25.5539,25.775,25.9966,26.2186,26.4412,26.6642,26.8878,27.1118,27.3363,27.5613,27.7868,28.0128,28.2392,28.4662,28.6936,28.9215,29.1499,29.3787,29.6081,29.8379,30.0682,30.299,30.5302,30.7619,30.9941,31.2268,31.4599,31.6936,31.9276,32.1622,32.3972,32.6327,32.8686,33.105,33.3419,33.5793,33.8171,34.0553,34.2941,34.5332,34.7729,35.013,35.2536,35.4946,35.7361};

	/*	TF1 f("CSDA", "CSDA_integrand(x)", 0.105, El);
		ROOT::Math::WrappedTF1 wf1(f);

	// Create the Integrator
	ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);


	// Set parameters of the integration
	ig.SetFunction(wf1, false);
	ig.SetRelTolerance(0.001);   
	ans = ig.Integral(0.105,El);
	//	std::cout << "integral result is " << ans <<std::endl;

	return ans*100; //for meters
	*/
	double y0 = 100*list[(int)floor((El-0.11)/0.02)];
	double y1 = 100*list[(int)floor((El-0.11)/0.02)+1];
	double x0 = 0.11+floor((El-0.11)/0.02)*0.02;
	double x1 = 0.11+(floor((El-0.11)/0.02)+1)*0.02;

	ans = y0+(y1-y0)*(El-x0)/(x1-x0);
	return ans;

}

double SBNdet::pion_track_length(double El){
	double	ans = 0.0;
	std::vector<double> list ={0.000925764,0.022103,0.0648654,0.122844,0.192052,0.269917,0.354701,0.445186,0.540483,0.639929,0.743014,0.849338,0.958581,1.07048,1.18483,1.30145,1.42017,1.54089,1.66347,1.78783,1.91388,2.04155,2.17076,2.30146,2.43359,2.5671,2.70195,2.83811,2.97552,3.11416,3.25399,3.39499,3.53713,3.68038,3.82471,3.97012,4.11657,4.26405,4.41253,4.56201,4.71246,4.86386,5.01621,5.16949,5.32368,5.47878,5.63476,5.79162,5.94935,6.10794,6.26737,6.42763,6.58873,6.75063,6.91335,7.07687,7.24118,7.40627,7.57213,7.73877,7.90617,8.07432,8.24321,8.41285,8.58323,8.75433,8.92615,9.09869,9.27194,9.4459,9.62056,9.79591,9.97195,10.1487,10.3261,10.5042,10.6829,10.8623,11.0424,11.2232,11.4046,11.5866,11.7693,11.9526,12.1366,12.3212,12.5065,12.6923,12.8788,13.0659,13.2537,13.442,13.631,13.8205,14.0107,14.2015,14.3929,14.5848,14.7774,14.9705,15.1643,15.3586,15.5535,15.749,15.9451,16.1417,16.339,16.5368,16.7351,16.934,17.1335,17.3336,17.5342,17.7354,17.9371,18.1394,18.3422,18.5456,18.7495,18.954,19.159,19.3646,19.5707,19.7773,19.9845,20.1922,20.4004,20.6092,20.8185,21.0283,21.2387,21.4495,21.6609,21.8729,22.0853,22.2983,22.5117,22.7257,22.9402,23.1553,23.3708,23.5868,23.8034,24.0205,24.238,24.4561,24.6747,24.8937,25.1133,25.3334,25.5539,25.775,25.9966,26.2186,26.4412,26.6642,26.8878,27.1118,27.3363,27.5613,27.7868,28.0128,28.2392,28.4662,28.6936,28.9215,29.1499,29.3787,29.6081,29.8379,30.0682,30.299,30.5302,30.7619,30.9941,31.2268,31.4599,31.6936,31.9276,32.1622,32.3972,32.6327,32.8686,33.105,33.3419,33.5793,33.8171,34.0553,34.2941,34.5332,34.7729,35.013,35.2536,35.4946,35.7361};

	/*	TF1 f("CSDA", "CSDA_integrand(x)", 0.105, El);
		ROOT::Math::WrappedTF1 wf1(f);

	// Create the Integrator
	ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);


	// Set parameters of the integration
	ig.SetFunction(wf1, false);
	ig.SetRelTolerance(0.001);   
	ans = ig.Integral(0.105,El);
	//	std::cout << "integral result is " << ans <<std::endl;

	return ans*100; //for meters
	*/

	if(floor((El-0.11)/0.02)+2 >= list.size()){
	
		ans = list.back()*100;

		return ans;
	}



	double y0 = 100*list.at((int)floor((El-0.11)/0.02));
	double y1 = 100*list.at((int)floor((El-0.11)/0.02)+1);
	double x0 = 0.11+floor((El-0.11)/0.02)*0.02;
	double x1 = 0.11+(floor((El-0.11)/0.02)+1)*0.02;

	ans = y0+(y1-y0)*(El-x0)/(x1-x0);
	


	return ans;

}



double SBNdet::photon_conversion_length(double ep, TRandom3 * r){

	double ans= 0.0;
	double att_len = 0;

	if(ep < 0.01)
	{
		att_len = 25.0;
	}
	else if (ep <= 1.0)
	{
		att_len = exp(0.321799*ep)*(-391.977+402.232/(pow(ep,0.00797561)));
	}
	else { 
		att_len = 15.0;
	}	


	ans = r->Exp(att_len);		
	//std::cout<<ep<<" "<<att_len<<" "<<ans/100.0<<std::endl;

	return ans;

}



double bethe_bloch(double betaV){

	if(betaV>1 || betaV <0){std::cout<<"ERROR: beta must be 0<b<1 in bethe!"<<std::endl;}

	double me = 0.0005;
	double mu = 0.105;
	double K  = 0.307075; // MeV mol^-1 cm^2
	double Z = 18;
	double A = 40;
	double I = 10*Z;
	double d = 0;
	double z = -1;
	double gam = 1/sqrt(1-betaV*betaV);
	double wmax = 2*me*betaV*betaV*gam*gam/(1+2*gam*me/mu+pow(me/mu,2));
	double ans = (1/1000.0)*fabs(K*z*z*Z/(A*betaV*betaV)*(0.5*log(2*me*pow(betaV*gam,2)*wmax/(I*I))-betaV*betaV-d/2));

	return ans ;

}


double CSDA_integrand(double Emu){
	double betaV = sqrt(Emu*Emu-0.105*0.105)/Emu;
	double rho = 1.3954;
	double ans = 1.0/(bethe_bloch(betaV)*rho);
	return ans;

}

double pion_containment(double posX, double posY, double posZ, TRandom3 * r)
{



	return 0;
}

int SBNdet::get_endpoint(double *vertex,double track_L,double * pl,double *  endpoint){
	double pnorm = sqrt(pl[0]*pl[0]+pl[1]*pl[1]+pl[2]*pl[2]);


	endpoint[0] = vertex[0]+track_L*pl[0]/pnorm;
	endpoint[1] = vertex[1]+track_L*pl[1]/pnorm;
	endpoint[2] = vertex[2]+track_L*pl[2]/pnorm;


	return 1;	
}


