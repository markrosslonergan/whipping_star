#include "prob.h"
#include <complex>
#include <vector>
#include <iostream>

using namespace sbn;
/*********************************************
*  Constructors for neutrinoModel
* *******************************************/

/******	NULL constructor ******/
neutrinoModel::neutrinoModel(){
	zero();
	difference();
}


/******	3+3 constructor ******/
neutrinoModel::neutrinoModel(double * mn, double * ue, double * um, double * ph){
	zero();
	for(int i = 0; i < 3; i ++){
			mNu[i] = mn[i]; Ue[i] = ue[i];
			Um[i] = um[i];  phi[i] =ph [i];
	}
	difference();
	if(mNu[2]==0&&mNu[1]==0&&mNu[0]!=0){
		numsterile =1 ;
	} else if(mNu[2]==0){
		numsterile = 2;
	} else if(mNu[1]!=0&&mNu[2]!=0&&mNu[0]!=0){
		numsterile = 3;
	}
}

/******	3+1 constructor ******/
neutrinoModel::neutrinoModel(double  mn, double  ue4, double  um4){
	zero();
	mNu[0] = mn;
	Ue[0]=ue4;
	Um[0]=um4;

	numsterile = 1;

	difference();
}




/*********************************************
*  	Other Models
* *******************************************/

void neutrinoModel::printall(){
	
	std::cout<<"m4: "<<mNu[0]<<" m5: "<<mNu[1]<<" m6: "<<mNu[2]<<std::endl;
	std::cout<<"Ue4: "<<Ue[0]<<" Ue5: "<<Ue[1]<<" Ue6: "<<Ue[2]<<std::endl;
	std::cout<<"Uu4: "<<Um[0]<<" Uu5: "<<Um[1]<<" Uu6: "<<Um[2]<<std::endl;
	std::cout<<"phi45: "<<phi[0]<<" phi46: "<<phi[1]<<" phi56: "<<phi[2]<<std::endl;
	std::cout<<"NumSterile: "<<numsterile<<std::endl;		
	
	std::cout<<"dm41sq: "<<dm41Sq<<" dm51sq: "<<dm51Sq<<" dm61sq: "<<dm61Sq<<std::endl;
	std::cout<<"log10 dm41sq: "<<log10(dm41Sq)<<" log10 dm51sq: "<<log10(dm51Sq)<<" log10 dm61sq: "<<log10(dm61Sq)<<std::endl;

	std::cout<<"dm454q: "<<dm54Sq<<" dm64sq: "<<dm64Sq<<" dm65sq: "<<dm65Sq<<std::endl;

}


void neutrinoModel::zero(){
	for(int i = 0; i < 3; i ++){
			mNu[i] = 0; Ue[i] = 0;
			Um[i] = 0;  phi[i] = 0;
	}
	UUee=0;
	UUmm=0;
	UUem=0;
	UUme=0;
	dm54Sq = 0;
	dm64Sq = 0;
	dm65Sq = 0;
	dm41Sq = 0;
	dm51Sq = 0;
	dm61Sq = 0;
	numsterile = 1;


}

void neutrinoModel::difference(){
		dm41Sq = pow(mNu[0],2);
		dm51Sq = pow(mNu[1],2);
		dm61Sq = pow(mNu[2],2);
		dm54Sq = dm51Sq - dm41Sq;
		dm64Sq = dm61Sq - dm41Sq;
		dm65Sq = dm61Sq - dm51Sq;
	}


double neutrinoModel::oscProbSin(double Ev, double L)
{
	if(dm41Sq==0 && dm51Sq==0 && dm61Sq==0)
	{ 
		return 1.0;
	} 
	else 
	{
		return sin(2*1.26711*dm41Sq*L/Ev);
	}

}

double neutrinoModel::oscProbSinSq(double Ev, double L)
{
	if(dm41Sq==0 && dm51Sq==0 &&dm61Sq==0)
	{ 
		return 1.0;
	} 
	else 
	{
		return pow(sin(1.26711*dm41Sq*L/Ev),2.0);
	}


}

double neutrinoModel::oscProb(int a, int b, double Ev, double L){


	if(a == b)
	{
		return oscProb_dis(a, Ev, L);
	} 
	else
	{
		return oscProb_app(a,b, Ev, L);
	}
}




double neutrinoModel::oscProb_app(int a, int b, double Ev, double L){
	
	double nubarmod= 1.0;

	if(a < 0 && b < 0){
		nubarmod = -1.0;
	}

	double Ua4=0,Ua5=0,Ua6=0,Ub4=0,Ub5=0,Ub6=0;

	switch(a)
	{
		case 1:
		case -1:
			Ua4 = Ue[0];
			Ua5 = Ue[1];
			Ua6 = Ue[2];
			break;
		case 2:
		case -2:
			Ua4 = Um[0];
			Ua5 = Um[1];
			Ua6 = Um[2];
			break;
		case 3:
		case -3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}
	switch(b)
	{
		case 1:
		case -1:
			Ub4 = Ue[0];
			Ub5 = Ue[1];
			Ub6 = Ue[2];
			break;
		case 2:
		case -2:
			Ub4 = Um[0];
			Ub5 = Um[1];
			Ub6 = Um[2];
			break;
		case 3:
		case -3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}

	
	double phi54 = nubarmod*phi[0];
	double phi64 = nubarmod*phi[1];
	double phi65 = nubarmod*phi[2];
	
	
	double ans =0.0;

	ans =  -4.0*fabs(Ua5*Ub5*Ua4*Ub4)*cos(phi54)*pow(sin(1.27*dm54Sq*L/Ev),2.0);
	ans += -4.0*fabs(Ua6*Ub6*Ua4*Ub4)*cos(phi64)*pow(sin(1.27*dm64Sq*L/Ev),2.0);
	ans += -4.0*fabs(Ua5*Ub5*Ua6*Ub6)*cos(phi65)*pow(sin(1.27*dm65Sq*L/Ev),2.0);


	ans += 4.0*(fabs(Ua4*Ub4)+fabs(Ua5*Ub5)*cos(phi54)+fabs(Ua6*Ub6)*cos(phi64))*fabs(Ua4*Ub4)*pow(sin(1.27*dm41Sq*L/Ev),2.0);
	ans += 4.0*(fabs(Ua4*Ub4)*cos(phi54)+fabs(Ua5*Ub5)+fabs(Ua6*Ub6)*cos(phi65))*fabs(Ua5*Ub5)*pow(sin(1.27*dm51Sq*L/Ev),2.0);
	ans += 4.0*(fabs(Ua4*Ub4)*cos(phi64)+fabs(Ua5*Ub5)*cos(phi65)+fabs(Ua6*Ub6))*fabs(Ua6*Ub6)*pow(sin(1.27*dm61Sq*L/Ev),2.0);


	ans += (2.0*sin(phi54))*fabs(Ub5*Ua5*Ub4*Ua4)*sin(2.53*dm54Sq*L/Ev);
	ans += (2.0*sin(phi64))*fabs(Ub6*Ua6*Ub4*Ua4)*sin(2.53*dm64Sq*L/Ev);
	ans += (2.0*sin(phi65))*fabs(Ub6*Ua6*Ub5*Ua5)*sin(2.53*dm65Sq*L/Ev);


	ans += 2.0*(fabs(Ua5*Ub5)*sin(phi54)+fabs(Ua6*Ub6)*sin(phi64))*fabs(Ua4*Ub4)*sin(2.53*dm41Sq*L/Ev);
	ans += 2.0*(-fabs(Ua4*Ub4)*sin(phi54)+fabs(Ua6*Ub6)*sin(phi65))*fabs(Ua5*Ub5)*sin(2.53*dm51Sq*L/Ev);
	ans += 2.0*(-fabs(Ua4*Ub4)*sin(phi64)-fabs(Ua5*Ub5)*sin(phi65))*fabs(Ua6*Ub6)*sin(2.53*dm61Sq*L/Ev);

	if(ans <0 || ans >1){
		std::cout<<"#ERROR: Lets preserve probability shall we?\n#ERROR neutrinoModel::oscProb_dis @ prob.c\n#ERROR Prob: "<<ans<<" L: "<<L<<" Ev: "<<Ev<<std::endl;
	}

	return ans;
}


double neutrinoModel::oscProb_dis(int a, double Ev, double L){
	double ans =0.0;

	double Ua4 =0,Ua5=0,Ua6=0;

	switch(a)
	{
		case 1:
			Ua4 = Ue[0];
			Ua5 = Ue[1];
			Ua6 = Ue[2];
			break;
		case 2:
			Ua4 = Um[0];
			Ua5 = Um[1];
			Ua6 = Um[2];
			break;
		case 3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}

	ans = 1.0 - 4.0*pow(Ua4*Ua5*sin(1.27*dm54Sq*L/Ev),2.0);
	ans += -4.0*pow(Ua4*Ua6*sin(1.27*dm64Sq*L/Ev),2.0);
        ans += -4.0*pow(Ua5*Ua6*sin(1.27*dm65Sq*L/Ev),2.0);
	ans += -4.0*(1.0-Ua4*Ua4-Ua5*Ua5-Ua6*Ua6)*( Ua4*Ua4*pow(sin(1.27*dm41Sq*L/Ev),2.0)+ Ua5*Ua5*pow(1.27*dm51Sq*L/Ev,2.0)+Ua6*Ua6*pow(sin(1.27*dm61Sq*L/Ev),2.0) );

	if(ans <0 || ans >1){
		std::cout<<"#ERROR: Lets preserve probability shall we?\n#ERROR neutrinoModel::oscProb_dis @ prob.c\n#ERROR Prob: "<<ans<<" L: "<<L<<" Ev: "<<Ev<<std::endl;
	}


	return ans;
}

double neutrinoModel::oscAmp(int a, int b, int which_dm, int sqornot){


	if(a == b)
	{
		return oscAmp_dis(a, which_dm);
	} 
	else
	{
		return oscAmp_app(a,b, which_dm, sqornot);
	}
}





double neutrinoModel::oscAmp_app(int a, int b, int which_dm, int sqornot)
{

	double nubarmod= 1.0;

	if(a < 0 && b < 0){
		nubarmod = -1.0;
	}



	double Ua4=0,Ua5=0,Ua6=0,Ub4=0,Ub5=0,Ub6=0;

	switch(a)
	{
		case 1:
		case -1:
			Ua4 = Ue[0];
			Ua5 = Ue[1];
			Ua6 = Ue[2];
			break;
		case 2:
		case -2:
			Ua4 = Um[0];
			Ua5 = Um[1];
			Ua6 = Um[2];
			break;
		case 3:
		case -3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}
	switch(b)
	{
		case 1:
		case -1:
			Ub4 = Ue[0];
			Ub5 = Ue[1];
			Ub6 = Ue[2];
			break;
		case 2:
		case -2:
			Ub4 = Um[0];
			Ub5 = Um[1];
			Ub6 = Um[2];
			break;
		case 3:
		case -3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}

	
	double phi54 = nubarmod*phi[0];
	double phi64 = nubarmod*phi[1];
	double phi65 = nubarmod*phi[2];
	
	
	double ans =0.0;

	switch(which_dm)
	{
		case 41:
			if(sqornot==2)
			{
				ans = 4.0*(fabs(Ua4*Ub4)+fabs(Ua5*Ub5)*cos(phi54)+fabs(Ua6*Ub6)*cos(phi64))*fabs(Ua4*Ub4);
			}
			else if(sqornot ==1)
			{
				ans = 2.0*(fabs(Ua5*Ub5)*sin(phi54)+fabs(Ua6*Ub6)*sin(phi64))*fabs(Ua4*Ub4);
			}
			break;
		case 51:
			if(sqornot==2)
			{
				ans =4.0*(fabs(Ua4*Ub4)*cos(phi54)+fabs(Ua5*Ub5)+fabs(Ua6*Ub6)*cos(phi65))*fabs(Ua5*Ub5);
			}
			else if(sqornot ==1)
			{
				ans = 2.0*(-fabs(Ua4*Ub4)*sin(phi54)+fabs(Ua6*Ub6)*sin(phi65))*fabs(Ua5*Ub5);
			}

			break;
		case 61:
			if(sqornot==2)
			{
				ans =4.0*(fabs(Ua4*Ub4)*cos(phi64)+fabs(Ua5*Ub5)*cos(phi65)+fabs(Ua6*Ub6))*fabs(Ua6*Ub6);
			}
			else if(sqornot ==1)
			{
				ans = 2.0*(-fabs(Ua4*Ub4)*sin(phi64)-fabs(Ua5*Ub5)*sin(phi65))*fabs(Ua6*Ub6);
			}

			break;
		case 54:
			if(sqornot==2)
			{
				ans = -4.0*fabs(Ua5*Ub5*Ua4*Ub4)*cos(phi54);
			}
			else if(sqornot ==1)
			{
				ans = (2.0*sin(phi54))*fabs(Ub5*Ua5*Ub4*Ua4);
			}

			break;
		case 64:
			if(sqornot==2)


			{
				ans =-4.0*fabs(Ua6*Ub6*Ua4*Ub4)*cos(phi64);
			}
			else if(sqornot ==1)
			{
				ans = (2.0*sin(phi64))*fabs(Ub6*Ua6*Ub4*Ua4);
			}
 
			break;
		case 65:
			if(sqornot==2)
			{
				ans = -4.0*fabs(Ua5*Ub5*Ua6*Ub6)*cos(phi65);
			}
			else if(sqornot ==1)
			{
				ans = (2.0*sin(phi65))*fabs(Ub6*Ua6*Ub5*Ua5);
			}
 
			break;


	}

	return ans;


}


double neutrinoModel::oscAmp_dis(int a, int which_dm){


	double Ua4 =0,Ua5=0,Ua6=0;

	switch(a)
	{
		case 1:
		case -1:
			Ua4 = Ue[0];
			Ua5 = Ue[1];
			Ua6 = Ue[2];
			break;
		case 2:
		case -2:
			Ua4 = Um[0];
			Ua5 = Um[1];
			Ua6 = Um[2];
			break;
		case 3:
		case -3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}

	double ans = 0.0;

	switch(which_dm)
	{
		case 41:
			ans = -4.0*(1.0-Ua4*Ua4-Ua5*Ua5-Ua6*Ua6)*Ua4*Ua4;
			break;
		case 51:
		
			ans = -4.0*(1.0-Ua4*Ua4-Ua5*Ua5-Ua6*Ua6)*Ua5*Ua5;
			break;
		case 61:
		
			ans = -4.0*(1.0-Ua4*Ua4-Ua5*Ua5-Ua6*Ua6)*Ua6*Ua6;
			break;
		case 54:
			ans = -4.0*pow(Ua4*Ua5,2.0);
			break;
		case 64:
		 
			ans = -4.0*pow(Ua4*Ua6,2.0);
			break;
		case 65:
		 
			ans = -4.0*pow(Ua6*Ua5,2.0);
			break;


	}


	return ans;
}




/*********************************************
*  Other non- neutrinoModel functions, might be useful for arbitrary complex matrix U, non-unitarity?
* *******************************************/





double oscProb(int a, int b, double Ev, double L,  std::vector< std::vector < std::complex<double> > >    U, std::vector<std::vector<double> > dm ){

	double del = 0;
	if(a == b){ del = 1.0;}
	
	double ans = 0;

	//Allocate a 6x6 vector	
//	std::vector< std::vector < std::complex<double> > >  U(6,std::vector<std::complex<double>>(6) );
//	std::vector<std::vector < double > >  dm;


//	U[1][2] = std::complex<double>(10.0,1.0); 
//	std::cout<<real(U[1][2])<<" "<<arg(U[1][2])<<std::endl;


	for (int i=0; i<6; i++)
	{
		for(int j = 0; j < i; j++)
		{
		std::complex<double > temp = U[b][i]*conj(U[a][i])*conj(U[b][j])*U[a][j];
		ans += 4.0*real(temp)*pow(sin(1.27*dm[i][j]*L/Ev),2.0);
		ans += -2.0*imag(temp)*sin(2.53*dm[i][j]*L/Ev);		

		}
	}
	
	

	return del-ans;
}



double Pmue(double L, double E, double dm, double sin2)
{
	return sin2*pow( sin(1.27*dm*(L/1000.0)/E),2.0);
}


double Pmm(double L, double Ev, double Dm, double sinSq2thmm){

	double ans = 1.0-sinSq2thmm*pow( sin(Dm*(L/1000.0)*1.27/Ev), 2.0);

	if(ans <0 || ans >1){
		std::cout<<"frack. violatinf prob? "<<ans<<" l "<<L<<" Ev "<<Ev<<" Dm "<<Dm<< " sinSq2thmm "<<sinSq2thmm<<std::endl;
	}
	return ans;
}



