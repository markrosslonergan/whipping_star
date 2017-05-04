#include "SBNconfig.h"
#include "TH1D.h"
#include "TFile.h"
#include <sstream>
#include <string>

SBNconfig::SBNconfig(std::string whichxml): xmlname(whichxml) {

	scname.resize(100);
	scBool.resize(100);
	char *end;

		TiXmlDocument doc( whichxml.c_str() );
		bool loadOkay = doc.LoadFile();
	    	TiXmlHandle hDoc(&doc);
	        TiXmlElement *pMode, *pDet, *pChan, *pCov;
		
		pMode = doc.FirstChildElement("mode");
		pDet =  doc.FirstChildElement("detector");
		pChan = doc.FirstChildElement("channel");
		pCov  = doc.FirstChildElement("covariance");

		while(pCov){
			CorrMatRoot = pCov->Attribute("file");
			CorrMatName = pCov->Attribute("name");
			pCov = pCov->NextSiblingElement("covariance");
		}
		while(pMode)
			{
			//std::cout<<"Mode: "<<pMode->Attribute("name")<<" "<<pMode->Attribute("use")<<std::endl;
			mname.push_back(pMode->Attribute("name"));	
			mBool.push_back(strtod(pMode->Attribute("use"),&end));	

			pMode = pMode->NextSiblingElement("mode");
		}

		pDet = doc.FirstChildElement("detector");
		while(pDet)
			{
			//std::cout<<"Detector: "<<pDet->Attribute("name")<<" "<<pDet->Attribute("use")<<std::endl;
			dname.push_back(pDet->Attribute("name"));
			dBool.push_back(strtod(pDet->Attribute("use"),&end));

			pDet = pDet->NextSiblingElement("detector");	
		}
		int nchan = 0;
		while(pChan)
			{
			//std::cout<<"Channel: "<<pChan->Attribute("name")<<" "<<pChan->Attribute("use")<<" Bins: "<<pChan->Attribute("numbins")<<std::endl;
			cname.push_back(pChan->Attribute("name"));
			cBool.push_back(strtod(pChan->Attribute("use"),&end));
			Bins.push_back(strtod(pChan->Attribute("numbins"), &end));		


			TiXmlElement *pBin = pChan->FirstChildElement("bins");
			//std::cout<<"Bin Edges: "<<pBin->Attribute("edges")<<" widths: "<<pBin->Attribute("widths")<<std::endl;
		        TiXmlElement *pSubChan;

			pSubChan = pChan->FirstChildElement("subchannel");
			int nsubchan=0;
			while(pSubChan){
				//std::cout<<"Subchannel: "<<pSubChan->Attribute("name")<<" use: "<<pSubChan->Attribute("use")<<std::endl;
				scname[nchan].push_back(pSubChan->Attribute("name"));
				scBool[nchan].push_back(strtod(pSubChan->Attribute("use"),&end));

				nsubchan++;
				pSubChan = pSubChan->NextSiblingElement("subchannel");	
			}
			Chan.push_back(nsubchan);

			std::stringstream iss(pBin->Attribute("edges"));
			std::stringstream pss(pBin->Attribute("widths"));

			double number;
			std::vector<double> binedge;
			std::vector<double> binwidth;
			while ( iss >> number ) binedge.push_back( number );
			while ( pss >> number ) binwidth.push_back( number );

			binEdges.push_back(binedge);
			binWidths.push_back(binwidth);
	
			nchan++;
			pChan = pChan->NextSiblingElement("channel");	
		}
		
		Nchan = cname.size();
		Nmode = mname.size();
		Ndet  = dname.size();


		if(false){

			std::cout<<" Covariance root path: "<<CorrMatRoot<<" and matrix name: "<<CorrMatName<<std::endl;
			std::cout<<"Modes: ";
			for(auto b: mname){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Modes Bools: ";
			for(auto b: mBool){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Dets: ";
			for(auto b: dname){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Dets Bools: ";
			for(auto b: dBool){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Bin Edges: ";
			for(auto b: binEdges[0]){
				std::cout<<b<<" ";
			}
			std::cout<<". Total#: "<<binEdges[0].size()<<std::endl;

			std::cout<<"Bin Widths: ";
			for(auto b: binWidths[0]){
				std::cout<<b<<" ";
			}
			std::cout<<". Total#: "<<binWidths[0].size()<<std::endl;
	}





	std::string tempn;
	int indexcount = 0;
	for(int im = 0; im < Nmode; im++){

		for(int id =0; id < Ndet; id++){

			for(int ic = 0; ic < Nchan; ic++){

				for(int sc = 0; sc < Chan[ic]; sc++){

					tempn = mname[im] +"_" +dname[id]+"_"+cname[ic]+"_"+scname[ic][sc];
				
			//std::cout<<tempn<<" from: "<<indexcount<<" to "<<indexcount+Bins[ic]-1<<" mBool: "<<mBool[im]<<" dBool: "<<dBool[id]<<" cBool: "<<cBool[ic]<<std::endl;
					// This is where you choose NOT to use some fields	
					if(mBool[im] && dBool[id] && cBool[ic] && scBool[ic][sc]){					
						fullnames.push_back(tempn);

						for(int k = indexcount; k < indexcount+Bins[ic]; k++){
								useBins.push_back(k);
							
						}
					}
					
					std::vector<int> tvec = {indexcount, indexcount+Bins[ic]-1}; 

					mapIndex[tempn] = 	tvec;				
					indexcount = indexcount + Bins[ic];
	

					/*				
					if(Bins[ic]==19){
						TH1D * tem =  ((TH1D*)f2.Get("dis_muon_sin")); 
						tem->SetName(tempn.c_str());
						tem->Write();
					}
					else if(Bins[ic]==11){

						TH1D * tem =  ((TH1D*)f2.Get("fullosc_nue_sin")); 
						tem->SetName(tempn.c_str());
						tem->Write();
					}
					*/
				

				}	
			}
		}
	}

	//For here on down everything is derivable, above is just until I actually get config working.
	Nmode = 0;
	for(bool i:mBool){	if(i) Nmode++;	}

	Ndet = 0;
	for(bool i:dBool){	if(i) Ndet++;	}
	
	Nchan = 0;
	for(bool i:cBool){	if(i) Nchan++;	}

	for(int i=0; i< Nchan; i++){
		Chan[i] = 0;
		for(bool j: scBool[i]){ if(j) Chan[i]++;}
	}



	Tdet = 0;
	TdetComp = 0;
	for(int i =0; i< Nchan; i++){
		Tdet += Chan[i]*Bins[i];
		TdetComp += Bins[i];	
	}		
	Tmode = Tdet*Ndet;
	TmodeComp = TdetComp*Ndet;
	Tall = Nmode*Tmode;
	TallComp = Nmode*TmodeComp;


	//std::cout<<"Ndet: "<<Ndet<<" Nmode: "<<Nmode<<" Nchan: "<<Nchan<<" Nsc[0]: "<<Chan[0]<<" Nsc[1]: "<<Chan[1]<<std::endl;
	//std::cout<<"Tdet: "<<Tdet<<" Tmode: "<<Tmode<<" Tall: "<<Tall<<" Tcomp: "<<TallComp<<std::endl;




/*
	f->Close();
	f2.Close();
*/
	





//std::cout<<"Map test: "<<mapIndex["nubar_SBND_elike_fulloscnu"][1]<<" "<<mapIndex["nubar_SBND_elike_fulloscnu"][0]<<std::endl;
}//end constructor


