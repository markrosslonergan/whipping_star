#include "SBNconfig.h"
using namespace sbn;


SBNconfig::SBNconfig(std::string whichxml): xmlname(whichxml) {

	subchannel_names.resize(2);
	subchannel_bool.resize(2);
	char *end;

		TiXmlDocument doc( whichxml.c_str() );
		bool loadOkay = doc.LoadFile();
	    	TiXmlHandle hDoc(&doc);
	        TiXmlElement *pMode, *pDet, *pChan, *pCov, *pMC, *pData;
		
		pMode = doc.FirstChildElement("mode");
		pDet =  doc.FirstChildElement("detector");
		pChan = doc.FirstChildElement("channel");
		pCov  = doc.FirstChildElement("covariance");
		pMC   = doc.FirstChildElement("MCevents");
	
		pData   = doc.FirstChildElement("data");

		while(pData){
			data_path = pData->Attribute("path");
			pData = pData->NextSiblingElement("data");
		}



		while(pCov){
			correlation_matrix_rootfile = data_path + pCov->Attribute("file");
			correlation_matrix_name = pCov->Attribute("name");
			pCov = pCov->NextSiblingElement("covariance");
		}
		while(pMode)
			{
			//std::cout<<"Mode: "<<pMode->Attribute("name")<<" "<<pMode->Attribute("use")<<std::endl;
			mode_names.push_back(pMode->Attribute("name"));	
			mode_bool.push_back(strtod(pMode->Attribute("use"),&end));	

			pMode = pMode->NextSiblingElement("mode");
		}

		pDet = doc.FirstChildElement("detector");
		while(pDet)
			{
			//std::cout<<"Detector: "<<pDet->Attribute("name")<<" "<<pDet->Attribute("use")<<std::endl;
			detector_names.push_back(pDet->Attribute("name"));
			detector_bool.push_back(strtod(pDet->Attribute("use"),&end));

			pDet = pDet->NextSiblingElement("detector");	
		}
		int nchan = 0;
		while(pChan)
			{
			channel_names.push_back(pChan->Attribute("name"));
			channel_bool.push_back(strtod(pChan->Attribute("use"),&end));
			num_bins.push_back(strtod(pChan->Attribute("numbins"), &end));		


			TiXmlElement *pBin = pChan->FirstChildElement("bins");
		        TiXmlElement *pSubChan;

			pSubChan = pChan->FirstChildElement("subchannel");
			int nsubchan=0;
			while(pSubChan){
				//std::cout<<"Subchannel: "<<pSubnum_subchannels->Attribute("name")<<" use: "<<pSubnum_subchannels->Attribute("use")<<std::endl;
				subchannel_names[nchan].push_back(pSubChan->Attribute("name"));
				subchannel_bool[nchan].push_back(strtod(pSubChan->Attribute("use"),&end));

				nsubchan++;
				pSubChan = pSubChan->NextSiblingElement("subchannel");	
			}
			num_subchannels.push_back(nsubchan);

			std::stringstream iss(pBin->Attribute("edges"));
			std::stringstream pss(pBin->Attribute("widths"));

			double number;
			std::vector<double> binedge;
			std::vector<double> binwidth;
			while ( iss >> number ) binedge.push_back( number );
			while ( pss >> number ) binwidth.push_back( number );

			bin_edges.push_back(binedge);
			bin_widths.push_back(binwidth);
	
			nchan++;
			pChan = pChan->NextSiblingElement("channel");	
		}
	
		while(pMC)
			{
			num_multisim = strtod(pMC->Attribute("multisim"),&end);
			multisim_name = pMC->Attribute("name");


		        TiXmlElement *pBranch;
			pBranch = pMC->FirstChildElement("branch");
			while(pBranch){
				//std::cout<<pBranch->Attribute("name")<<" num_multisim"<<num_multisim<<std::endl;
				branch_names.push_back(pBranch->Attribute("name"));

				pBranch = pBranch->NextSiblingElement("branch");	
			}
			pMC=pMC->NextSiblingElement("MCevents");
		}
	



		num_channels = channel_names.size();
		num_modes = mode_names.size();
		num_detectors  = detector_names.size();


		if(false){

			std::cout<<" Covariance root path: "<<correlation_matrix_rootfile<<" and matrix name: "<<correlation_matrix_name<<std::endl;
			std::cout<<"Modes: ";
			for(auto b: mode_names){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Modes Bools: ";
			for(auto b: mode_bool){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Dets: ";
			for(auto b: detector_names){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Dets Bools: ";
			for(auto b: detector_bool){
				std::cout<<b<<" ";
			}
			std::cout<<std::endl;

			std::cout<<"Bin Edges: ";
			for(auto b: bin_edges[0]){
				std::cout<<b<<" ";
			}
			std::cout<<". Total#: "<<bin_edges[0].size()<<std::endl;

			std::cout<<"Bin Widths: ";
			for(auto b: bin_widths[0]){
				std::cout<<b<<" ";
			}
			std::cout<<". Total#: "<<bin_widths[0].size()<<std::endl;
	}





	std::string tempn;
	int indexcount = 0;
	for(int im = 0; im < num_modes; im++){

		for(int id =0; id < num_detectors; id++){

			for(int ic = 0; ic < num_channels; ic++){

				for(int sc = 0; sc < num_subchannels[ic]; sc++){

					tempn = mode_names[im] +"_" +detector_names[id]+"_"+channel_names[ic]+"_"+subchannel_names[ic][sc];
				
					// This is where you choose NOT to use some fields	
					if(mode_bool[im] && detector_bool[id] && channel_bool[ic] && subchannel_bool[ic][sc]){					
						fullnames.push_back(tempn);

						for(int k = indexcount; k < indexcount+num_bins[ic]; k++){
								used_bins.push_back(k);
							
						}
					}
					
					std::vector<int> tvec = {indexcount, indexcount+num_bins[ic]-1}; 

					mapIndex[tempn] = 	tvec;				
					indexcount = indexcount + num_bins[ic];
	

				}	
			}
		}
	}

	//For here on down everything is derivable, above is just until I actually get config working.
	num_modes = 0;
	for(bool i:mode_bool){	if(i) num_modes++;	}

	num_detectors = 0;
	for(bool i:detector_bool){	if(i) num_detectors++;	}
	
	num_channels = 0;
	for(bool i:channel_bool){	if(i) num_channels++;	}

	for(int i=0; i< num_channels; i++){
		num_subchannels[i] = 0;
		for(bool j: subchannel_bool[i]){ if(j) num_subchannels[i]++;}
	}



	num_bins_detector_block = 0;
	num_bins_detector_block_compressed = 0;
	for(int i =0; i< num_channels; i++){
		num_bins_detector_block += num_subchannels[i]*num_bins[i];
		num_bins_detector_block_compressed += num_bins[i];	
	}		
	num_bins_mode_block = num_bins_detector_block*num_detectors;
	num_bins_mode_block_compressed = num_bins_detector_block_compressed*num_detectors;
	num_bins_total = num_modes*num_bins_mode_block;
	num_bins_total_compressed = num_modes*num_bins_mode_block_compressed;


}//end constructor

