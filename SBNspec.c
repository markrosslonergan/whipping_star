#include "SBNspec.h"
using namespace SBNFIT;


SBNspec::SBNspec(const char * name, std::string whichxml) : SBNconfig(whichxml) {


	char namei[200];
	sprintf(namei,"%s.root",name);	
	TFile f(namei);



	for(auto fn: fullnames){ 	//Loop over all filenames that should be there, and load up the histograms.
		//std::cout<<"Attempting to load: "<<fn.c_str()<<" from: "<<namei<<std::endl;
		hist.push_back(*((TH1D*)f.Get(fn.c_str()))); 
	}


	f.Close();
}//end constructor



int SBNspec::Add(SBNspec *in){

	if(xmlname != in->xmlname){ std::cout<<"ERROR: SBNspec::Add, trying to add differently configured SBNspecs!"<<std::endl; exit(EXIT_FAILURE);}

	for(int i=0; i< hist.size(); i++){
		hist[i].Add( &(in->hist[i]));
	}	


	return 0;
}


int SBNspec::poissonScale(){
	TRandom3 *rangen = new TRandom3(0);
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
		}	
	}
return 0;
}


int SBNspec::randomScale(){
	TRandom3 *rangen    = new TRandom3(0);

	for(auto& h: hist){
			h.Scale(rangen->Uniform(0,2));

	}
return 0;
}


int SBNspec::Scale(std::string name, TF1 * func){
	for(auto& h: hist){
		std::string test = h.GetName();
			if(test.find(name)!=std::string::npos ){
				for(int b=0; b<=h.GetNbinsX(); b++){
					//std::cout<<h.GetBinContent(b)<<" "<<h.GetBinCenter(b)<<" "<<func->Eval(h.GetBinCenter(b) )<<std::endl; 
					h.SetBinContent(b, h.GetBinContent(b)*func->Eval(h.GetBinCenter(b) ) );
				}
			}

	}
return 0;
}


int SBNspec::ScaleAll(double sc){
	for(auto& h: hist){
		h.Scale(sc,"nosw2");
	}
return 0;
}

int SBNspec::Scale(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();
		
			if(test.find(name)!=std::string::npos ){
				//std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
				h.Scale(val,"nosw2" );
			}

	}
return 0;
}

int SBNspec::NormAll(double n){

	for(auto& h: hist) {
		h.Scale(n/h.GetSumOfWeights());
	}
	return 0;
}

int SBNspec::Norm(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();
		
			if(test.find(name)!=std::string::npos ){
				//std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
				h.Scale(val/h.GetSumOfWeights());
			}

	}
return 0;
}

int SBNspec::calcFullVector(){
	fullVec.clear();

	for(auto& h: hist){
		//std::cout<<"Hist size: "<<h.GetSize()-2<<std::endl;
		for(int i = 1; i <= h.GetSize()-2; i++){
			//std::cout<<h.GetBinContent(i)<<" ";
			fullVec.push_back(h.GetBinContent(i));
		}	
	}

return 0;
}

int SBNspec::compressVector(){
	compVec.clear();
	//This needs to be confirmed and checked. Looks good, mark 24th april
	calcFullVector();
	for(int im = 0; im < Nmode; im++){
		for(int id =0; id < Ndet; id++){
			int edge = id*Tdet + Tmode*im; // This is the starting index for this detector


			for(int ic = 0; ic < Nchan; ic++){
				int corner=edge;
				for(int j=0; j< Bins[ic]; j++){
					double tempval=0;

					for(int sc = 0; sc < Chan[ic]; sc++){
						tempval += fullVec[j+sc*Bins[ic]+corner];
						edge +=1;	//when your done with a channel, add on every bin you just summed
					}
					compVec.push_back(tempval);
				}
				
				
				
									


			}
		}
	}
return 0;
}

int SBNspec::printFullVec(){
	for(double d: fullVec){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
return 0;
}

int SBNspec::printCompVec(){ 
	for(double d: compVec){
		std::cout<<d<<" ";
	}
	std::cout<<std::endl;
return 0;
}


int SBNspec::writeOut(std::string filename){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

	std::vector<int> mycol = {416-6, 800+3, 616+1, 632-7, 600-7, 432+1, 900}; 				

	TFile *f = new TFile(filename.c_str(),"RECREATE" ); 
	f->cd();

	std::vector<TH1D> temp = hist;

	for(int i =0; i < Nchan; i++)
	{
	TCanvas* Cstack= new TCanvas(cname[i].c_str(),cname[i].c_str());
	Cstack->cd();
	THStack * hs 	   = new THStack(cname[i].c_str(),  cname[i].c_str());
	TLegend * legStack = new TLegend(0.6,0.35,0.875,0.875);
		int n=0;
		for(auto &h : temp){
			std::string test = h.GetName();
			if(test.find(cname[i])!=std::string::npos ){
				h.Scale(1,"width,nosw2");
				h.GetYaxis()->SetTitle("Events/GeV");
				h.SetMarkerStyle(20);
				h.SetMarkerColor(mycol[n]);
				h.SetFillColor(mycol[n]);
				h.SetLineColor(kBlack);
				h.SetTitle(h.GetName());
				h.Write();

				legStack->AddEntry(&h,h.GetName() , "f");
	
				hs->Add(&h);
				n++;

			}
		}
	/****Not sure why but this next line seg faults...******
	*	hs->GetYaxis()->SetTitle("Events/GeV");
	******************************************************/

	hs->Draw();
	Cstack->Update();
	legStack->Draw();	
	Cstack->Write();
	}

	f->Close();

	return 0;
}


