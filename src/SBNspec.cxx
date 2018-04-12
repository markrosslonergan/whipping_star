#include "SBNspec.h"
using namespace sbn;

SBNspec::SBNspec(std::string whichxml): SBNspec(whichxml,-1){
}


SBNspec::SBNspec(std::string whichxml, int which_universe) : SBNconfig(whichxml){

	//Initialise all the things
	//for every multisim, create a vector of histograms, one for every subchannel we want 
	int ctr=0;
	for(auto fn: fullnames){
		for(int c=0; c<channel_names.size(); c++){
			if(fn.find(channel_names[c])!=std::string::npos ){
				double * tbins =&bin_edges[c][0];
				std::string thisname;
				if(which_universe<0){
					thisname = fn;
				}else{
					thisname = fn+"_MS"+std::to_string(which_universe);
				}
				TH1D thischan(thisname.c_str(),"",num_bins[c], tbins );
				hist.push_back(thischan);
				//auto it = hist.begin()+ctr;
				//map_hist[fn] = &(*it);
				map_hist[fn] = ctr;

				ctr++;
			}

		}
	}

	has_been_scaled = false;


}

SBNspec::SBNspec(const char * name, std::string whichxml, bool isverbose) : SBNconfig(whichxml, isverbose) {
	//Contruct from a prexisting histograms!

	char namei[200];
	sprintf(namei,"%s.root",name);	
	TFile f(namei);

	//Loop over all filenames that should be there, and load up the histograms.
	for(auto fn: fullnames){
		//std::cout<<"Attempting to load: "<<fn.c_str()<<" from: "<<namei<<std::endl;
		hist.push_back(*((TH1D*)f.Get(fn.c_str()))); 
	}
	has_been_scaled=false;


	f.Close();


}//end constructor



SBNspec::SBNspec(const char * name, std::string whichxml) : SBNspec(name, whichxml, true)  {

}//end constructor


int SBNspec::Clear(){
	for(auto &h: hist){
		h.Reset();
	}	

	return 0;
}



int SBNspec::Add(SBNspec *in){
	//Addes all hists together
	if(xmlname != in->xmlname){ std::cout<<"ERROR: SBNspec::Add, trying to add differently configured SBNspecs!"<<std::endl; exit(EXIT_FAILURE);}

	for(int i=0; i< hist.size(); i++){
		hist[i].Add( &(in->hist[i]));
	}	


	return 0;
}

int SBNspec::setAsGaussian(double mean, double sigma, int ngen){
	TRandom3 *seedgetter = new TRandom3(0);
	int seed = seedgetter->Integer(1000000);

	for(auto &h: hist){
		TRandom3 *rangen = new TRandom3(seed);
		h.Reset();
		for(int i=0; i<ngen; i++){
			double eve = rangen->Gaus(mean,sigma); 
			h.Fill( eve ); 
		}	
	}

	return 0;

}

int SBNspec::setAsFlat(double val){
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, val );
		}	
	}
}



//All scaling functions are quite self explanatory
int SBNspec::poissonScale(){
	TRandom3 *rangen = new TRandom3(0);
	for(auto &h: hist){
		for(int i=0; i<h.GetSize(); i++){
			h.SetBinContent(i, rangen->Poisson( h.GetBinContent(i)    ));
		}	
	}
	return 0;
}

int SBNspec::poissonScale(TRandom3* rangen){
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
		h.Scale(sc);
	}
	this->compressVector();

	return 0;
}

int SBNspec::Scale(std::string name, double val){
	for(auto& h: hist){
		std::string test = h.GetName();

		if(test.find(name)!=std::string::npos ){
			//	std::cout<<name<<". found in: "<<test<<" at "<<test.find(name)<<std::endl;
			h.Scale(val);
		}

	}

	has_been_scaled = true;
	scale_hist_name =name;
	scale_hist_val = val;

	this->compressVector();
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
	//std::cout<<"num_modes: "<<num_modes<<" num_detectors: "<<num_detectors<<" num_channels: "<<num_channels<<std::endl;

	for(int im = 0; im < num_modes; im++){
		for(int id =0; id < num_detectors; id++){
			int edge = id*num_bins_detector_block + num_bins_mode_block*im; // This is the starting index for this detector

			for(int ic = 0; ic < num_channels; ic++){
				int corner=edge;

				for(int j=0; j< num_bins.at(ic); j++){

					double tempval=0;


					for(int sc = 0; sc < num_subchannels.at(ic); sc++){

						//std::cout<<im<<"/"<<num_modes<<" "<<id<<"/"<<num_detectors<<" "<<ic<<"/"<<num_channels<<" "<<j<<"/"<<num_bins[ic]<<" "<<sc<<"/"<<num_subchannels[ic]<<std::endl;
						tempval += fullVec.at(j+sc*num_bins.at(ic)+corner);
						edge +=1;	//when your done with a channel, add on every bin you just summed
					}
					compVec.push_back(tempval);
				}






			}
		}
	}
	return 0;
}

double SBNspec::getTotalEvents(){
	double ans =0;
	this->calcFullVector();

	for(double d: fullVec){
		ans+=d;
	}

	return ans;


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

int SBNspec::writeSpec(std::string filename){

	std::string fn1= filename+".SBNspec.root";
	TFile *f2 = new TFile(fn1.c_str(),"RECREATE" ); 

	for(auto& h: hist){
		h.Write();
	}
	f2->Close();	

}


int SBNspec::writeOut(std::string filename){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
	std::vector<int> mycol = {kRed-7, kRed+1, kYellow-7, kOrange-3, kBlue+3, kBlue,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};

	std::string fn1= "SBN_"+filename;
	TFile *f2 = new TFile(fn1.c_str(),"RECREATE" ); 

	for(auto& h: hist){
		h.Write();
	}
	f2->Close();	


	TFile *f = new TFile(filename.c_str(),"RECREATE" ); 
	f->cd();

	std::vector<TH1D> temp = hist;



	for(auto m: mode_names){
		for(auto d: detector_names){
			for(auto c: channel_names){

				std::string canvas_name = m+"_"+d+"_"+c;

				bool this_run = false;

				TCanvas* Cstack= new TCanvas(canvas_name.c_str(),canvas_name.c_str());
				Cstack->cd();
				THStack * hs 	   = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.6,0.35,0.875,0.875);
				int n=0;
				for(auto &h : temp){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos ){
						double total_events = h.GetSumOfWeights();
						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
						h.GetYaxis()->SetTitle("Events/GeV");
						h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol[n]);
						h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						h.Write();

						std::ostringstream out;
						out << std::setprecision(6) << total_events;
						std::string hmm = " , Total : ";
						std::string tmp = h.GetName() +hmm+ out.str();
						legStack.AddEntry(&h, tmp.c_str() , "f");

						hs->Add(&h);
						n++;

						this_run=true;

					}
				}
				/****Not sure why but this next line seg faults...******
				 *	hs->GetYaxis()->SetTitle("Events/GeV");
				 ******************************************************/
				if(this_run){
					hs->Draw();
					Cstack->Update();
					legStack.Draw();	
					Cstack->Write("hist");
				}

			}
		}
	}

	f->Close();

	return 0;
}



int SBNspec::getLocalBinNumber(double invar, int which_hist)
{
	int localbin = hist.at(which_hist).GetXaxis()->FindBin(invar);
	double bin = localbin-1;

	if(localbin==0 || localbin > hist.at(which_hist).GetNbinsX() ){bin = -99;}
	return bin;
}


int SBNspec::getGlobalBinNumber(double invar, int which_hist)
{
	int localbin = hist.at(which_hist).GetXaxis()->FindBin(invar);
	double bin = localbin-1;

	for(int i=0; i<which_hist; i++){
		bin += hist.at(i).GetNbinsX();	
	}

	if(localbin==0 || localbin > hist.at(which_hist).GetNbinsX() ){bin = -99;}
	return bin;
}


int SBNspec::compareSBNspecs(SBNspec * compsec, std::string filename){
	//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
	//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
	//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
	std::vector<int> mycol = {kRed-7, kRed+1, kYellow-7, kOrange-3, kBlue+3, kBlue,  kGreen+1,kBlue-7, kPink, kViolet, kCyan,kMagenta,kAzure};


	TFile *f = new TFile(filename.c_str(),"RECREATE" ); 
	f->cd();


	std::vector<TH1D> temp = hist;
	std::vector<TH1D> temp_comp = compsec->hist;

	for(int k=0; k< fullnames.size(); k++){
		TCanvas *ctmp = new TCanvas((std::to_string(k)+"_"+fullnames.at(k)).c_str(), (std::to_string(k)+"_"+fullnames.at(k)).c_str());
		ctmp->cd();
		TH1D * h1 = (TH1D*) temp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_1").c_str());
		TH1D * h2 = (TH1D*) temp_comp.at(map_hist[fullnames.at(k)]).Clone((std::to_string(k)+fullnames.at(k)+"_2").c_str());
	
		h1->Scale(1,"width");	
		h2->Scale(1,"width");	

		h1->SetLineColor(kRed);
		h1->SetLineWidth(2);
		h1->Draw("hist");
		
		h2->SetLineColor(kBlue);
		h2->SetLineStyle(5);
		h2->SetLineWidth(2);
		h2->Draw("hist same");		

		h1->SetMaximum(std::max(h1->GetMaximum(), h2->GetMaximum())*1.10);
		

		ctmp->Write();
	}
		



	for(auto m: mode_names){
		for(auto d: detector_names){
			for(auto c: channel_names){

				std::string canvas_name = m+"_"+d+"_"+c;

				bool this_run = false;
				bool this_run_comp = false;

				TCanvas* Cstack= new TCanvas(canvas_name.c_str(),canvas_name.c_str());
				Cstack->cd();
				THStack * hs = new THStack(canvas_name.c_str(),  canvas_name.c_str());
				TLegend legStack(0.6,0.35,0.875,0.875);
				int n=0;
				int nc=0;
				TH1D * hcomp;				
				TH1D *hsum;
	
		

				for(auto &h : temp_comp){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos ){
						double total_events = h.GetSumOfWeights();


						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
						//h.GetYaxis()->SetTitle("Events/GeV");
						//h.SetMarkerStyle(20);
						//h.SetMarkerColor(mycol[n]);
						//h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						//h.SetTitle(h.GetName());
						
				
						if(!this_run_comp){
							hcomp = (TH1D*)h.Clone(("comp_"+canvas_name).c_str());
							hcomp->Reset();
						}

						std::ostringstream out;
						out << std::setprecision(6) << total_events;
						std::string hmm = " , Total : ";
						std::string tmp = h.GetName() +hmm+ out.str();


						hcomp->Add(&h);
						nc++;

						this_run_comp=true;

					}
				}


				for(auto &h : temp){
					std::string test = h.GetName();
					if(test.find(canvas_name)!=std::string::npos ){
						double total_events = h.GetSumOfWeights();
						h.Sumw2(false);
						h.Scale(1,"width,nosw2");
						h.GetYaxis()->SetTitle("Events/GeV");
						h.SetMarkerStyle(20);
						h.SetMarkerColor(mycol[n]);
						h.SetFillColor(mycol[n]);
						h.SetLineColor(kBlack);
						h.SetTitle(h.GetName());
						//h.Write();

						if(!this_run){
							hsum = (TH1D*)h.Clone(("sum_"+canvas_name).c_str());
							hsum->Reset();
						}

						std::ostringstream out;
						out << std::setprecision(6) << total_events;
						std::string hmm = " , Total : ";
						std::string tmp = h.GetName() +hmm+ out.str();
						legStack.AddEntry(&h, tmp.c_str() , "f");
						hsum->Add(&h);
						hs->Add(&h);
						n++;

						this_run=true;

					}
				}
		

				/****Not sure why but this next line seg faults...******
				 *	hs->GetYaxis()->SetTitle("Events/GeV");
				 ******************************************************/
				if(this_run && this_run_comp){
					double plot_pot=5e19;

					double title_size_ratio=0.1;
					double label_size_ratio=0.1;
					double title_offset_ratioY = 0.3 ;
					double title_offset_ratioX = 1.1;

					double title_size_upper=0.15;
					double label_size_upper=0.05;
					double title_offset_upper = 1.45;




					Cstack->cd();	
					TPad *pad0top = new TPad(("pad0top_"+canvas_name).c_str(), ("pad0top_"+canvas_name).c_str(), 0, 0.35, 1, 1.0);
					pad0top->SetBottomMargin(0); // Upper and lower plot are joined
					pad0top->Draw();             // Draw the upper pad: pad2top
					pad0top->cd();               // pad2top becomes the current pad
					hs->Draw();
					hcomp->SetLineColor(kBlack);
					hcomp->SetMarkerColor(kBlack);
					hcomp->SetMarkerStyle(20);

					hcomp->Draw("E1 same");
					Cstack->Update();
					legStack.Draw();

					Cstack->cd();

					TPad *pad0bot = new TPad(("padbot_"+canvas_name).c_str(),("padbot_"+canvas_name).c_str(), 0, 0.05, 1, 0.35);
					pad0bot->SetTopMargin(0);
					pad0bot->SetBottomMargin(0.351);
					pad0bot->SetGridx(); // vertical grid
					pad0bot->Draw();
					pad0bot->cd();       // pad0bot becomes the current pad

					TH1* ratpre = (TH1*)hsum->Clone(("ratio_"+canvas_name).c_str());
					ratpre->Divide(hcomp);		
					ratpre->Draw("hist");	
					ratpre->SetFillColor(kWhite);
					ratpre->SetFillStyle(0);
					ratpre->SetLineWidth(2);

					TLine *line = new TLine(ratpre->GetXaxis()->GetXmin(),1.0,ratpre->GetXaxis()->GetXmax(),1.0 );
					line->Draw("same");
					ratpre->SetLineColor(kBlack);
					ratpre->SetTitle("");
					ratpre->GetYaxis()->SetTitle("Ratio");
					ratpre->GetXaxis()->SetTitleOffset(title_offset_ratioX);
					ratpre->GetYaxis()->SetTitleOffset(title_offset_ratioY);
					ratpre->SetMinimum(0);	
					ratpre->SetMaximum(2);
					ratpre->GetYaxis()->SetTitleSize(title_size_ratio);
					ratpre->GetXaxis()->SetTitleSize(title_size_ratio);
					ratpre->GetYaxis()->SetLabelSize(label_size_ratio);
					ratpre->GetXaxis()->SetLabelSize(label_size_ratio);
					ratpre->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");

					Cstack->Write("hist");
				}

			}
		}
	}

	f->Close();

	return 0;
}


