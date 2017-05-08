void makeMC(){

	TFile f1("test.root","RECREATE");


    	 TTree t("MCtree","");
	 int numsim =5;
	 double Enu; 
	 double Etype ;
	 
	 std::vector<double>  weights;

	 t.Branch("weights", &weights);
         t.Branch("Enu",&Enu,"Enu/D");
	 t.Branch("Etype",&Etype,"Etype/D");


	TRandom3 *rangen = new TRandom3();

	 for(int i=0; i<20000; i++){
		 weights.clear();
		 weights.push_back(1);
		for(int k=1; k<numsim; k++){
			weights.push_back(rangen->Uniform(0,2));
		}
		double tt=rangen->Gaus(1,0.5);
		while(tt<0){
			tt=rangen->Gaus(1,0.5);
		}
		Enu = tt;
		Etype =1;	
		t.Fill();
	 }

	 t.Write();
	f1.Close();

	return;
}
