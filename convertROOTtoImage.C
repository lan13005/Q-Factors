// SIMPLE SCRIPT TO CONVERT ALL THE ROOT FILES IN THE HISTOGRAM FOLDER INTO A GIANT PDF

//gROOT->SetBatch(kTRUE);

void gatherAllImagesInSpecificFolder(string inDir){
	cout << "Looking in directory: " << inDir << endl;
	char* dir = gSystem->ExpandPathName(inDir.c_str());
	void* dirp = gSystem->OpenDirectory(dir);
        
        gStyle->SetLineScalePS(1);
	
	const char* entry;
	static const int maxFiles = 500;
	const char* subfolder[maxFiles];
	Int_t n = 0;
	TString str;
	
	while((entry = (char*)gSystem->GetDirEntry(dirp))) {
		str = entry;
		if(str.EndsWith(".root")){
		 	subfolder[n++] = gSystem->ConcatFileName(dir, entry);
			cout << "Adding file: " << str << endl;
		}
		if (n>maxFiles) {break;}
	}
	cout << "------------------\nThere are " << n << " files added" << endl; 
	
	TFile* histFile;
	TH1F* hist;
        int counter=0;
	for (Int_t i = 0; i < n; i++){
	  	Printf("Opening file -> %s", subfolder[i]);
	        string fileString = (string)subfolder[i];
		histFile = TFile::Open(subfolder[i]);
                if(histFile){ // make sure the file loaded, sometimes the root file can be corrupted if you exited the program improperly
		    TIter keyList(histFile->GetListOfKeys());
		    TKey *key;
   		    TCanvas *h;
   		    while ((key = (TKey*)keyList())) {
   		       	TClass *cl = gROOT->GetClass(key->GetClassName());
                        string keyName = (string)key->GetName();
                        if (keyName.find("BS") == std::string::npos) {
   		       	    if (cl->InheritsFrom("TCanvas")){
                                h = (TCanvas*)key->ReadObj();
                                //cout << key->GetName() << endl;
                                //h->SetTitle(key->GetName());
		    	    	Printf("Drawing canvas: %s", h->GetName());
                                // --- FOR SOME REASON I NEED TO DRAW CLONE OR DO SOME COUT STATEMENTS TO FULLY COMPLETE THE LOOP. AN ERROR STILL OCCURS SAYING THE LIST IS TRYING TO ACCESS A DELETED VARIABLE
                                // --- BUT IT ALLOWS ME TO COMPLETE THE PROGRAM.
		    	    	h->Draw();//Clone();
		    	    	//h->SaveAs((inDir+"/"+key->GetName()+".png").c_str());

                                if ( counter==0 ){
		    	    	    h->Print((inDir+"/massHists.pdf(").c_str(),"pdf");
                                }
                                // This only works if there is only one key inside each root file. The counter wont match properly otherwise
                                else if ( counter == n-1 ){
		    	    	    h->Print((inDir+"/massHists.pdf)").c_str(),"pdf");
                                }
                                else{
		    	    	    h->Print((inDir+"/massHists.pdf").c_str(),"pdf");
                                }
                                ++counter;
		    	    }
                        }
		    }
		    histFile->Close();	
                }
	}

	cout << "\n\n" << endl;
	//cout << "rm " << inDir << "/*.root" << endl;
	//gSystem->Exec(("rm "+inDir+"/*.root").c_str());
}

void convertROOTtoImage(){
	gatherAllImagesInSpecificFolder("histograms/flatEtapi");
}
