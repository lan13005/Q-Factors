#include "main.h"

using namespace RooFit;
int extra=0;

void QFactorAnalysis::initialize(string rootFileLoc, string rootTreeName){
	cout << "Loading root file and tree" << endl;
        if (saveMemUsage)
            outputMemUsage(pinfo,"Before QFactorAnalysis initialization: ");
        dataFile=new TFile((rootFileLoc).c_str());
	dataFile->GetObject((rootTreeName).c_str(),dataTree);
	// Get the total number of entries and potentially overwrite it if we want to have a shorter run
	total_nentries = (Long64_t)dataTree->GetEntries();
	if (!override_nentries){
		nentries=dataTree->GetEntries();
	}
	cout << "Chosen Total Entries: " << nentries << endl;
        
        cout << "Reserving space in vectors..." << endl;
	// Use these vectors to import all the data to RAM instead of reading from root file
        // First we should reserve the space so there is no resizing of the vectors
        std::vector<float> emptyVec;
        for (auto iVar=0; iVar<phaseSpaceDim; ++iVar){
            phaseSpaceVars.push_back(emptyVec);
            phaseSpaceVars[iVar].reserve(nentries);
        }
	for (int iVar=0; iVar<discrimVarDim; ++iVar){
            // discrimVars will always have dimension 2 (at most 2D fits can be done). We can reserve and fill 1 or 2 dimensions
            discrimVars.push_back(emptyVec);
            discrimVars[iVar].reserve(nentries);
        }
	accWeights.reserve(nentries);
	// will hold all the ids of the unique combos
	phasePoint2PotentialNeighbor.reserve(nentries);

        if (redistributeBkgSigFits) { sigFracs={0,0.5,1}; }
        else { sigFracs={0.5}; }

        if (saveMemUsage)
            outputMemUsage(pinfo,"After QFactorAnalysis initialization: ");
}


void QFactorAnalysis::loadData(){
        if (saveMemUsage)
            outputMemUsage(pinfo,"\n\nBefore loading data: ");
	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         LOAD IN THE DATA
	// -----------------------------------------------------
	// -----------------------------------------------------
	// Create variables to hold the data as we read in the data from the tree
        // We MUST match the data type of the input branch. There are three cases I will consider: {float, double, Long64} 
        //   We will keep track of what data type a specific branch used and fill accordingly. Only these variables we
        //   use to SetBranchAddress needs to match. We can fill a double or a float into a vector<float>. Using floats
        //   saves memory so its worth the effort
        double phaseSpaceVar[phaseSpaceDim];
        double discrimVar[discrimVarDim];
        float phaseSpaceVar_f[phaseSpaceDim];
        float discrimVar_f[discrimVarDim];
        Long64_t phaseSpaceVar_l[phaseSpaceDim];
        Long64_t discrimVar_l[discrimVarDim];
        string phaseSpaceDataTypes[phaseSpaceDim];
        string discrimVarDataTypes[discrimVarDim];  
	double accWeight;
	float accWeight_f;
        Long64_t accWeight_l;
        string weightType;

	// Set branch addresses so we can read in the data
        string typeName;

        dataTree->SetBranchStatus("*",0);
	parsePhaseSpace.parseString(s_phaseVar);
        if ( parsePhaseSpace.varStringSet.size() != phaseSpaceDim ) { cout << "Uh-oh something went wrong. varString size not right size" << endl; }
	for (int iVar=0; iVar<phaseSpaceDim; ++iVar){
            dataTree->SetBranchStatus(parsePhaseSpace.varStringSet[iVar].c_str(),1);
            typeName=setBranchAddress(dataTree, parsePhaseSpace.varStringSet[iVar], &phaseSpaceVar_l[iVar], &phaseSpaceVar_f[iVar], &phaseSpaceVar[iVar]);
            phaseSpaceDataTypes[iVar]=typeName;
        }
	parseDiscrimVars.parseString(s_discrimVar);
        if ( parseDiscrimVars.varStringSet.size() != discrimVarDim ) { cout << "Uh-oh something went wrong. discrimVar string size not right size" << endl; }
	for (int iVar=0; iVar<discrimVarDim; ++iVar){
            dataTree->SetBranchStatus(parseDiscrimVars.varStringSet[iVar].c_str(),1);
            typeName=setBranchAddress(dataTree, parseDiscrimVars.varStringSet[iVar], &discrimVar_l[iVar], &discrimVar_f[iVar], &discrimVar[iVar]);
            discrimVarDataTypes[iVar]=typeName;
        }
        
        //////////////////////////////////////////////////////////
        // SET UP WEIGHTS
        //////////////////////////////////////////////////////////
        // Setting up tracking of accidental weights
        if (!s_accWeight.empty()){ // if string is not empty we will set the branch address
            dataTree->SetBranchStatus(s_accWeight.c_str(),1);
            weightType=setBranchAddress(dataTree, s_accWeight, &accWeight_l, &accWeight_f, &accWeight);
            cout << "Using weights in branch: "+s_accWeight << endl;
        }
        else{
            accWeight_f=1;
            weightType="Float_t"; // we wont be using weights, just keep as a value of 1 with datatype float
            cout << "No additional weights used" << endl;
        }

        /////////////////////////////////////////////////////////////
        // LOAD THE DATA
        /////////////////////////////////////////////////////////////
	// We will use a ientry to keep track of which entries we will get from the tree. We will simply use ientry when filling the arrays.  
        if (saveMemUsage)
            outputMemUsage(pinfo,"Begin loading data: ");
	for (Long64_t ientry=0; ientry<nentries; ientry++)
	{
		dataTree->GetEntry(ientry);
	        for (int iVar=0; iVar<phaseSpaceDim; ++iVar){
                    if (phaseSpaceDataTypes[iVar]=="Float_t")
                        phaseSpaceVars[iVar].push_back(phaseSpaceVar_f[iVar]);
                    if (phaseSpaceDataTypes[iVar]=="Double_t")
                        phaseSpaceVars[iVar].push_back(phaseSpaceVar[iVar]);
                    if (phaseSpaceDataTypes[iVar]=="Long64_t")
                        phaseSpaceVars[iVar].push_back(phaseSpaceVar_l[iVar]);
                    //if (saveMemUsage)
                    //    outputMemUsage(pinfo,"loaded phaseSpace iVar "+to_string(iVar)+": ");
                }
	        for (int iVar=0; iVar<discrimVarDim; ++iVar){
                    if (discrimVarDataTypes[iVar]=="Float_t")
                        discrimVars[iVar].push_back(discrimVar_f[iVar]);
                    if (discrimVarDataTypes[iVar]=="Double_t")
                        discrimVars[iVar].push_back(discrimVar[iVar]);
                    if (discrimVarDataTypes[iVar]=="Long64_t")
                        discrimVars[iVar].push_back(discrimVar_l[iVar]);
                    //if (saveMemUsage)
                    //    outputMemUsage(pinfo,"loaded discrimVar iVar "+to_string(iVar)+": ");
                }
                if (weightType=="Float_t")
	            accWeights.push_back(accWeight_f);
                if (weightType=="Double_t")
	            accWeights.push_back(accWeight);
                if (weightType=="Long64_t")
	            accWeights.push_back(accWeight_l);
                //if (saveMemUsage)
                //    outputMemUsage(pinfo,"loaded accWeight: ");
	}
        if (saveMemUsage)
            outputMemUsage(pinfo,"After loading data: ");

	if ( verbose_outputDistCalc ) {
	    cout << "Before standarization the first nentries of phaseSpaceVar[0]" << endl;
	    for ( int ientry=0 ; ientry < nentries; ientry++){
	        cout << phaseSpaceVars[0][ientry] << endl;
	    }
	}

        if (saveMemUsage)
            outputMemUsage(pinfo,"After standardization: ");
        

        dataFile->Close(); // dataTree is owned by dataFile. So it will automatically free memory
        gSystem->GetProcInfo(&pinfo);
        if (saveMemUsage)
            outputMemUsage(pinfo,"After closing dataFile: ");

	// -----------------------------------------------------
	// -----------------------------------------------------
	//                         STANDARDIZE THE DATA
        //                         Range vs StdDev standardization
        //                         Not much difference so far but the option is left here
	// -----------------------------------------------------
	// -----------------------------------------------------
	//
        cout << "standardizationType flag was set to : " << standardizationType << endl;
        if(standardizationType=="range")
            cout << "Using range standardization" << endl;
        else if (standardizationType=="std")
            cout << "Using std standardization" << endl;
        else
            cout << "Not standardizing phase space" << endl;
	standardizeArray standarizationClass;
	for (int iVar=0; iVar<phaseSpaceDim; ++iVar){
            phaseSpaceVars[iVar].push_back(phaseSpaceVar[iVar]);

            if(standardizationType=="range"){
                standarizationClass.rangeStandardization(phaseSpaceVars[iVar],nentries);
            }
            else if (standardizationType=="std"){
                standarizationClass.stdevStandardization(phaseSpaceVars[iVar],nentries);
            }
        }

	if ( verbose_outputDistCalc ) {
	    cout << "After standarization the first nentries of phaseSpaceVar[0]" << endl;
	    for ( int ientry=0 ; ientry < nentries; ientry++){
	        cout << phaseSpaceVars[0][ientry] << endl;
	    }
	}
        gSystem->GetProcInfo(&pinfo);
        if (saveMemUsage)
            outputMemUsage(pinfo,"After standardization: ");

	
	// phasePoint1 will consider all events from lowest to largest since these will be our attached q values. 
        // phasePoint2 can use a random subset (hopefully representative of the entire dataset), or accept/reject an entry given a certain criteria, i.e. if some unique indentifier has been seen before
        //      Currently just considers all other entries as potential neighbors
        if(nRndRepSubset>0 && nRndRepSubset<nentries){
            std::set<Int_t> rndRepSubset; //probably should make the subset have unique elements. So sampling without replacement
            srand (time(NULL)); // set a random seed
            while(rndRepSubset.size() < nRndRepSubset){
                rndRepSubset.insert(rand()%nentries);
            }
            phasePoint2PotentialNeighbor.assign(rndRepSubset.begin(),rndRepSubset.end());
            rndRepSubset.clear();
        }
        else {
	    for (Int_t ientry=0; ientry<nentries; ientry++){ 
	        phasePoint2PotentialNeighbor.push_back(ientry);
	    }
        }
        cout << "\n\n--------------------------------------" << endl;
	cout << phasePoint2PotentialNeighbor.size() << "/" << nentries << " are used as potential neighbors" << endl;
        cout << "Some indicies that will be used as potential neighbors: " << endl;
        for (auto iPhasePoint2=0; iPhasePoint2 < 20; ++iPhasePoint2){
            cout << phasePoint2PotentialNeighbor[iPhasePoint2] << " ";
        }
        cout << endl;
        gSystem->GetProcInfo(&pinfo);
        if (saveMemUsage)
            outputMemUsage(pinfo,"After determining representative subset: ");

}

// This method was originally designed for multithreading. Turns out RooFit is not thread safe we we had to resort back to spawning multiple root processes. The orignially uses a lambda functon which contains the code to extract the q-values in batches. Some extra code to spawn the threads and waits for all of them to execute is there also
void QFactorAnalysis::runQFactorThreaded(int iProcess){
	// make sure the global variables are read in correctly
	cout << "kDim: " << kDim << endl;
	cout << "numberEventsToSavePerProcess: " << numberEventsToSavePerProcess << endl;
	cout << "nProcess: " << nProcess << endl;
	cout << "nentries: " << nentries << endl;
	cout << "override_nentries: " << override_nentries << endl;

	// [=] refers to a capture list which is used by this lambda expression. The lambda gets a copy of all the local variables that it uses when it is created. If we
	// just use [] we will get an error since the lambda will have no idea what these variables are	
	//ROOT::EnableThreadSafety();
	//auto f = [=](int iProcess){
	// ----------------------------------------
	// Open up a root file to save the q-factors and other diagnostics to it
	// ---------------------------------------
        // Initializing some variables we can track during the q-value extraction
        ULong64_t flatEntryNumber;
        float NLL;
        float bestNLL;
	float worstNLL;
        float qvalue;
        float best_qvalue;
        int fitStatus;
        int best_fitStatus;
        float worst_qvalue;
        float eff_nentries;
        float qvalueBS_std=0;
        float best_nsig;
        float best_nbkg;
        float best_ntot;
        float effNentriesMinusTotal;
        int neighbors[kDim];

        // Saving the results along with some diagnostics
        TFile *resultsFile = new TFile((cwd+"/logs"+runTag+"/"+fileTag+"/results"+to_string(iProcess)+".root").c_str(),"RECREATE");
        TTree* resultsTree = new TTree(rootTreeName.c_str(),"results");
        resultsTree->Branch("flatEntryNumber",&flatEntryNumber,"flatEntryNumber/l");
        resultsTree->Branch("qvalue",&best_qvalue,"qvalue/F");
        resultsTree->Branch("fitStatus",&best_fitStatus,"status/I");
        resultsTree->Branch("worst_qvalue",&worst_qvalue,"worst_qvalue/F");
        resultsTree->Branch("qvalueBS_std",&qvalueBS_std,"qvalueBS_std/F");
        resultsTree->Branch("bestNLL",&bestNLL,"bestNLL/F");
        resultsTree->Branch("worstNLL",&worstNLL,"worstNLL/F");
        resultsTree->Branch("best_nsig",&best_nsig,"best_nsig/F");
        resultsTree->Branch("best_nbkg",&best_nbkg,"best_nbkg/F");
        resultsTree->Branch("best_ntot",&best_ntot,"best_ntot/F");
        resultsTree->Branch("eff_nentries",&eff_nentries,"eff_nentries/F");
        resultsTree->Branch("effNentriesMinusTotal",&effNentriesMinusTotal,"effNentriesMinusTotal/F");
        resultsTree->Branch("kDim",&kDim,"kDim/I"); // 32 bit integer. Could have used unsigned but not really worth the change...
        if (saveBranchOfNeighbors){
            resultsTree->Branch("neighbors",neighbors,"neighbors[kDim]/I"); // I = 32 bit integer
        }
        cout << "Set up branch addresses" << endl;

	// Define some needed variables like canvases, histograms, and legends
	cout << "Creating canvas " << iProcess << endl;
    	TCanvas *allCanvases = new TCanvas(("anyHists"+to_string(iProcess)).c_str(),"",1440,900);
        TH1F* dHist_qvaluesBS = new TH1F(("qvaluesBS"+to_string(iProcess)).c_str(),"Bootstrapped Q-Factors",100,0,1);

        // Variables to keep track of phase space neighbors
	float phasePoint1[phaseSpaceDim];
	float phasePoint2[phaseSpaceDim];
	float distance;
        distSort_kNN distKNN(kDim);
        pair<float,int> newPair;
        
        // Track current discriminating variable values
        float currentValues[discrimVarDim];
        float neighborValues[discrimVarDim];

	// opening a file to write my log data to
    	ofstream logFile;
    	logFile.open((cwd+"/logs"+runTag+"/"+fileTag+"/processLog"+to_string(iProcess)+".txt").c_str());
	
	// Determine what events each thread should run
	int batchEntries = (int)nentries/nProcess; // batchEntries the size of the batch
	int lowest_nentry = iProcess*batchEntries;
	int largest_nentry;
        if (iProcess!=(nProcess-1)) {
            largest_nentry  = (iProcess+1)*batchEntries;
        }
        else {
            largest_nentry = nentries; 
        }
	cout << "nentries we will use for this process: " << lowest_nentry << ", " << largest_nentry << endl;

	// randomly select some events (but deterministically since we migth want to double check) to write histograms for 
        parseEventsToSave.parseString(alwaysSaveTheseEvents);
        TFile *qHistsFile;
	set<int> selectRandomIdxToSave;
        bool saveAllHistograms=false;
        if (numberEventsToSavePerProcess>=0){
	    int randomEvent;
	    srand(iProcess+seedShift);
	    for (int i=0; i<numberEventsToSavePerProcess; i++){
	    	// basically randomly sample a uniform number between (lowest_nentry, largest_nentry)
	    	randomEvent = rand() % (int)batchEntries; // batchEntries is the size of the batch
	    	randomEvent += lowest_nentry; // shift by the lowest entry of the batch
	    	selectRandomIdxToSave.insert( randomEvent );
	    }
            cout << "randomly selected some events to save" << endl;
            cout << "now gathering the events we are told to always save:";
            for (int i=0; i<(int)parseEventsToSave.varStringSet.size(); ++i){
                cout << " " << parseEventsToSave.varStringSet[i] << endl;
                selectRandomIdxToSave.insert(stoi(parseEventsToSave.varStringSet[i]));
            }
            cout << "\n" << endl;
        }
        else {
            saveAllHistograms=true;
        }

        // START FIT MANAGER
        string sThread = to_string(iProcess);
        fitManager fm(sThread);
        // ---------------------------

        // IF WE WANT TO BOOTSTRAP WE NEED TO SAVE THE DATA
        if (nBS<0) {
            cout << "nBS must be >= 0. Set to 0 if you do not want to do any bootstrapping!" << endl;
            exit(0);
        }
        vector<float> qvalues; qvalues.reserve(nBS);

	// the main loop where we loop through all events in a double for loop to calculate dij. Iterating through all j we can find the k nearest neighbors to event i.
	// Then we can plot the k nearest neighbors in the discriminating distribution
        // Finally calculate the q-value by fitting and getting signal fraction
	//logFile << std::fixed << std::setprecision(6);
	int randomEntry;
        int counter=-1;
        RooTrace::mark();
        for (int ientry=lowest_nentry; ientry<largest_nentry; ientry++){
                fm.rooData->reset();
                dHist_qvaluesBS->Reset();
                if (saveBranchOfNeighbors){
                    memset(neighbors,-1,sizeof(neighbors)); // last argument sets that number of bytes to the specified value. Dont want just kDim here since it is 4 Bytes per
                }
                qvalues.clear();
        	if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end() || saveAllHistograms) {
                    qHistsFile = new TFile((cwd+"/histograms"+runTag+"/"+fileTag+"/qValueHists_"+to_string(ientry)+".root").c_str(),"RECREATE");
                }

                // Outputting the progress of each thread
		flatEntryNumber=ientry;
                if (ientry==lowest_nentry)
                    start = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
		auto duration_beginEvent = std::chrono::high_resolution_clock::now();
		if(saveEventLevelProcessSpeed) { logFile << "Starting event " << ientry << "/" << largest_nentry << " ---- Global Time: " << duration << "ms" << endl; }
		cout << "Starting event " << ientry << "/" << largest_nentry << " ---- Global Time: " << duration << "ms" << endl; 
		
                // Phase space Definitions
		for ( int iVar=0; iVar<phaseSpaceDim; ++iVar ){
			phasePoint1[iVar] = phaseSpaceVars[iVar][ientry];
		}
                int nPotentialNeighbors=(int)phasePoint2PotentialNeighbor.size();
                vector<int> phasePoint2PotentailNeighbor_BS;
                phasePoint2PotentailNeighbor_BS.reserve(nPotentialNeighbors);

                // Grab current discriminating variable values
                for ( int iVar=0; iVar<discrimVarDim; ++iVar)
                    currentValues[iVar] = discrimVars[iVar][ientry];

                // --------------------------------------------
                // Random generator for resampling of the input data, in this case the input data will be the set of potential neighbors for phasePoint2.
                // --------------------------------------------
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
                // -----------------------------------------------------------------------------------------------------------------------------
                // RNG using distributions are really cheap. It is also threadsafe also whereas rand() is not. 
                //      You will see all the threads are not operating at max potential if you using rand
                // https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
                // Mersenne twistor requires large state storage but has the highest quality + period. 
                // Lagged Fibonoci Generator liek ranlux24 has a much smaller period but should be much faster
                // https://books.google.com/books?id=zi08DQAAQBAJ&pg=PA237&lpg=PA237&dq=ranlux24+%22period%22&source=bl&ots=N4HXolKGDr&sig=ACfU3U1JfZcQczw-xoAv-lyHTOCOf4Xpjw&hl=en&sa=X&ved=2ahUKEwiosZP8-_npAhW2RTABHatDDe8Q6AEwAnoECAsQAQ#v=onepage&q=ranlux24%20%22period%22&f=false
                // mt19937 has a period around 10^6000 whereas ranlux24 has a period around 10^171. ranlux24 produces 24 bit integers. 
                //      24 bits ~ 17M which is much larger than the number of potential neighbors. so in our case ranlux is probably just fine
                // -----------------------------------------------------------------------------------------------------------------------------
                //if we ask for a range of 0 to 10 it will include 10. phasePoint2PotentialNeighbor has a zero index. Need to subtract by 1
                std::uniform_int_distribution<int> distribution(0,nPotentialNeighbors-1); 
                static thread_local std::ranlux24_base generator(seed); 

                // ---------------------
                // NOW FIND NEIGHBORS AND EXTRACT Q-FACTORS.
                // - The last itertion is also the full data. This is because we will draw the histograms on the last iteration to 
                //   skip intermediate saving of the histograms
                //  if nBS > 0 then we will rerun and resample to get the bootstrapped q-factors
                // ---------------------
                for (int iBS=0; iBS<nBS+1; ++iBS){ 
		    bestNLL=DBL_MAX;
		    worstNLL=-1*DBL_MAX; // DBL_MIN is basically 0. We want a very negative number
		    duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(saveEventLevelProcessSpeed){logFile << "\tBegin bootstrapping potential neighbors: +" << duration << "ms" << endl; }
                    phasePoint2PotentailNeighbor_BS.clear();
                    if (iBS!=nBS){ // resampling neighbors with replacement if we want to do bootstrapping
                        for(int iNeighbor=0; iNeighbor<nPotentialNeighbors; ++iNeighbor){
                            phasePoint2PotentailNeighbor_BS.push_back(phasePoint2PotentialNeighbor[distribution(generator)]);
                        }
                    }
                    else{
                        phasePoint2PotentailNeighbor_BS = phasePoint2PotentialNeighbor;
                    }
		    duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(saveEventLevelProcessSpeed){logFile << "\tBegin finding neighbors: +" << duration << "ms" << endl; }

                    // What if we wanted to look at k random neighbors?
                    if (doKRandomNeighbors){
                        cout << "\tDoing k random neighbors" << endl;
		        for (int jentry=0; jentry<kDim;++jentry) {  
		              randomEntry = rand() % nentries;
		              distKNN.insertPair(make_pair(1, randomEntry) ); //just using 1 as a distance. Doesnt matter anyways
		        }
                    }
                    else {
                        cout << "\tFinding nearest neighbors" << endl;
		        for (int jentry : phasePoint2PotentailNeighbor_BS) {  
                            if (jentry == ientry){ continue; } 
		            for ( int iVar=0; iVar<phaseSpaceDim; ++iVar ){
		               	phasePoint2[iVar] = phaseSpaceVars[iVar][jentry];
		            }
		            distance = calc_distance(phaseSpaceDim,phasePoint1,phasePoint2,verbose_outputDistCalc);
		            if ( verbose_outputDistCalc ) { cout << "event (i,j)=(" << ientry << "," << jentry << ") has distance=" << distance << endl;} 
		            distKNN.insertPair(make_pair(distance,jentry));
		        }
                    }
		    duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(saveEventLevelProcessSpeed){logFile << "\tFound neighbors: +" << duration << "ms" << endl; }
		    if (distKNN.kNN.size() != kDim){ cout << "size of distKNN is not equal to kDim! size,kDim="<< distKNN.kNN.size() << "," << kDim 
		        << "\n    -- if size is 1 less than kDim it is probably because kDim=nentries and event i cannot be a neighbor to itself" 
                        << "\n    -- if size != kDim it could also mean that the number of spectroscopically unique neighbors reduces the number of poential neighbors below kDim" << endl;}
                    if(verbose_outputDistCalc){ cout << "These are our neighbors" << endl; }

                    int iNeighbor=-1;
		    while ( distKNN.kNN.empty() == false ){
		        newPair = distKNN.kNN.top();
		        distKNN.kNN.pop();
                        weight=accWeights[newPair.second];

                        for ( int iVar=0; iVar<discrimVarDim; ++iVar)
                            neighborValues[iVar] = discrimVars[iVar][newPair.second];
                        bool keptNeighbor=fm.insert(neighborValues,weight);

                        if (keptNeighbor && saveBranchOfNeighbors){
                            neighbors[++iNeighbor]=newPair.second;
                        }
                        if (verbose_outputDistCalc){
                            if (discrimVarDim==1)
		                cout << "(distance, idx)=(" << newPair.first << ", " << newPair.second << ") with x = " << discrimVars[0][newPair.second] << 
                                    " and weight= " << weight << endl; 
                            if (discrimVarDim==2)
		                cout << "(" << newPair.first << ", " << newPair.second << ") with x,y = " << discrimVars[0][newPair.second] 
                                     << ", " << discrimVars[1][newPair.second] << " and weight= " << weight << endl; 
                        }
		    }
		    
		    duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		    if(saveEventLevelProcessSpeed){logFile <<	"\tFilled neighbors: +" << duration << "ms" << endl;}
		    
		    // /////////////////////////////////////////
		    // Calcuclate q-value
		    // /////////////////////////////////////////
                    RooArgSet* savedParams; 
                    RooArgSet* params; // intermedate parameter values for the pdfs
		    for ( auto initSigFrac : sigFracs ){
		        duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(saveEventLevelProcessSpeed){logFile <<	"\tPrepping 2D fits:	 +" << duration << "ms" << endl;}
                        // intitialize the variables for the fit PDF
                        fm.reinitialize(initSigFrac);

		        duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(saveEventLevelProcessSpeed){logFile <<	"\tBeginning Fit:	 +" << duration << "ms" << endl;}
                        ////////////////////////////////////////////////////////////////////////////////////////////
                        //  BEGIN FITTING
                        ////////////////////////////////////////////////////////////////////////////////////////////
                        gSystem->GetProcInfo(&pinfo);
                        if (saveMemUsage)
                            outputMemUsage(pinfo,"\tBefore fitTo: ");
                        NLL = fm.fit();
                        //status = 1    : Covariance was made pos defined
                        //status = 2    : Hesse is invalid
                        //status = 3    : Edm is above max 
                        //status = 4    : Reached call limit
                        //status = 5    : Any other failure
                        fitStatus=fm.fitStatus;

                        params=fm.getParameters();
                        params->Print("v");
                        gSystem->GetProcInfo(&pinfo);
                        if (saveMemUsage)
                            outputMemUsage(pinfo,"\tAfter fitTo: ");
                        
		        duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(saveEventLevelProcessSpeed){logFile <<	"\tCompleted Fit: +" << duration << "ms" << endl;}
                        
                        // setting parameters for bkg/sig and extracting q-value
                        qvalue = fm.calculate_q(currentValues);
                        
                        if ( isnan(qvalue) ){
                            params=fm.getParameters();
                            params->Print("s");
                            cout << "QVALUE IS NAN!" << endl;
                        }

		        duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		        if(saveEventLevelProcessSpeed){logFile <<	"\tExtracted q-value (" << qvalue << "): +" << duration << "ms" << endl;}
                        
		    	if (NLL < bestNLL){
		    		best_qvalue = qvalue;
                                best_fitStatus = fitStatus;
                                best_nsig = fm.nsig->getVal();
                                best_nbkg = fm.nbkg->getVal();
                                best_ntot = best_nsig + best_nbkg;
                                eff_nentries=fm.eff_nentries;
                                effNentriesMinusTotal = eff_nentries-best_ntot;
                                //cout << effNentriesMinusTotal << ", " << best_nsig << ", " << best_nbkg
                                //     << ", " << nsig.getError() << ", " << nbkg.getError() << 
                                //    ", " << best_ntot << ", " << eff_nentries << endl;
		    		bestNLL=NLL;
                                params=fm.getParameters();//fm.rooSigPlusBkg->getParameters(RooArgList(*fm.roo_Mpi0,*fm.roo_Meta));
                                savedParams = (RooArgSet*)params->snapshot();
		    	} 
		    	if (NLL > worstNLL){
		    		worstNLL = NLL;
                                worst_qvalue = qvalue;
		    	}
                        //delete roo_result;
                    
                        // Different ways to calculate chiSq: https://nbviewer.jupyter.org/gist/wiso/443934add13fd7226e4b
                        // We can calculate chiSq by binning the data and using RooChi2Var. Example calculation is at:
                        // https://hep.lancs.ac.uk/~ajf/root/RooChi2MCSModule_8cxx_source.html
                        // Think we want to use nDataPts - nConstraints as in https://ned.ipac.caltech.edu/level5/Leo/Stats7_2.html
                        //RooDataHist* binnedData = rooData.binnedClone();
                        //RooChi2Var chiSq("chiSq","chiSq",rooSigPlusBkg,*binnedData);//,DataError(RooAbsData::SumW2));
                        //RooPlot* roo_Mpi0_Meta_frame = new RooPlot(roo_Mpi0,roo_Meta);
                        //rooData.plotOn(roo_Mpi0_Meta_frame);
                        //rooSigPlusBkg.plotOn(roo_Mpi0_Meta_frame);
                        //double chiSq = roo_Mpi0_Meta_frame->chiSquare();
                        //RooArgSet* floatPars = (RooArgSet*)rooSigPlusBkg.getParameters(rooData)->selectByAttrib("Constant",kFALSE);
                        //int nParams = floatPars->getSize();
                        //int nBins = 900;
                        //double chiSqPerDOF = chiSq/(nBins-nParams);
		    	// now that the q-value is found we can get the NLL and save the parameters with the best NLL
                        // for more information look at RooFitResult class https://root.cern.ch/doc/master/classRooFitResult.html
                        //cout << "NLL, chiSq, reduced ChiSq, nParams: " << NLL << ", " << chiSq << ", " << chiSqPerDOF << ", " << nParams << endl;
		    } // close redistribute initialization fit loop

                    // Save the qvalues for this iteration
                    if(iBS<nBS){  
                        qvalues.push_back(best_qvalue);
                        dHist_qvaluesBS->Fill(best_qvalue);
                    }

                    // EITHER WE POTENTIALLY SAVE THE HISTOGRAM AFTER ALL BS IS DONE OR WE SAVE IT EVERY ITERATION DURING THE BS. ALLOWS US TO SEE HOW THINGS DEVELOP
                    // THERE IS STILL AN INNER CONDITION THAT ONLY SELECTS A SUBSET OF THESE. SINCE WE DONT WANT 100K+ HISTOGRAMS
                    if (iBS==nBS || saveBShistsAlso){
        	        // /////////////////////////////////////////
        	        // Drawing histogram 
        	        // /////////////////////////////////////////
        	        if ( selectRandomIdxToSave.find(ientry) != selectRandomIdxToSave.end() || saveAllHistograms) {
        	            // Here we draw the histograms that were randomly selected
        	            allCanvases->Clear();

		            duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - duration_beginEvent).count();
		            cout <<	"\tSaving diagnostic histogram: " << duration << "ms" << endl;
		            if(saveEventLevelProcessSpeed){logFile <<	"\tSaving histogram: +" << duration << "ms" << endl;}
                            /// -- Might need to use this code section if RooFit ever becomes thread safe. 
                            //         If we draw rooSigPlusBkg with the various Components we get an error when using plotOn
                            ///        when using multi threads. Projecting each sub PDF makes it work fine. 
                            //RooAbsPdf* rooGaus2DProjMpi0 = rooGaus2D.createProjection(roo_Mpi0);
                            //rooGaus2DProjMpi0->setNormRange(("roo_fitRangeMpi0"+sThread).c_str());
                            //rooGaus2DProjMpi0->plotOn(roo_Meta_frame);
                            //RooAbsPdf* rooBkgProjMpi0 = rooBkg.createProjection(roo_Mpi0);
                            //rooBkgProjMpi0->setNormRange(("roo_fitRangeMpi0"+sThread).c_str());
                            //rooBkgProjMpi0->plotOn(roo_Meta_frame, LineStyle(kDashed),LineColor(kOrange));
                            //RooAbsPdf* rooSigPlusBkgProjMpi0 = rooSigPlusBkg.createProjection(roo_Mpi0);
                            //rooSigPlusBkgProjMpi0->setNormRange(("roo_fitRangeMpi0"+sThread).c_str());
                            //rooSigPlusBkgProjMpi0->plotOn(roo_Meta_frame);
                            //rooSigPlusBkgProjMpi0->plotOn(roo_Meta_frame, Components("rooBkg*"),LineStyle(kDashed),LineColor(kOrange));
                            //rooSigPlusBkgProjMpi0->paramOn(roo_Meta_frame);
                            // OMG FOUND HOW TO FIX THE BUG! WE NEED TO SET THE NORMRANGE AND RANGE HERE IF WE ARE GOING TO FIT 
                            // USING THE SAME PDF ON THE SAME DATA MULTIPLE TIMES. THERE IS AN ERROR ABOUT THE NORMALIZATIONS 
                            // AND SOMEHOW MORE FUNCTIONS ARE CREATED AND WE GET A BUNCH OF SAME OF THE SAME OBJECTS THAT GETS ADDED

                            // reload the best params
                            params=fm.getParameters();//fm.rooSigPlusBkg->getParameters(RooArgList(*fm.roo_Mpi0,*fm.roo_Meta));
                            *params = *savedParams;
                            //params->Print("s");
                            fm.NLL=bestNLL;

                            /////////////////////////////////////////////////////////////////////////
                            ////////////////////////////////////// BOOTSTRAP HISTOGRAM OF Q-FACTORS
                            if (iBS==nBS){
                                qvalueBS_std=calculateStd(nBS,&qvalues[0]);
                                dHist_qvaluesBS->SetTitle(("STD: "+to_string(qvalueBS_std)).c_str());
                            }
                            else{
                                dHist_qvaluesBS->SetTitle("Bootstrapped Q-Factors");
                            }
                            /////////////////////////////////////////////////////////////////////////
                            /////////////////////////////////////////////////////////////////////////

                            fm.drawFitPlots(currentValues, ientry, dHist_qvaluesBS, 
                                            best_qvalue, iBS, allCanvases, qHistsFile);
                            
        	            if(saveEventLevelProcessSpeed){
        	                 duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() 
                                             - duration_beginEvent).count();
        	                 logFile << "\tSaved this histogram since it was randomly selected: +" << duration <<  "ms" << endl;
        	            }
		        } 
                    } // closes condition to save histograms
                    delete params;
                    delete savedParams;
                } // finishes nBS loop

		if(saveEventLevelProcessSpeed){logFile << "\tCurrent Best NLL = " << to_string(bestNLL) << ": " << duration << "ms" << endl; }
		resultsTree->Fill();

                gSystem->GetProcInfo(&pinfo);
                if (ientry==lowest_nentry){
                    startMemResident=pinfo.fMemResident;
                }
                if (ientry==(largest_nentry-1)){
                    endMemResident=pinfo.fMemResident;
                }
                if (saveMemUsage)
                    outputMemUsage(pinfo,"\tEnd of current iteration - mem usuage: ");
	}
        resultsFile->cd();
        resultsTree->Write();
        int deltaMem=endMemResident-startMemResident;
        float avgMem=(float)deltaMem/(largest_nentry-lowest_nentry);
        cout << endl << "Memory increase = " << deltaMem << "KB" << endl;
        cout << "Avg increase per iteration = " << avgMem << "KB" << endl;
        
	// Finish the log files by including an elapsed time and finally closing the file
	if (saveEventLevelProcessSpeed){
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
		logFile << "Total time: " << duration << " ms" << endl;
		logFile << "Time Per Event: " << duration/(largest_nentry-lowest_nentry) << " ms" << endl;
	}
	logFile.close();

        // OLD CODE FOR SPAWNING THREADS INSTEAD OF WHAT WE USE NOW, WHICH IS JUST PROCESSES
	//};
        //// Now that we have the lambda function we can start to spawn threads
	//cout << "Launching " << nProcess << " threads 1 second apart!" << endl;
	//vector<thread> threads;
	//for ( int iProcess=0; iProcess<nProcess; ++iProcess){
	//	cout << "(Thread " << iProcess << ") is starting" << endl;
	//	threads.emplace_back( [f, iProcess] { f(iProcess); } );
        //        sleep(1);
	//	//threads[iProcess] = std::thread(QFactorAnalysis::staticEntryPoint, this, iProcess);
	//}
	//for (auto&& t : threads) t.join(); // join waits for completion
	//cout << "Threads have completed running!" << endl;
}

int main( int argc, char* argv[] ){
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // For thread safety we need this
	TH1::AddDirectory(kFALSE);

	// This suppresses all the "info" messages from roofit. Like saying that it will use numeric integrator for the generic pdf we define
	// https://root-forum.cern.ch/t/suppressing-info-messages/14642/6
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
        RooTrace::active(kTRUE);
    	gStyle->SetOptFit(1111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);

	// -----------------------------------------------------
	// -----------------------------------------------------
	//                          PARSE THE COMMAND LINE ARGUMENTS
	// -----------------------------------------------------
	// -----------------------------------------------------
        int iProcess=atoi(argv[1]);

	QFactorAnalysis analysisControl;
	analysisControl.initialize(rootFileLoc, rootTreeName);
	analysisControl.loadData();
        analysisControl.runQFactorThreaded(iProcess);
        return 0;
}


