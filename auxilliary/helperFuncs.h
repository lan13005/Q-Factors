#ifndef HELPERFUNCS_H
#define HELPERFUNCS_H

#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <queue>
#include <TTree.h>
#include <TLeaf.h>

using namespace std;

// ------------------------------------------------------
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
#define MAXVALUE 1
std::vector<float> binRangeEta={200,0.34,0.8};
std::vector<float> binRangePi0={200,0.075,0.21};
float binWidthEta=(binRangeEta[2]-binRangeEta[1])/binRangeEta[0];
float binWidthPi0=(binRangePi0[2]-binRangePi0[1])/binRangePi0[0];
std::vector<float> fitRangeY={0.36,0.75};
std::vector<float> fitRangeX={0.085,0.185};
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// DO NOT MODIFY THIS SECTION. GETS OVERWRITTEN BY RUN.PY
// ------------------------------------------------------


// The following two classes will be used to keep track of pairs of distances and index, sorted by distance. priority_queue in stl requires three arguments
// which are type, container for type, and a comparator. we will use pair as the type which is held in a vector container. compareDist is the comparator
// which compares the first elements of two pairs. The first element is the distance, the second is the index j (in dij when calculating the distances). 
// distSort_kNN will setup our priority_queue that keeps a maximum of kDim elements simply by popping and pushing data. 
// This type of sorting should have k*log(k) sorting, I think. If we do this N times then the complexit ~ N*k*log(k). Probably have to check me on this
class compareDist
{
    public:
        bool operator() (const std::pair<float,int>& p1, const std::pair<float,int>& p2)
        {
            // < for min heap
            // maybe > for max heap?
            return p1.first < p2.first;
        }
};
class distSort_kNN
{
    public:
        // constructor with basic output and setting kDim
        distSort_kNN ( UInt_t kDim ) {
            //std::cout << "Priority Queue is set up to sort the {distance,index} pairs keeping " << kDim << " neighbors" << std::endl;
            _kDim = kDim;
        }

        std::priority_queue <std::pair<float,int>, vector<std::pair<float,int>>, compareDist > kNN;
        
        // depending on the size and the distanceis we will push, or pop+push
        void insertPair ( std::pair<float,int> newPair ){
            if ( kNN.size() >= _kDim ){
                // > for min heap and maybe < for max heap 
                if (kNN.top().first > newPair.first) {
                    //std::cout << "\nPOPPING OUT " << kNN.top().first << std::endl;
                    kNN.pop();
                    kNN.push( newPair );
                }
            } 
            else {      
                //std::cout << "\nPUSHING " << newPair.first << std::endl;
                kNN.push( newPair );
            }
        }

    private:
        UInt_t _kDim;
        std::pair<float,int> _pair;
};

// This class will be used to standardize the phase space variable. 
// Can either do stddev or range standardization
class standardizeArray{
	public:
    		float max_inputVector;
         	float min_inputVector;
		
		void rangeStandardization(std::vector<float> &inputVector, int varIdx, long long nentries, int nbranches){
    			max_inputVector = DBL_MIN;
         		min_inputVector = DBL_MAX;
                        std::cout << "\nStarting range std" << std::endl;
			for (int ientry=0; ientry<nentries; ++ientry){
				if (inputVector[ientry*nbranches+varIdx] > max_inputVector){
					max_inputVector = inputVector[ientry*nbranches+varIdx];
				}
				if (inputVector[ientry*nbranches+varIdx] < min_inputVector){
					min_inputVector = inputVector[ientry*nbranches+varIdx];
				}
			}
			for (int ientry=0; ientry<nentries; ++ientry){
				inputVector[ientry*nbranches+varIdx] = (inputVector[ientry*nbranches+varIdx]-min_inputVector)/(max_inputVector-min_inputVector);
			}
                        std::cout << "Max,min: " << max_inputVector << "," << min_inputVector << std::endl;
                        std::cout << "--Finished Range standardizing " << std::endl;
		}
		
		float calcStd(std::vector<float> &inputVector, int varIdx, long long nentries, int nbranches){
			float local_std=0;
			float diff=0;
			float mean=0;
			for (int ientry=0; ientry<nentries; ++ientry){
			    mean+=inputVector[ientry*nbranches+varIdx];
			} 
			mean/=nentries;
			for (int ientry=0; ientry<nentries; ++ientry){
                                //std::cout << "mean: " << mean << std::endl;
				diff = (inputVector[ientry*nbranches+varIdx]-mean);
                                //std::cout << "diff: " << diff << std::endl;
				local_std += diff*diff;
			}
			local_std /= nentries;
                        std::cout << "STD: " << local_std << std::endl;
                        std::cout << "--Finished Std standardizing " << std::endl;
			return sqrt(local_std);
		}
		void stdevStandardization(std::vector<float> &inputVector, int varIdx, long long nentries, int nbranches){
			float std = calcStd(inputVector, varIdx, nentries, nbranches);
			for (int ientry=0; ientry < nentries; ++ientry){
				inputVector[ientry*nbranches+varIdx] = inputVector[ientry*nbranches+varIdx]/std; 
			} 
                        std::cout << "Finished Stdev Standardization" << endl;
		}
};

// our distance calculation between two phase points
float calc_distance( int dim, float* phaseSpace_1, float* phaseSpace_2, bool verbose_outputDistCalc){
	float sum = 0;
	float diff=0;
        if(verbose_outputDistCalc){std::cout << "New event, new sum = " << sum << std::endl;}
	for (int i=0; i<dim; ++i){
		diff = phaseSpace_1[i]-phaseSpace_2[i];
		sum += abs(diff);//diff*diff;
                if(verbose_outputDistCalc){
                    std::cout << "phasePoint1["<<i<<"]="<<phaseSpace_1[i]<<", phasePoint2["<<i<<"]="<<phaseSpace_2[i]<<" --- abs distance=" << abs(diff) << "---- total so far="<<sum<<std::endl;
		}
	}
	return sum;
}


// Class to parse the string of phase space variable to consider
class parseVarString{
    private:
        int nVars=0;
        std::string delimiter=";";
        std::string token;
        std::string varString;
        size_t pos=0;
    public:
        std::vector<std::string> varStringSet;
        parseVarString(){}; 
        void parseString(std::string inputString){
            varString =  inputString;
            varStringSet.clear();
            if (varString!=""){
                while ((pos = varString.find(delimiter)) != std::string::npos) {
                    token = varString.substr(0, pos);
                    varStringSet.push_back(token);
                    //std::cout << token << std::endl;
                    varString.erase(0, pos + delimiter.length());
                    ++nVars;
                }
                varStringSet.push_back(varString);
                //std::cout << varString << std::endl;
            }
        }
};

float calculateStd(int nentries, float* input){
    float mean=0;
    for(int ientry=0; ientry<nentries; ++ientry){ 
        mean += input[ientry];
    }
    mean /= nentries;
    float diff;
    float std=0;
    for(int ientry=0; ientry<nentries; ++ientry){ 
        diff = input[ientry]-mean;
        std += diff*diff;
    }
    std /= nentries-1;
    return sqrt(std);
}

// Not used anymore. I thought it would interseting to check the stdev of the k nearest neighbors. 
// Used to calculate the current std as we stream/insert more data
class cumulativeStd{
    public:
        std::vector<float> inputVector;
        // kDim would be the typical size of the calculation. But sometimes we will choose nentries < kDim when testing quickly (this is never the case in real example though)
        cumulativeStd ( int kDim ){ _kDim=800; }
        void insertValue ( float value ) {
            inputVector[_timesCalled] = value;
            ++_timesCalled;
        }
        float calcStd(){
            for (int i=0; i<_timesCalled; ++i){
                _sum+=inputVector[i];
            }
            _sum /= _timesCalled;
            _mean=_sum; _sum=0;
            //std::cout << "mean: " << _mean << std::endl;
            for (int i=0; i<_timesCalled; ++i){
                _diff = (inputVector[i]-_mean);
                _sum += _diff*_diff;
                //std::cout << "sum: " << _sum << std::endl;
            }
            _sum /= _timesCalled;
            return sqrt(_sum);
        }


    private:
        int _timesCalled=0;
        float _sum=0;
        float _mean=0;
        float _diff;
        UInt_t _kDim;

};

void outputMemUsage(ProcInfo_t pinfo, string contextString){
    gSystem->GetProcInfo(&pinfo);
    cout << contextString << pinfo.fMemResident << "KB" << endl; 
}

string setBranchAddress(TTree* tree, string variable, Long64_t* value_i, float* value_f, double* value){
    //tree->Print();
    cout << "setBranchAddress: " << variable << endl;
    string typeName=tree->GetLeaf(variable.c_str())->GetTypeName();
    if (typeName=="Long64_t"){
        cout << "Setting Branch with data type Long64 " << variable << endl; 
        tree->SetBranchAddress(variable.c_str(),value_i);
    }
    else if (typeName=="Float_t"){
        cout << "Setting Branch with data type float " << variable << endl; 
        tree->SetBranchAddress(variable.c_str(),value_f);
    }
    else if (typeName=="Double_t"){
        cout << "Setting Branch with data type double " << variable << endl; 
        tree->SetBranchAddress(variable.c_str(),value);
    }
    else { cout << "Branch " << variable << " has unlisted datatype: " << typeName << ". exiting" << endl; exit(0); }
   return typeName; 
}

// Find the intersection of two vector of strings
std::vector<std::string> intersection(std::vector<std::string> &v1,
                                      std::vector<std::string> &v2){
    std::vector<std::string> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}

// Get index of string in vector<string>
int getIndex(vector<string> v, string K)
{
    auto it = find(v.begin(), v.end(), K);
    if (it != v.end())
        return it-v.begin();
    else
        return -1;
}

string replace_str(string s, string search_str, string replace_str){
    std::string::size_type n = 0;
    while ( ( n = s.find( search_str, n ) ) != std::string::npos ){
        s.replace( n, search_str.size(), replace_str );
        n += replace_str.size();
    }
    return s;
}

vector<bool> checkConditions(vector<int> neighbors, vector<float> values, map<string,int> &nameToIdx, int n_branches, vector<string> conditions){
    bool verbose=false;
    
    // The final output: checks if the condition is passed for this neighbor
    vector<bool> isins(neighbors.size(),1);

    for (auto condition: conditions){
        // The two conditions we will check for 
        bool foundLT = condition.find("<") != string::npos;
        bool foundGT = condition.find(">") != string::npos;

        if (verbose) cout << "condition -- passed?" << endl;
        if (foundLT){
            string var=condition.substr(0,condition.find("<"));
            float threshold=stof(condition.substr(condition.find("<")+1, condition.length()));
            if ( nameToIdx.find(var)==nameToIdx.end() ){
                cout << "checkConditions function could not find " << var << " branch in nameToIdx. exiting..." << endl;
                exit(0);
            }
            int varIdx = nameToIdx[var];
            for (int i=0; i<neighbors.size(); ++i){
                bool passed = values[i*n_branches+varIdx]<threshold;
                isins[i] = isins[i]*passed;
                if (verbose) cout << var << "=" << values[neighbors[i]*n_branches+varIdx] << "<" << threshold << " -- " << isins[i] << endl;
            }
        }
        else if(foundGT){
            string var=condition.substr(0,condition.find(">"));
            float threshold=stof(condition.substr(condition.find(">")+1, condition.length()));
            if ( nameToIdx.find(var)==nameToIdx.end() ){
                cout << "checkConditions function could not find " << var << " branch in nameToIdx. exiting..." << endl;
                exit(0);
            }
            int varIdx = nameToIdx[var];
            for (int i=0; i<neighbors.size(); ++i){
                bool passed = values[neighbors[i]*n_branches+varIdx]>threshold;
                isins[i] = isins[i]*passed;
                if (verbose) cout << var << "=" << values[i*n_branches+varIdx] << ">" << threshold << " -- " << isins[i] << endl;
            }
        }
        else{
            cout << "Received unexpected comparator in neighborReqs! Only accepts {<,>} signs. exiting..." << endl;
            exit(0);
        }
    }
    return isins;
}

#endif
