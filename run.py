import subprocess
import os
import sys
import time
from itertools import combinations
from termcolor import colored

start_time = time.time()
#############################################################################
###################  DEFINING ENVIRONMENT VARIABLES #########################
#############################################################################
### NOTE: The code calculates the q-factor for all events in the branch, make sure you feed only those events of interest
## --- STANDARD ----
# rootFileLocs: Root file to be analyzed, followed by the name of the tree and a name tag to avoid overwriting runs
# fitWeights: semicolon separated names of the branches for fit weighting. Fitted histograms will be weighted by the product of the branch values. "none" to set weights to 1
# sigWeights: semicolon separated names of the branches to draw with the makePlots program. These are the "Q" plots. Q_tot will be weighted by the product of the branch values in sigWeights. Q_sig will be weighted by sigWeights AND the newly computed q-factor. Q_tot is weighted by sigWeights. "none" to set weights to 1
# altWeights: semicolon separted names of the branches to draw with the makePlots program. These are the "SB" plots. SB_sig will be weighted by the product of the branch values in altWeights. SB_tot is not weighted. "none" to set weights to 1
# varStringBase: semicolon separated branch names to get phase space variables for distance calculation
# discrimVars: semicolon separated branch names to get discriminating/reference variables
# extraVars: semicolon separated branch names to get extra variables. One use for this is to load extra variables to make a selection on the set of neighbors (i.e. fitrange). Then the variable can be used in neighborReqs, for instance.
# neighborReqs: requirements for the neighbors. Useful to set a fit range. This can be useful if you are trying to compare to sideband subtraction where sideband subtraction might require a larger range than you want q-factors to consider
# nProcess: Number of processes to spawn
# kDim: Number of neighbors
# nentries: Number of combos to run over. Set to -1 to run over all. This should be significantly larger than kDim or we might get errors.
# numberEventsToSavePerProcess: how many event level fit histograms (root files) we want to save. -1 = Save all histograms
## --- ADVANCED ----
# standardizationType: {range,std} what type of standardization to apply when normalizing the phase space variables. Do nothing if any other string 
# redistributeBkgSigFits: {1,0} should we do the 3 different fits where there is 100% bkg, 50/50, 100% signal initialization. 1 to do all 3 and 0 to do only 50/50 
# nRndRepSubset: size of the random subset of potential neighbors. If <=0 or > nentries then we will not consider random subsets. This flag will limit the number of neighbors to consider. The usefulness of this flag will depend on the density of the phase space. If there is sufficient density then it can be argued that the set of potential neighbors can be reduced. 
# doKRandomNeighbors: should we use k random neighbors as a test instead of doing k nearest neighbors? Useful check to decouple fit performance and neighbor finding
# nBS: number of times we should bootstrap the set of neighbors to calculate q-factors with. Used to extract an error on the q-factors. 0 = no bootstrapping
# runTag: 3 folders are outputs of this set of programs {fitResults/diagnosticPlots/histograms}. Append a runTag to the names - i.e. fitResults_newTag to prevent overwriting
# seedShift: in case we dont want to save the same q-value histogram we can choose another random seed. This flag is related to the numberEventsToSavePerProcess flag.
# saveBShistsAlso: should we save every bootstrapped histogram also?
# alwaysSaveTheseEvents: A histogram of the fit will always be saved for these semicolon separated events which allows for studying a single fit
# saveBranchOfNeighbors: Should we save a branch containing all the neighbors per event. Size increase ~ (4Bytes per int)*kDim*nentries
# saveMemUsage: should we output the memory usage into the logs file? Useful to debug source of memory usuage
# saveEventLevelProcessSpeed:  (LEAVE ON) include info on process speed into processLogX.log files. 
# emailWhenFinished: we can send an email when the code is finished, no email sent if empty string
# runBatch: partly production ready - (default=0) 0=run on a single computer, 1=submit to condor for batch processing
# runAllPhaseCombos: (bool) whether we should run over all possible subsets of the phase space. Number fits = Sum_i (n choose i); i=1 to n+1; n=size of varStringBase. If you have no idea what variables to use you can just search it all if you have enough time and computing power. 
# extraLibs: Include extra libraries to compile main with. Intended for loading custom PDFs

# -------- STANDARD ---------
rootFileLocs=[
#        ("degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol000_090.root",
#            "degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat", "000_090")
#        ,("degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol045.root",
#            "degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat", "045")
#        ,("degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol090.root",
#            "degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat", "090")
#        ,("degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol135.root",
#            "degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat", "135")
#        ,("degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_polAMO.root",
#            "degALL_data_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat", "AMO")

        ("degALL_data_phase1_mEllipse_8288_tLT1_treeFlat_DSelector.root",
            "tree_4g_flat", "phase1")

#        ("zb1_plus_etapi_as_4g_dataset/b1_and_etapi_mEllipse_8288_chi13_tpLT05_omegacut_treeFlat_subset.root",
#            "tree_4g_flat", "flatEtapi_b1_Mpi0VsMeta")
#        ("/d/grid17/ln16/q-values-3/logs/flatEtapi_b1_Mpi0_1111/postQVal_flatTree_flatEtapi_b1_Mpi0_1111.root",
#            "tree_4g_flat", "flatEtapi_b1_Mpi0_MetaIndep")
#        ("zb1_plus_etapi_as_4g_dataset/b1_and_etapi_mEllipse_8288_chi13_tpLT05_omegacut_treeFlat_subset.root",
#            "tree_4g_flat", "test")

#        ("degALL_data_2017_mEllipse_8288_tLT1_chi13_omegacut_treeFlat_DSelector.root", 
#            "tree_4g_flat", "2017_2D")
        ]

_SET_fitWeights="AccWeight" 
_SET_sigWeights="AccWeight" # Just for makePlots to draw with appropriate weights
_SET_altWeights="AccWeight;weightBS" # Just for makePlots to draw with appropriate weights 
#_SET_varStringBase="cosTheta_eta_gj;phi_eta_gj;cosTheta_X_cm;Phi;Mpi0eta;Mpi0g3;Mpi0g4;ph124Rest_angle_g34;mandelstam_teta;ph123Rest_angle_g34"#phi_X_lab" 
_SET_varStringBase="cosTheta_eta_gj;phi_eta_gj;Mpi0eta;Mpi0g3" 
_SET_discrimVars="Meta"
_SET_extraVars=""
_SET_neighborReqs=""#"Meta>0.36;Meta<0.75;Mpi0>0.085;Mpi0<0.185"
_SET_nProcess=36
_SET_kDim=300
_SET_nentries=-1
_SET_numberEventsToSavePerProcess=2
# -------- ADVANCED ---------
_SET_standardizationType="range" 
_SET_redistributeBkgSigFits=0 
_SET_nRndRepSubset=0 
_SET_doKRandomNeighbors=0
_SET_nBS=0 
_SET_runTag="" 
_SET_seedShift=1341 
_SET_saveBShistsAlso=0 
_SET_alwaysSaveTheseEvents=""
_SET_saveBranchOfNeighbors=0
_SET_saveMemUsage=1 
_SET_saveEventLevelProcessSpeed=1 #do not turn off 
_SET_emailWhenFinished="" 
_SET_runBatch=0 
_SET_runAllPhaseCombos=0
_SET_extraLibs=["./auxilliary/customPDFs/bivariateGaus/bivariateGaus_cxx.so"]

#############################################################################
###################  DEALING WITH CMDLINE ARGS   #########################
#############################################################################
args=sys.argv
def showHelp():
    print("\n-help\n")
    print("run.py accepts only one argument. A 2 digit binary code is needed:")
    print("1. 0/1 given to run the q factor analysis")
    print("2. 0/1 given to make all the diagnostic histograms applyng the q-factors")
    print("i.e. if we want to run the q-factor analysis without make the graphs then use 10 as the argument")


def parseCmdArgs():
    if(len(args)!=2):
        showHelp()
        exit()
    if(_SET_kDim >= _SET_nentries and _SET_nentries>0):
        print("We cannot have kDim >= nentries. Since the combo in question cannot be a neighbor to itself we will never have kDim neighbors!")
        exit()
    for arg in args[1:]:
        arg=str(arg)
        if(len(arg)==2):
            print("\n--------------")
            if arg[0]=="1":
                _SET_runQFactor=True
                print("Running Q-Factor Analysis")
            elif arg[0]=="0":
                _SET_runQFactor=False
                print("Skipping Q Factor Analysis")
            else:
                showHelp()
                exit()
            if arg[1]=="1":
                _SET_runMakeHists=True
                print("Running Make Diagnostic Hists")
            elif arg[1]=="0":
                _SET_runMakeHists=False
                print("Skipping Make Diagnostic Hists")
            else:
                showHelp()
                exit()
            print("--------------\n")
            return _SET_runQFactor,_SET_runMakeHists
        else:
            showHelp()
            exit()

_SET_runQFactor,_SET_runMakeHists = parseCmdArgs()

#############################################################################
###################  FIRST DEFINE SOME FUNCTIONS   #########################
#############################################################################

varVec=_SET_varStringBase.rstrip().split(";")

def spawnProcessChangeSetting(varName,varValue,fileName,isString):
    '''
    Replace the value of varName to varValue in the file called fileName. Depending on the value type we can include quotes or not
    '''
    if isString:
        sedArgs=["sed","-i",'s@'+varName+'=".*";@'+varName+'="'+varValue+'";@g',fileName]
    else:
        sedArgs=["sed","-i",'s@'+varName+'=.*;@'+varName+'='+str(varValue)+';@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

def changeDim(search,value,fail_if_empty_str):
    dim=len(value.split(";"))
    if value=="":
        if fail_if_empty_str:
            print(value+" cannot be an empty string. Exiting...")
            exit(0)
        else:
            dim=dim-1
    changeDims=["sed","-i","s@const int "+search+"=.*;@const int "+search+"="+str(dim)+";@g","configSettings.h"]
    print(" ".join(changeDims))
    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...

def reconfigureSettings(fileName, combo, _SET_rootFileLoc, _SET_rootTreeName, _SET_fileTag):
    '''
    Settings we need to set in helpFuncs
    '''
    spawnProcessChangeSetting("rootFileLoc",_SET_rootFileLoc,fileName,True)
    spawnProcessChangeSetting("rootTreeName",_SET_rootTreeName,fileName,True)
    spawnProcessChangeSetting("fileTag",_SET_fileTag,fileName,True)
    spawnProcessChangeSetting("runTag",_SET_runTag,fileName,True)
    spawnProcessChangeSetting("cwd",os.getcwd(),fileName,True)
    spawnProcessChangeSetting("standardizationType",_SET_standardizationType,fileName,True)
    spawnProcessChangeSetting("s_discrimVar",_SET_discrimVars,fileName,True)
    spawnProcessChangeSetting("s_extraVar",_SET_extraVars,fileName,True)
    spawnProcessChangeSetting("s_neighborReqs",_SET_neighborReqs,fileName,True)
    spawnProcessChangeSetting("s_fitWeight",_SET_fitWeights,fileName,True)
    spawnProcessChangeSetting("s_sigWeight",_SET_sigWeights,fileName,True)
    spawnProcessChangeSetting("s_altWeight",_SET_altWeights,fileName,True)
    spawnProcessChangeSetting("standardizationType",_SET_standardizationType,fileName,True)
    spawnProcessChangeSetting("alwaysSaveTheseEvents",_SET_alwaysSaveTheseEvents,fileName,True)

    # Need to update the phaseVar to match the phase space combination we are looking at
    selectedVar = [varVec[ele] for ele in combo]
    _SET_varString=";".join(selectedVar)
    spawnProcessChangeSetting("s_phaseVar",_SET_varString,fileName,True)

    spawnProcessChangeSetting("nProcess",_SET_nProcess,fileName,False)
    spawnProcessChangeSetting("kDim",_SET_kDim,fileName,False)
    spawnProcessChangeSetting("ckDim",_SET_kDim,fileName,False)
    spawnProcessChangeSetting("redistributeBkgSigFits",_SET_redistributeBkgSigFits,fileName,False)
    spawnProcessChangeSetting("doKRandomNeighbors",_SET_doKRandomNeighbors,fileName,False)
    spawnProcessChangeSetting("numberEventsToSavePerProcess",_SET_numberEventsToSavePerProcess,fileName,False)
    spawnProcessChangeSetting("seedShift",_SET_seedShift,fileName,False)
    spawnProcessChangeSetting("nentries",_SET_nentries,fileName,False)
    spawnProcessChangeSetting("nBS",_SET_nBS,fileName,False)
    spawnProcessChangeSetting("saveBShistsAlso",_SET_saveBShistsAlso,fileName,False)
    spawnProcessChangeSetting("saveEventLevelProcessSpeed",_SET_saveEventLevelProcessSpeed,fileName,False)
    spawnProcessChangeSetting("saveBranchOfNeighbors",_SET_saveBranchOfNeighbors,fileName,False)
    spawnProcessChangeSetting("saveMemUsage",_SET_saveMemUsage,fileName,False)
    if _SET_nentries == -1:
        _SET_override_nentries=0
    else:
        _SET_override_nentries=1
    spawnProcessChangeSetting("override_nentries",_SET_override_nentries,fileName,False)


def runOverCombo(combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag):
    '''
    This is the main driver section that runs the "main" program with appropriate settings. 
    Here we can accept a combo flag that selects the combination of varString variables to use. This
    is useful if you are trying to do a scan over all possible combinations of phase space variables
    of all set sizes. You should probably set nentries to be smaller or else you will have to wait a long time....
    '''
    # clean up the outputs of the "main" program
    print("Cleaning folders that are used by main")
    os.system("rm -f main")
    os.system("rm -rf logs"+_SET_runTag+"/"+_SET_fileTag)
    os.system("rm -rf histograms"+_SET_runTag+"/"+_SET_fileTag)
    os.system("mkdir -p logs"+_SET_runTag+"/"+_SET_fileTag)
    os.system("cp configSettings.h logs"+_SET_runTag+"/"+_SET_fileTag) 
    os.system("cp configPDFs.h logs"+_SET_runTag+"/"+_SET_fileTag) 
    os.system("cp main.C logs"+_SET_runTag+"/"+_SET_fileTag+"/main.C")
    os.system("cp run.py logs"+_SET_runTag+"/"+_SET_fileTag+"/run.py")
    os.system("mkdir -p histograms"+_SET_runTag+"/"+_SET_fileTag)


    # We use this setup to make it easy to loop through all combinations for phase space variables to determine which set of variables are the best. 
    selectedVar = [varVec[ele] for ele in combo]
    _SET_varString=";".join(selectedVar)
    print("PHASE SPACE VARIABLES: {}".format(_SET_varString))
    numVarInChosenVarStr=len(_SET_varString.split(";"))

        
    changeDim("phaseSpaceDim",_SET_varString,True)
    changeDim("discrimVarDim",_SET_discrimVars,True)
    changeDim("extraVarDim",_SET_extraVars,False)
    changeDim("fitWeightsDim",_SET_fitWeights,False)
    changeDim("sigWeightsDim",_SET_sigWeights,False)
    changeDim("altWeightsDim",_SET_altWeights,False)

#    # Updated dimensionality of phase space and discriminating variables and fitWeights and altWeights
#    changeDims=["sed","-i","s@const int phaseSpaceDim=.*;@const int phaseSpaceDim="+str(numVarInChosenVarStr)+";@g","configSettings.h"]
#    print(" ".join(changeDims))
#    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
#    changeDims=["sed","-i","s@const int discrimVarDim=.*;@const int discrimVarDim="+str(len(_SET_discrimVars.split(";")))+";@g","configSettings.h"]
#    print(" ".join(changeDims))
#    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
#    n_extraVars = len(_SET_extraVars.split(";")) if _SET_extraVars!="" else 0
#    changeDims=["sed","-i","s@const int extraVarDim=.*;@const int extraVarDim="+str(n_extraVars)+";@g","configSettings.h"]
#    print(" ".join(changeDims))
#    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
#    changeDims=["sed","-i","s@const int fitWeightsDim=.*;@const int fitWeightsDim="+str(len(_SET_fitWeights.split(";")))+";@g","configSettings.h"]
#    print(" ".join(changeDims))
#    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
#    changeDims=["sed","-i","s@const int sigWeightsDim=.*;@const int sigWeightsDim="+str(len(_SET_sigWeights.split(";")))+";@g","configSettings.h"]
#    print(" ".join(changeDims))
#    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
#    changeDims=["sed","-i","s@const int altWeightsDim=.*;@const int altWeightsDim="+str(len(_SET_altWeights.split(";")))+";@g","configSettings.h"]
#    print(" ".join(changeDims))
#    subprocess.Popen(changeDims, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...

    # the distance calculation needs dim to know the dimension of phase space. Before we compile it the script needs to know so we have replace before compilation 
    print("Deleting main and recompiling it")
    subprocess.Popen("rm main", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    compileMain=["g++","-o","main","main.C"]
    try:
        rootFlags = subprocess.check_output(["root-config","--cflags","--glibs", "--libs"])
    except:
        raise Exception("ROOT not loaded!")
    rootFlags = subprocess.check_output(["root-config","--cflags","--glibs", "--libs"])
    rootFlags = rootFlags.decode(encoding="utf-8").rstrip().split(" ") #python3 requires decoding first, which python2 doesnt. But it doesn hurt
    rooFitFlags = ["-lRooStats","-lRooFitCore", "-lRooFit"]
    compileMain.extend(rootFlags)
    compileMain.extend(rooFitFlags)
    if ( any([loc=="" for loc in _SET_extraLibs]) ):
       print("Why would you want to import a library with location as an empty string? If you do not want to import any libraries leave _SET_extraLibs=[]")
       print("Exiting...")
       exit(0)
    compileMain.extend(_SET_extraLibs)
    print("\nStarting new compilation\n----------------------")
    print(" ".join(compileMain))
    out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    print(out)
    

    print('Number of threads: '+str(_SET_nProcess))
    if _SET_runBatch==1:
        # ----- BATCH 
        os.system("rm -rf condor"+_SET_runTag)
        for i in range(_SET_nProcess):
            os.system("mkdir -p condor"+_SET_runTag+"/job"+str(i))
        sedArgs=["sed","-i","s/queue.*/queue "+str(_SET_nProcess)+"/g","submit.main"]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        sedArgs=["sed","-i","s@Arguments      = .*@Arguments      =  $(Process)@","submit.main"]
        #sedArgs=["sed","-i","s@Arguments      = .*@Arguments      = "+str(_SET_kDim)+" "+_SET_varString+" "+_SET_standardizationType+" "+_SET_fitLocation+" "+str(_SET_redistributeBkgSigFits)+" "+str(_SET_doKRandomNeighbors)+" "+str(_SET_numberEventsToSavePerProcess)+" $(Process) "+str(_SET_nProcess)+" "+str(_SET_seedShift)+" "+str(_SET_nentries)+" "+str(_SET_nRndRepSubset)+" "+str(_SET_nBS)+" "+str(_SET_saveBShistsAlso)+" "+str(_SET_override_nentries)+" "+str(_SET_saveEventLevelProcessSpeed)+" "+_SET_cwd+"@","submit.main"]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        sedArgs=["sed","-i","s/counts=$((.*-1))/counts=$(("+str(_SET_nProcess)+"-1))/g","submit.sh"]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        os.system("./submit.sh") # call the submit program to launch our condor_submit program
    else:
        # ----- Dumb Multiprocessing 
        _SET_iProcess=0;
        print("Launching processes one second apart...")
        time.sleep(1)
        outLogs=[]
        errLogs=[]
        openProcesses=[]
        pids=[]
        for _SET_iProcess in range(_SET_nProcess):
            print("Launching process "+str(_SET_iProcess))
            outLog = open("logs"+_SET_runTag+"/"+_SET_fileTag+"/out"+str(_SET_iProcess)+".txt","w")
            errLog = open("logs"+_SET_runTag+"/"+_SET_fileTag+"/err"+str(_SET_iProcess)+".txt","w")
            outLogs.append(outLog)
            errLogs.append(errLog)
            executeMain=["./main",str(_SET_iProcess),"&"]
            # print out the commands with the appropriate quotations so it can be directly run if needed
            validCommand=executeMain[:-1]
            for i in range(1,len(validCommand)):
               validCommand[i]='"'+validCommand[i]+'"' 
            validCommand=" ".join(validCommand)
            print(validCommand)
            outLog.write("\n--------------------\nCMD:\n"+validCommand+"\n--------------------")
            openProcess = subprocess.Popen(executeMain,stdout=outLog,stderr=errLog)
            pids.append(openProcess.pid)
            openProcesses.append(openProcess)

        # wait for the progress checks to fully complete and have found completion before querying the exit codes of the `main` program
        #   The time the completion of repeatProgressChecks means the `main` process have completed or very close to completion
        #   Good to make sure though
        print("\n\nRunning progress checks...")
        time.sleep(0.5)
        subprocess.Popen([sys.executable,"repeatProgressChecks.py","-t","10","-n",str(_SET_nProcess),"-f",_SET_fileTag]).wait() 
        exit_codes = [proc.wait() for proc in openProcesses]

def mergeResults(_SET_fileTag):
    '''
    After running the multi process q-factors there will be a bunch of resultsX.root files
    We will use the program "hadd" to add all the resulting root files together and then merge the final result file with the input tree
    so we can have a tree that has everything in it
    '''
    print("\n\n")
    concatRootCmd="hadd -f logs"+_SET_runTag+"/"+_SET_fileTag+"/resultsMerged_"+_SET_fileTag+".root"
    for iProcess in range(_SET_nProcess):
        concatRootCmd=concatRootCmd+" logs"+_SET_runTag+"/"+_SET_fileTag+"/results"+str(iProcess)+".root"
    print(concatRootCmd)
    subprocess.Popen(concatRootCmd,shell=True).wait()
    # Finally if everything successfully exists we can merge out q-factor results with the original input tree
    print("All processes completed without error codes")
    print("Begin merging qfactor results...")
    os.system("root -l -b -q mergeQresults.C")
            
       
def runMakeGraphs(_SET_fileTag,_SET_emailWhenFinished):
    '''
    Clean and draw the histograms
    '''
    # clean up the files that are created by makeDiagnosticHists before we rerun it
    print("Cleaning diagnosticPlots"+_SET_runTag+" folder")
    os.system("rm -rf diagnosticPlots"+_SET_runTag+"/"+_SET_fileTag)
    os.system("mkdir -p diagnosticPlots"+_SET_runTag+"/"+_SET_fileTag)
    # ------------------------------------
    # run the makeDiagnosticHists program
    # ------------------------------------
    cmd='root -l -b -q "makePlots.C(false)"'
    print("running cmd: "+cmd)
    subprocess.Popen(cmd,shell=True).wait()

def combineAllGraphs(_SET_fileTag):
    '''
    After running over all the different dataset we will just add all the histograms together and make one final plot
    postQVal_hists_TAG.root are outputs of makePlots.C
    postQVal_flatTree_TAG.root are outputs of mergeQresults.C
    '''
    tags = [rootFileLoc[2] for rootFileLoc in rootFileLocs]
    haddHistCmd="hadd -f diagnosticPlots"+_SET_runTag+"/postQVal_hists.root"
    haddTreeCmd="hadd -f logs"+_SET_runTag+"/postQVal_flatTree.root"
    for tag in tags:
        haddHistCmd = haddHistCmd+" diagnosticPlots"+_SET_runTag+"/"+_SET_fileTag+"/postQVal_hists_"+_SET_fileTag+".root"
        haddTreeCmd = haddTreeCmd+" logs"+_SET_runTag+"/"+_SET_fileTag+"/postQVal_flatTree_"+_SET_fileTag+".root"
    print("\n\n")
    print(haddHistCmd)
    print(haddTreeCmd)
    os.system("rm -f diagnosticPlots/postQVal_hists.root")
    os.system(haddHistCmd)
    os.system("rm -f diagnosticPlots/postQVal_flatTree.root")
    os.system(haddTreeCmd)

    cmd='root -l -b -q "makePlots.C(true)"'
    print("running cmd: "+cmd)
    subprocess.Popen(cmd,shell=True).wait()


##############################################################################
####################    BEGIN RUNNING THE PROGRAM    #########################
##############################################################################
for _SET_rootFileLoc, _SET_rootTreeName, _SET_fileTag in rootFileLocs:
    print("\n\n-------------------------------------")
    print("Starting running of {0}".format(_SET_rootFileLoc))
    print("TreeName: {0}".format(_SET_rootTreeName))
    print("FileTag: {0}".format(_SET_fileTag))

    numVar=len(varVec)

    if not _SET_runAllPhaseCombos:
        combo=range(numVar)
        tagVec=["0" for i in range(len(varVec))]
        for ele in combo:
            tagVec[ele]="1"
        paddedCombo="".join(tagVec)
        reconfigureSettings("configSettings.h",combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag+"_"+paddedCombo)
        if _SET_runQFactor:
            runOverCombo(combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag+"_"+paddedCombo)
        if _SET_runMakeHists:
            mergeResults(_SET_fileTag+"_"+paddedCombo)
            runMakeGraphs(_SET_fileTag+"_"+paddedCombo,_SET_emailWhenFinished)
    #if _SET_runMakeHists:
    #    combineAllGraphs(_SET_fileTag+"_"+paddedCombo)
    else:
        for numVar in range(1,len((varVec))+1):
            combos=combinations(range(len(varVec)),numVar)
            for combo in combos:
                tagVec=["0" for i in range(len(varVec))]
                for ele in combo:
                    tagVec[ele]="1"
                paddedCombo="".join(tagVec)
                reconfigureSettings("configSettings.h",combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag+"_"+paddedCombo)
                if _SET_runQFactor:
                    runOverCombo(combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag+"_"+paddedCombo)
                if _SET_runMakeHists:
                    mergeResults(_SET_fileTag+"_"+paddedCombo)
                    runMakeGraphs(_SET_fileTag+"_"+paddedCombo,_SET_emailWhenFinished)
    if _SET_emailWhenFinished:
        print("Sending program finished email")
        subprocess.Popen("sendmail "+_SET_emailWhenFinished+" < defaultEmail.txt",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()

# Once all the datasets have been run over we can combine all the results. This assumes the datasets should be combined...
    
# OLD: Use this code block to run over all possible combinations of variables
# We are going pass as arugment a list of lists known as combo. This combo list contains all the lists of combos with numVar elements from the list varVec. If we use the command comboinations(range(3),2) we would get something like [ [1,2], [2,3], [1,3] ]. We can use these as indicies to index a a string of 0's to fill in whether a variable will be in use. i.e. if [1,3] is chosen then the string would be 101 with the second var turnedo off. This is useful when we are doing a scan of which variables we should use. Bruteforce style. 

print("--- %s seconds ---" % (time.time() - start_time))







