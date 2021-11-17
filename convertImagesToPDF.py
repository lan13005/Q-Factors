#!/usr/bin/env python

######################
# Useful to see results when running over all potential combinations of phase space variables, run.py's runAllPhaseCombos variable.  
#   Useful to grab a specific plot determined from imageFile in diagnosticPlots baseFolder and aggregates all occurences over all
#   potential combos.
######################

import os

baseFolder="diagnosticPlots"
imageFile="Mpi0g3.png"
annotate=True
gather=True
makeBackup=True

subfolders=os.listdir(baseFolder)
fileTags=[folder.split("_")[1] for folder in subfolders]
#print(fileTags)

searchStr="_SET_varStringBase="
with open("run.py","r") as cfgFile:
    cfgFileLines=cfgFile.readlines()
    cfgFileLines=[line for line in cfgFileLines if line.startswith(searchStr)] # search for correct line
    assert(len(cfgFileLines)==1) # there shoudlnt be more than one version of this line
    cfgFileLines="".join(cfgFileLines[0].split()) # remove all whitespaces, in case you are Jason and include extra white space
    cfgFileLines=cfgFileLines.rstrip().lstrip() # remove whitespaces
    cfgFileLines=cfgFileLines.split("#")[0] # remove comments
    cfgFileLines=cfgFileLines.split(searchStr)[1] # remove searchStr
    cfgFileLines=cfgFileLines.strip('"')
    variables=cfgFileLines.split(";")

#print(variables)
if makeBackup:
    print("Making a backup just in case...")
    os.system("rsync -r diagnosticPlots diagnosticPlots_backup")

if annotate:
    for subfolder in subfolders:
        print("annotating "+imageFile+" in "+subfolder)
        fileTag=subfolder.split("_")[1]
        variable_combo=[]
        for idx,ch in enumerate(fileTag):
            if ch=="1":
                variable_combo.append(variables[idx])
    
        variable_combo="\\n".join(variable_combo)
        variable_combo="'"+variable_combo+"'"
        
        mogrify_baseCmd="mogrify -font Liberation-Sans -fill white -undercolor '#00000080' -pointsize 26 -gravity NorthEast -annotate +20+500"
        command_args=[mogrify_baseCmd,variable_combo,baseFolder+"/"+subfolder+"/"+imageFile]
        #print(" ".join(command_args))
        os.system(" ".join(command_args))

if gather:
    print("Gathering images")
    imgPrefix=imageFile.split(".")[0]
    os.system("convert "+baseFolder+"/*/"+imageFile+" "+imgPrefix+"-gathered.pdf")
    
