import subprocess
import time
from termcolor import colored
import argparse
from os import path

def checkProgress(nProcess, tag, percentages):
    folder="logs/"+tag
    nentriesPerProc=0
    proc_status=[False]*nProcess
    current_percentages=[0 for i in range(nProcess)] # Useful to compare previous percentages with current to see if progress appear stalled
    for iproc in range(nProcess):
        if not proc_status[iproc]:
            #######################
            # Load the last n lines of process log
            #######################
            p=subprocess.Popen(["tail", "-n20",folder+"/"+"processLog"+str(iproc)+".txt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output=p.communicate()[0].decode("utf-8")
            if path.exists(folder+"/"+"processLog"+str(iproc)+".txt") and output!="": # make sure we have some lines the log file
                allLines=output.split("\n") 
                newestLine=""
                allLines.reverse() # reverse to get the newest line
                for line in allLines:
                    if line[:14]=="Starting event": # Check the progress
                        newestLine=line
                    if line[:10]=="Total time": # check if completed
                        proc_status[iproc]=True 

                #######################
                # calculate percentages
                #######################
                current, end = newestLine.split(" ")[2].rstrip().lstrip().split("/")
                current, end = int(current), int(end)
                if iproc==0:
                    nentriesPerProc=end
                    perc_complete = int(100.0*current/end)
                else:
                    perc_complete = int(100.0*(current-(end-nentriesPerProc))/nentriesPerProc)
                current_percentages[iproc]=perc_complete

                #######################
                # Check if complete, if not show status
                #######################
                if proc_status[iproc]: # Process completed this iteration?
                    print(colored("process {:d} is complete".format(iproc),"green"))
                else: # show progress or stalled
                    if current_percentages[iproc]==percentages[iproc]:
                        print(colored("process {:d} did not show progress since last check: ({:d}/{:d})  ---- {:d}%".format(iproc,current,end,perc_complete),"red"))
                    else:
                        print(colored("process {:d} percent completed: ({:d}/{:d})  ---- {:d}%".format(iproc,current,end,perc_complete),"yellow"))
            else: # if log file does not exist or is empty then it is still initializing
                print(colored("process {:d} was still initializing, lets give it some more time...".format(iproc),"red"))
                continue
        else: # this checks if the process was already complete before the current check
            print(colored("process {:d} is complete".format(iproc),"green"))
    print("----------------------")
    return proc_status, current_percentages

if __name__ == "__main__":
    parser=argparse.ArgumentParser("Checking completion of main programs")
    parser.add_argument("-t",help="time in seconds to wait between progress checks",default=10, type=int)
    parser.add_argument("-n",help="nProcess", type=int)
    parser.add_argument("-f",help="tag of the folder to look for the log information")
    args=parser.parse_args()
    time_step=args.t
    nProcess=args.n
    tag=args.f
    assert(nProcess!=None and tag!=None)
        
    proc_status, percentages=checkProgress(nProcess, tag, [0 for i in range(nProcess)])
    while not all(proc_status):
        time.sleep(time_step)
        proc_status, percentages=checkProgress(nProcess, tag, percentages)







