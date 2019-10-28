#!/usr/bin/env python

import os
import subprocess
import time
import argparse
from multiprocessing import Process

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config

def parser():
    parser = argparse.ArgumentParser(description = "Skim MINIAOD with crab", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--monitor", action = "store_true", help = "Check if jobs only should be monitored")

    return parser.parse_args()

def crabConfig(dataSet, setName, outDir):
    isSignal = "HPlus" in setName

    ##Crab config
    crabConf = config()

    crabConf.General.requestName = "MiniSkim_{}".format(setName)
    crabConf.General.workArea = outDir
    crabConf.General.transferOutputs = True
    crabConf.General.transferLogs = False

    crabConf.JobType.pluginName = "Analysis"
    crabConf.JobType.psetName = "ChargedSkimming/Skimming/python/miniskimmer.py"
    crabConf.JobType.pyCfgParams = ["outname={}.root".format(setName)]
    crabConf.JobType.outputFiles = ["{}.root".format(setName)]
    crabConf.JobType.maxMemoryMB = 3000
    crabConf.JobType.maxJobRuntimeMin = 1440

    crabConf.Data.inputDataset = dataSet
    crabConf.Data.inputDBS = "global" if not isSignal else "phys03"
    crabConf.Data.splitting = "EventAwareLumiBased" if not isSignal else "FileBased"
    crabConf.Data.unitsPerJob = 250000 if not isSignal else 1
    crabConf.Data.outLFNDirBase = "/store/user/dbrunner/skim"

    crabConf.Site.storageSite = "T2_DE_DESY"
    crabConf.User.voGroup = "dcms"

    return crabConf

def monitor(crabJob):
    ##Crab dir
    crabDir = "{}/crab_{}".format(crabJob.General.workArea, crabJob.General.requestName)

    ##Get status
    crabStatus = crabCommand("status", dir=crabDir)
    nFinished = 0
    
    ##If submission failed, resubmit and leave function
    if crabStatus["status"] == "SUBMITFAILED":
        os.system("command rm -r {}".format(crabDir))

        Process(target=crabCommand, args=("submit", crabJob))

        return False

    for status in crabStatus["jobsPerStatus"].keys():
        ##If one jobs failed, resubmit and leave loop
        if status == "failed" and (crabStatus["dbStatus"] == "SUBMITTED" or crabStatus["dbStatus"] == "RESUBMITFAILED"):
            if crabStatus["jobsPerStatus"][status] != 0:
                ##Check if you have to increase memory/run time
                exitCode = [crabStatus["jobs"][key]["Error"][0] for key in crabStatus["jobs"] if "Error" in crabStatus["jobs"][key]]

                runTime = "1440" if not 50664 in exitCode else "1700"
                memory = "3000" if not 50660 in exitCode else "3200"

                crabCommand("resubmit", dir=crabDir, maxmemory=memory, maxjobruntime=runTime)
                break

        if "finished" in crabStatus["jobsPerStatus"]:
            if crabStatus["jobsPerStatus"]["finished"] == len(crabStatus["jobs"].keys()):
                jobsNrs = crabStatus["jobs"].keys()

                ##Skip if output already retrieved
                if(len(os.listdir("{}/results/".format(crabDir))) == len(jobsNrs)):
                    return True

                ##Retrieve output
                for jobList in [jobsNrs[i:i + 500] for i in range(0, len(jobsNrs), 500)]:
                    while(True):
                        try:
                            crabCommand("getoutput", dir=crabDir, jobids=",".join(jobList))
                            break
                        except:
                            pass       
     
                ##If all jobs finished, retrieve output
                return True

    return False
        
def main():
    ##Parser arguments
    args = parser()

    ##Txt with dataset names
    filePath = "{}/src/ChargedSkimming/Skimming/data/filelists".format(os.environ["CMSSW_BASE"])

    fileLists = [
                filePath + "/filelist_bkg_2017_MINI.txt",
                filePath + "/filelist_data_2017_MINI.txt",
                filePath + "/filelist_signal_2017_MINI.txt",
    ]

    ##Create with each dataset a crab config
    crabJobs = []
    
    for fileName in fileLists:
        with open(fileName) as f:
            datasets = [dataset for dataset in f.read().splitlines() if dataset != ""]

            for dataset in datasets:
                if "SIM" in dataset or "HPlus" in dataset:
                    if "ext" in dataset:
                        name = dataset.split("/")[1] + "_ext"
                    else: 
                        name = dataset.split("/")[1]
                
                else: 
                    name = dataset.split("/")[1] + "_" + dataset.split("/")[2]


                crabJobs.append(crabConfig(dataset, "MiniSkim_{}".format(name), "{}/Skim/{}".format(os.environ["CHDIR"], name)))

    ##Submit all crab jobs
    if not args.monitor:
        processes = [Process(target=crabCommand, args=("submit", config)) for config in crabJobs]

        for process in processes:
            process.start()
            process.join()

    ##Just monitor crab jobs
    else:
        while(True):
            for crabJob in crabJobs:
                isFinished = monitor(crabJob)
                
                if isFinished:
                    crabJobs.remove(crabJob)    
                    continue
                
                time.sleep(30)
    
if __name__ == "__main__":
    main()
