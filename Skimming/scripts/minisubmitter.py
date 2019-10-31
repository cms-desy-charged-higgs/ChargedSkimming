#!/usr/bin/env python

import os
import subprocess
import time
import argparse
from multiprocessing import Process, Pool, cpu_count

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config

def parser():
    parser = argparse.ArgumentParser(description = "Skim MINIAOD with crab", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--submit", action = "store_true", help = "Submit all the jobs")
    parser.add_argument("--monitor", action = "store_true", help = "Check if jobs only should be monitored")
    parser.add_argument("--merge", action = "store_true", help = "Merge output together")

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
    crabConf.JobType.maxJobRuntimeMin = 1750

    crabConf.Data.inputDataset = dataSet
    crabConf.Data.inputDBS = "global" if not isSignal else "phys03"
    crabConf.Data.splitting = "EventAwareLumiBased" if not isSignal else "FileBased"
    crabConf.Data.unitsPerJob = 250000 if not isSignal else 1
    crabConf.Data.outLFNDirBase = "/store/user/dbrunner/skim"

    crabConf.Site.storageSite = "T2_DE_DESY"
    crabConf.User.voGroup = "dcms"

    return crabConf

def submit(crabJob):
    crabCommand("submit", config=crabJob)

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

                runTime = "1750" if not 50664 in exitCode else "1900"
                memory = "3000" if not 50660 in exitCode else "3200"

                crabCommand("resubmit", dir=crabDir, maxmemory=memory, maxjobruntime=runTime)
                break 

    if "finished" in crabStatus["jobsPerStatus"]:
        jobNrs = [nr for nr in crabStatus["jobs"].keys() if crabStatus["jobs"][nr]["State"] == "finished"]

        ##Skip if output already retrieved
        if(len(os.listdir("{}/results/".format(crabDir))) == len(jobNrs)):
            if len(os.listdir("{}/results/".format(crabDir))) == len(crabStatus["jobs"].keys()):
                return True

            return False

        ##Retrieve output
        for jobList in [jobNrs[i:i + 500] for i in range(0, len(jobNrs), 500)]:
            while(True):
                try:
                    crabCommand("getoutput", dir=crabDir, jobids=",".join(jobList))
                    break
                except:
                    pass

    return False

def merge(crabJob):
    resultDir = "{}/crab_{}/results".format(crabJob.General.workArea, crabJob.General.requestName)
    files = ["{}/{}".format(resultDir, f) for f in os.listdir(resultDir)]
    
    mergeOutput =  "{}/merged/{}.root".format(crabJob.General.workArea, "_".join(crabJob.General.requestName.split("_")[2:]))

    if(len(files) == 0):
        return None

    subprocess.call(["mkdir", "-p", "{}/merged/".format(crabJob.General.workArea)])
    subprocess.call(["hadd", "-n", "50", "-f", mergeOutput] + files)       

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
    if args.submit:
        processes = [Process(target=submit, args=(config,)) for config in crabJobs]

        for process in processes:
            process.start()
            process.join()

    ##Merge output files
    elif args.merge:
        for crabJob in crabJobs:
            merge(crabJob)

    ##Just monitor crab jobs
    elif args.monitor:
        while(True):
            for crabJob in crabJobs:
                isFinished = monitor(crabJob)
                
                if isFinished:
                    crabJobs.remove(crabJob)    
                    continue
                
                time.sleep(30)
    
if __name__ == "__main__":
    main()
