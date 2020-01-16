#!/usr/bin/env python

import os
import subprocess
import time
import argparse
import numpy as np
from multiprocessing import Process, Pool, cpu_count

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config

def parser():
    parser = argparse.ArgumentParser(description = "Skim MINIAOD with crab", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--submit", action = "store_true", help = "Submit all the jobs")
    parser.add_argument("--monitor", action = "store_true", help = "Check if jobs only should be monitored")
    parser.add_argument("--merge", action = "store_true", help = "Merge output together")
    parser.add_argument("--process", action = "store", default = "all", help = "Merge output together")

    return parser.parse_args()

def crabConfig(dataSet, setName, outDir, systematics):
    isSignal = "HPlus" in setName
    isData = "Single" in setName

    channels = ["mu4j", "e4j", "mu2j1fj", "e2j1fj", "mu2fj", "e2fj"]

    outFiles = []

    for systematic in systematics:
        if systematic == "":
            outFiles.append("{}.root".format(setName))
            continue

        if isData:
            break

        for shift in ["Up", "Down"]:
            outFiles.append("{}_{}{}.root".format(setName, systematic, shift))

    ##Crab config
    crabConf = config()

    crabConf.General.requestName = "Skim"
    crabConf.General.workArea = outDir
    crabConf.General.transferOutputs = True
    crabConf.General.transferLogs = False

    crabConf.JobType.allowUndistributedCMSSW = True
    crabConf.JobType.pluginName = "Analysis"
    crabConf.JobType.psetName = "ChargedSkimming/Skimming/python/miniskimmer.py"
    crabConf.JobType.pyCfgParams = ["outname={}.root".format(setName), "channel={}".format(",".join(channels))]
    crabConf.JobType.outputFiles = outFiles
    crabConf.JobType.maxMemoryMB = 3500

    crabConf.Data.inputDataset = dataSet
    crabConf.Data.inputDBS = "global" if not isSignal else "phys03"
    crabConf.Data.splitting = "EventAwareLumiBased" if not isSignal else "FileBased"
    crabConf.Data.unitsPerJob = 150000 if not isSignal else 2
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

                runTime = "1400" if not 50664 in exitCode else "1500"
                memory = "3500" if not 50660 in exitCode else "3600"

                crabCommand("resubmit", dir=crabDir, maxmemory=memory, maxjobruntime=runTime)
                break 

    if "finished" in crabStatus["jobsPerStatus"]:
        jobNrs = [nr for nr in crabStatus["jobs"].keys() if crabStatus["jobs"][nr]["State"] == "finished"]

        ##Skip if output already retrieved
        if(len(crabStatus["jobs"].keys()) == len(jobNrs)):
            return True

    time.sleep(10)

    return False

def merge(crabJob, systematics):
    ##Create output directories
    crabDir = "{}/crab_{}".format(crabJob.General.workArea, crabJob.General.requestName)
    subprocess.call(["mkdir", "-p", "{}/merged/".format(crabJob.General.workArea)])
    setName = crabJob.General.workArea.split("/")[-1]

    ##Get all output files
    allFiles = subprocess.check_output(["crab", "getoutput", "--xrootd", crabDir]).split("\n")[:-1]

    for systematic in systematics:
        for shift in ["Up", "Down"]:
            ##Get files for specific systematic + shift
            if systematic == "":
                files = [f for f in allFiles if not "Down" in f and not "Up" in f]
                if shift == "Down":
                    continue

            else:
                files = [f for f in allFiles if shift in f and systematic in f]

            systName = "" if systematic == "" else "_{}{}".format(systematic, shift) 

            ##Create directory with dag files
            subDir = "{}/Tmp/SkimMerge/{}{}".format(os.environ["CHDIR"], setName, systName)
            subprocess.call(["mkdir", "-p", subDir])

            ##Condor and dagman templates
            condorSub = [
                "universe = vanilla",
                "arguments = $(FILES)",
                "getenv = True",
                'requirements = (OpSysAndVer =?= "CentOS7")',
                "error = {}/err$(Process).txt".format(subDir),
                "output = {}/out$(Process).txt".format(subDir),
                "executable = {}/src/ChargedSkimming/Skimming/scripts/merge.sh".format(os.environ["CMSSW_BASE"]),
            ]

            dagmanSub = [
                "JOB A {}/subMerge.sub".format(subDir),
                "JOB B {}/merge.sub".format(subDir),
                "JOB C {}/rm.sub".format(subDir),
                "PARENT A CHILD B",
                "PARENT B CHILD C",
            ]

            ##Split files in chuncks
            if(len(files) == 0):
                return None
            files = np.array_split(files, 1 if len(files) < 40 else int(len(files)/40.))
            
            tmpOutput = []

            ##Write condor script for merge of chunks of files
            with open("{}/subMerge.sub".format(subDir), "w") as condFile:
                for line in condorSub:
                    condFile.write(line + "\n")
                    
                condFile.write("queue FILES from (\n")

                for index, filesBunch in enumerate(files):
                    tmpFile = "{}/merged/tmpFile{}_{}.root".format(crabJob.General.workArea, systName, index)
                    condFile.write(" ".join([tmpFile] + list(filesBunch)) + "\n")

                    tmpOutput.append(tmpFile)

                condFile.write(")")

            ##Write condor script for merge of chunks to complete file
            with open("{}/merge.sub".format(subDir), "w") as condFile:
                for line in condorSub:
                    condFile.write(line + "\n")

                mergeOutput =  "{}/merged/{}{}.root".format(crabJob.General.workArea, setName, systName)

                condFile.write("queue FILES from (\n")
                condFile.write(" ".join([mergeOutput] + tmpOutput) + "\n")
                condFile.write(")")

            ##Write condor script to clean up merge directory
            with open("{}/rm.sub".format(subDir), "w") as condFile:
                for line in condorSub[:-1]:
                    condFile.write(line + "\n")

                condFile.write("executable = /bin/rm\n")

                condFile.write("queue FILES from (\n")
                condFile.write(" ".join(tmpOutput) + "\n")
                condFile.write(")")

            ##Write and call dagman script
            with open("{}/merge.dag".format(subDir), "w") as dagFile:
                for line in dagmanSub:
                    dagFile.write(line + "\n")
            
            subprocess.call(["condor_submit_dag", "-f", "{}/merge.dag".format(subDir)])

def main():
    ##Parser arguments
    args = parser()

    ##Txt with dataset names
    filePath = "{}/src/ChargedSkimming/Skimming/data/filelists".format(os.environ["CMSSW_BASE"])

    if args.process == "all":
        fileLists = [
                    filePath + "/filelist_bkg_2017_MINI.txt",
                    filePath + "/filelist_data_2017_MINI.txt",
                    filePath + "/filelist_signal_2017_MINI.txt",
        ]

    elif args.process == "data":
        fileLists = [
                    filePath + "/filelist_data_2017_MINI.txt",
        ]

    elif args.process == "sig":
        fileLists = [
                    filePath + "/filelist_signal_2017_MINI.txt",
        ]

    elif args.process == "bkg":
        fileLists = [
                    filePath + "/filelist_bkg_2017_MINI.txt",
        ]

    ##Create with each dataset a crab config
    crabJobs = []

    systematics = ["", "energyScale", "energySigma", "JECTotal", "JER"]
    
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

                crabJobs.append(crabConfig(dataset, name, "{}/Skim/{}".format(os.environ["CHDIR"], name), systematics))

    ##Submit all crab jobs
    if args.submit:
        processes = [Process(target=submit, args=(config,)) for config in crabJobs]

        for process in processes:
            process.start()
            process.join()

    ##Merge output files
    elif args.merge:
        pool = Pool(processes=cpu_count())
        mergeJobs = []

        for crabJob in crabJobs:
            mergeJobs.append(pool.apply_async(merge, (crabJob, systematics)))

        [job.get() for job in mergeJobs]

    ##Just monitor crab jobs
    elif args.monitor:
        while(True):
            for crabJob in crabJobs:
                isFinished = monitor(crabJob)
                
                if isFinished:
                    crabJobs.remove(crabJob)    
                    continue
    
if __name__ == "__main__":
    main()
