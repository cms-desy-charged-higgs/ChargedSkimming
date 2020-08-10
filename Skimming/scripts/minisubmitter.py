#!/usr/bin/env python

import os
import sys
import subprocess
import time
import argparse
import numpy as np
import yaml
import math
from multiprocessing import Process, Pool, cpu_count

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from dbs.apis.dbsClient import DbsApi

def parser():
    parser = argparse.ArgumentParser(description = "Skim MINIAOD with crab", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("skim_dir", metavar='skim-dir', type=str, action="store", help="Working directory")
    parser.add_argument("--era", action = "store", help = "Run era")
    parser.add_argument("--submit", action = "store_true", help = "Submit all the jobs")
    parser.add_argument("--monitor", action = "store_true", help = "Check if jobs only should be monitored")
    parser.add_argument("--get-output", action = "store_true", help = "Save list of all output files")
    parser.add_argument("--process", action = "store", default = "all", help = "Merge output together")
    parser.add_argument("--channels", nargs='+', default = ["MuonIncl", "EleIncl"], help = "Merge output together")

    return parser.parse_args()

def crabConfig(dataSet, setName, outDir, systematics, channels, era):
    isSignal = "HPlus" in setName
    isData = "Single" in setName or "JetHT" in setName

    outFiles = []

    for systematic in systematics:
        if systematic == "":
            outFiles.append("{}.root".format(setName))
            continue

        if isData:
            break

        for shift in ["Up", "Down"]:
            outFiles.append("{}_{}{}.root".format(setName, systematic, shift))

    #Caculate number of files per job
    url="https://cmsweb.cern.ch/dbs/prod/{}/DBSReader".format("global" if not isSignal else "phys03")
    api=DbsApi(url=url)
    files = api.listFiles(dataset = dataSet, detail = 1)
    
    eventsPerFile = sum(f["event_count"] for f in files)/len(files)
    filesPerJob = int(math.ceil(300000./eventsPerFile))

    ##Crab config
    crabConf = config()

    crabConf.General.requestName = "Skim_{}".format(era)
    crabConf.General.workArea = outDir
    crabConf.General.transferOutputs = True
    crabConf.General.transferLogs = False

    crabConf.JobType.pluginName = "Analysis"
    crabConf.JobType.psetName = "ChargedSkimming/Skimming/python/miniskimmer.py"
    crabConf.JobType.pyCfgParams = ["outname={}.root".format(setName), "channel={}".format(",".join(channels)), "era={}".format(era)]
    crabConf.JobType.outputFiles = outFiles
    crabConf.JobType.maxJobRuntimeMin = 1440
    crabConf.JobType.maxMemoryMB = 2500
    crabConf.JobType.allowUndistributedCMSSW = True

    crabConf.Data.inputDataset = dataSet
    crabConf.Data.inputDBS = "global" if not isSignal else "phys03"
    crabConf.Data.splitting = "FileBased"
    crabConf.Data.unitsPerJob = filesPerJob
    crabConf.Data.outLFNDirBase = "/store/user/dbrunner/skim/{}".format(era)

    crabConf.Site.storageSite = "T2_DE_DESY"
    crabConf.User.voGroup = "dcms"

    return crabConf

def submit(crabJob):
    crabCommand("submit", config=crabJob)

def getFiles(jobNr, crabDir):
    fileNames =  subprocess.check_output(["crab", "getoutput", "--xrootd", "--jobids", ",".join(jobNr), crabDir]).split("\n")[:-1]
    print("Got {} xrood paths".format(len(jobNr)))

    return fileNames
    
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

                runTime = "1440" if not 50664 in exitCode else "1500"
                memory = "2500" if not 50660 in exitCode else "3000"

                crabCommand("resubmit", dir=crabDir, maxmemory=memory, maxjobruntime=runTime)
                break 

    if "finished" in crabStatus["jobsPerStatus"]:
        jobNrs = [nr for nr in crabStatus["jobs"].keys() if crabStatus["jobs"][nr]["State"] == "finished"]

        ##Skip if output already retrieved
        if(len(crabStatus["jobs"].keys()) == len(jobNrs)):
            return True

    time.sleep(10)

    return False

def getOutput(crabJob):
    pool = Pool(processes=cpu_count())

    files = []
    name = crabJob.General.workArea.split("/")[-1]

    print("Get xrootd paths: {}".format(name))
        
    crabDir = "{}/crab_{}".format(crabJob.General.workArea, crabJob.General.requestName)

    sys.stdout = open(os.devnull, 'w')
    crabStatus = crabCommand("status", dir=crabDir)
    sys.stdout = sys.__stdout__

    jobNrs = [nr for nr in crabStatus["jobs"].keys() if crabStatus["jobs"][nr]["State"] == "finished"]
    jobNrs = [jobNrs[i:i + 30] for i in range(0, len(jobNrs), 30)]

    jobs = [pool.apply_async(getFiles, (jobNr, crabDir)) for jobNr in jobNrs]

    for job in jobs:
        files.extend(job.get())

    files.sort()

    print("Got all xrootd paths: {}".format(crabJob.General.workArea.split("/")[-1]))

    with open("{}/outputFiles.txt".format(crabJob.General.workArea), "w") as fileList:
        for f in files:
            fileList.write(f)
            fileList.write("\n")

def main():
    ##Parser arguments
    args = parser()

    ##Txt with dataset names
    filePath = "{}/src/ChargedSkimming/Skimming/data/filelists/{}".format(os.environ["CMSSW_BASE"], args.era)
    fileName = "/filelist_{}_MINI.yaml"
    
    procType = ["bkg", "sig", "data"] if args.process == "all" else [args.process]
    systematics = ["", "energyScale", "energySigma", "JECTotal", "JER"]

    ##Create with each dataset a crab config
    crabJobs = {}
   
    for process in procType:
        with open(filePath + fileName.format(process)) as f:
            datasets = yaml.load(f, Loader=yaml.Loader)

            for dataset in datasets:
                if "SIM" in dataset or "HPlus" in dataset:
                    if "ext" in dataset:
                        name = dataset.split("/")[1] + "_ext"
                    else: 
                        name = dataset.split("/")[1]
                
                else: 
                    name = dataset.split("/")[1] + "_" + dataset.split("/")[2]

                if args.era not in ["2016", "2017", "2018"]:
                    raise RuntimeError("Unvalid value for era: {}".format(args.era))

                crabJobs.setdefault(process, []).append(crabConfig(dataset, name, "{}/{}".format(args.skim_dir, name), systematics, args.channels, args.era))

    ##Submit all crab jobs
    if args.submit:
        for proc in procType:
            processes = [Process(target=submit, args=(config,)) for config in crabJobs[proc]]

            for process in processes:
                process.start()
                process.join()

    ##Get output files xrood paths and save in txt files
    if args.get_output:
        for proc in procType:
            for crabJob in crabJobs[proc]:
                getOutput(crabJob)

    ##Just monitor crab jobs
    elif args.monitor:
        while(True):
            for proc in procType:
                for crabJob in crabJobs[proc]:
                    try:
                        isFinished = monitor(crabJob)

                    except Exception as e:
                        print(str(e))
                        isFinished = False
                    
                    if isFinished:
                        crabJobs[proc].remove(crabJob)    
                        continue
    
if __name__ == "__main__":
    main()
