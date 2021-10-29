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
from termcolor import colored

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from CRABClient.ClientExceptions import TaskNotFoundException, CachefileNotFoundException, ConfigException
from dbs.apis.dbsClient import DbsApi

def parser():
    parser = argparse.ArgumentParser(description = "Skim MINIAOD with crab", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("skim_dir", metavar='skim-dir', type=str, action="store", help="Working directory")
    parser.add_argument("--eras", nargs='+', default = ["2016", "2017", "2018"], help = "Run eras")
    parser.add_argument("--channels", nargs='+', default = ["MuonIncl", "EleIncl"], help = "Merge output together")
    parser.add_argument("--submit", action = "store_true", help = "Submit all the jobs")
    parser.add_argument("--monitor", action = "store_true", help = "Check if jobs only should be monitored")
    parser.add_argument("--get-output", action = "store_true", help = "Save list of all output files")

    return parser.parse_args()

def crabConfig(dataSet, setName, outDir, systematics, channels, era):
    isSignal = "HPlus" in setName
    isData = "Single" in setName or "JetHT" in setName or "EGamma" in setName

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
    url="https://cmsweb.cern.ch/dbs/prod/{}/DBSReader".format("global") # if not isSignal else "phys03")
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
    crabConf.JobType.psetName = "{}/src/ChargedSkimming/Skimming/python/miniskimmer.py".format(os.environ["CMSSW_BASE"])
    crabConf.JobType.pyCfgParams = ["outname={}.root".format(setName), "channel={}".format(",".join(channels)), "era={}".format(era)]
    crabConf.JobType.outputFiles = outFiles
    crabConf.JobType.maxJobRuntimeMin = 1440
    crabConf.JobType.maxMemoryMB = 2500
    crabConf.JobType.allowUndistributedCMSSW = True

    crabConf.Data.inputDataset = dataSet
    crabConf.Data.inputDBS = "global" # if not isSignal else "phys03"
    crabConf.Data.splitting = "FileBased"
    crabConf.Data.unitsPerJob = filesPerJob
    crabConf.Data.outLFNDirBase = "/store/user/dbrunner/skim/{}/{}".format("_".join([str(getattr(time.localtime(), "tm_" + t)) for t in ["mday", "mon", "year"]]), era)

    crabConf.Site.storageSite = "T2_DE_DESY"
    crabConf.User.voGroup = "dcms"

    return crabConf

def submit(era, name, config):
    crabDir = "{}/crab_{}".format(config.General.workArea, config.General.requestName)

    try:
        sys.stdout = open(os.devnull, 'w')
        crabCommand("submit", config=config)
        sys.stdout = sys.__stdout__

        print("{} ({}): Succesfully submitted".format(name, era))

    except ConfigException as e:
        try:
            crabCommand("status", dir=crabDir)

        except CachefileNotFoundException as e2:
            subprocess.check_output("command rm -rfv {}".format(crabDir), shell = True)
            crabCommand("submit", config=config)

            sys.stdout = sys.__stdout__

            print("{} ({}): Succesfully submitted".format(name, era))
            return
      
        if "already exists" in str(e):
            sys.stdout = sys.__stdout__
            print("{} ({}): Already submitted".format(name, era))

def monitor(era, name, config):
    out = "{} ({}): ".format(name, era)

    ##Crab dir
    crabDir = "{}/crab_{}".format(config.General.workArea, config.General.requestName)

    ##Get status
    try:
        sys.stdout = open(os.devnull, 'w')
        crabStatus = crabCommand("status", dir=crabDir)
        sys.stdout = sys.__stdout__

    except:
        sys.stdout = sys.__stdout__
        submit(era, name, config)

        return ""

    nFailed, nFinished, nRunning, nIdle = 0, 0, 0, 0
    
    ##If submission failed, resubmit and leave function
    if crabStatus["status"] == "SUBMITFAILED":
        os.system("command rm -rf {}".format(crabDir))
        submit(era, name, config)

        return ""

    for idx, job in crabStatus["jobs"].items():
        ##If one jobs failed, resubmit and leave loop
        if job["State"] == "failed":
            nFailed += 1

        elif job["State"] == "running":
            nRunning += 1

        elif job["State"] == "finished":
            nFinished += 1

        else:
            nIdle += 1

    out += "{}/{}/{}/{}".format(nIdle, colored(nRunning, "blue"), colored(nFinished, "green"), colored(nFailed, "red"))

    if nFailed != 0 and (crabStatus["dbStatus"] == "SUBMITTED" or crabStatus["dbStatus"] == "RESUBMITFAILED"):
        try:
            sys.stdout = open(os.devnull, 'w')
            crabCommand("resubmit", dir=crabDir)
            sys.stdout = sys.__stdout__

        except:
            pass

        out += " (Resubmitting...)"

    if nFinished == len(crabStatus["jobs"]):
        out += " (Finished)"

    return out

def getOutput(era, name, config):
    out = "{} ({}): ".format(name, era)
        
    crabDir = "{}/crab_{}".format(config.General.workArea, config.General.requestName)

    sys.stdout = open(os.devnull, 'w')
    crabStatus = crabCommand("status", dir=crabDir)
    sys.stdout = sys.__stdout__

    jobNrs = [nr for nr in crabStatus["jobs"].keys() if crabStatus["jobs"][nr]["State"] == "finished"]
    jobNrs.sort()

    while True:
        files = []

        try:
            for ids in np.array_split(jobNrs, np.arange(100, len(jobNrs), 100)):
                sys.stdout = open(os.devnull, 'w')
                files.extend(crabCommand("getoutput", dir = crabDir, dump = True, jobids= ",".join(ids))["lfn"])
                sys.stdout = sys.__stdout__

            break

        except:
            pass

    if config.Site.storageSite == "T2_DE_DESY":
        files = ["/pnfs/desy.de/cms/tier2/" + f for f in files]

    else:
        files = ["root://cms-xrd-global.cern.ch" + f for f in files]

    if "Run2018D" in name:
        splittedFile = np.array_split(files, 3)
        
        files = splittedFile[0]

        for index, fList in enumerate(splittedFile[1:]):
            newDir = config.General.workArea.replace("Run2018D", "Run2018D{}".format(index+2))
            subprocess.call(["mkdir", "-pv", newDir])

            with open("{}/outputFiles.txt".format(newDir), "w") as fileList:
                for f in fList:
                    fileList.write(f)
                    fileList.write("\n")

    with open("{}/outputFiles.txt".format(config.General.workArea), "w") as fileList:
        for f in files:
            fileList.write(f)
            fileList.write("\n")

    nSyst = len(files)/len(crabStatus["jobs"].keys())

    print("{} ({}): {} of {} extracted".format(name, era, len(files), len(crabStatus["jobs"].keys())*nSyst))
 
def main():
    ##Parser arguments
    args = parser()

    ##Txt with dataset names
    processes = ["sig"]
    systematics = ["", "energyScale", "energySigma", "JECTotal", "JER"]

    ##Create with each dataset a crab config
    crabConfigs = {}

    for era in args.eras:
        crabConfigs[era] = {}

        for process in processes:
            filePath = "{}/src/ChargedSkimming/Skimming/data/filelists/{}/filelist_{}_MINI.yaml".format(os.environ["CMSSW_BASE"], era, process)

            with open(filePath) as f:
                datasets = yaml.load(f, Loader=yaml.Loader)

                for dataset in datasets:
                    if "SIM" in dataset or "HPlus" in dataset:
                        name = dataset.split("/")[1]

                        nExt = len([n for n in crabConfigs[era].keys() if name in n])

                        if nExt != 0:
                            name += "_ext{}".format(nExt)
                
                    else: 
                        name = dataset.split("/")[1] + "_" + dataset.split("/")[2]

                    crabConfigs[era][name] = crabConfig(dataset, name, "{}/{}/{}/".format(args.skim_dir, era, name), systematics, args.channels, era)

    ##Submit all crab jobs
    if args.submit:
        print("\nStart submission of crab jobs\n")

        for era in args.eras:
            jobs = [Process(target=submit, args = (era, name, config)) for name, config in crabConfigs[era].items()]

            for job in jobs:
                job.start()
                job.join()

    if args.monitor:
        p = Pool(cpu_count())

        while crabConfigs:
            print("\nStatus of crab jobs ({})\n".format(time.asctime()))

            for era in args.eras:
                jobs = {name: p.apply_async(monitor, (era, name, config)) for name, config in crabConfigs[era].items()}

                for name, job in jobs.items():
                    out = job.get()
                    if out != "":
                        print(out)

                    if "Finished" in out:
                        crabConfigs[era].pop(name)

                time.sleep(60)

                if not crabConfigs[era]:
                    crabConfigs.pop(era)

    if args.get_output:
        p = Pool(cpu_count())

        for era in args.eras:
            jobs = {name: p.apply_async(getOutput, (era, name, config)) for name, config in crabConfigs[era].items()}
            [job.get() for job in jobs.values()]

if __name__ == "__main__":
    main()
