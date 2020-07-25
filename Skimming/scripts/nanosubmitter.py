#!/usr/bin/env python

import os
import subprocess
import argparse
import yaml

from dbs.apis.dbsClient import DbsApi

def parser():
    parser = argparse.ArgumentParser(description = "Skim MINIAOD with crab", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--submit", action = "store_true", help = "Submit all the jobs")
    parser.add_argument("--channel", nargs='+', default = ["mu4j", "e4j"], help = "Channel to be analyzed")
    parser.add_argument("--merge", action = "store_true", help = "Merge output together")
    parser.add_argument("--process", action = "store", default = "all", help = "Type of process to be analyzed")
    parser.add_argument("--skim-dir", action = "store", default = "Skim", help = "Type of process to be analyzed")

    return parser.parse_args()


def submit(jobInfo, channel, dirName):
    skimDir = "{}/{}".format(os.environ["CHDIR"], dirName)

    for d in ["out", "log", "err"]:
        subprocess.call(["mkdir", "-p", "{}/{}".format(skimDir, d)])

    condorSub = [
            "universe = vanilla",
            "executable = {}/src/ChargedSkimming/Skimming/batch/produceSkim.sh".format(os.environ["CMSSW_BASE"]),
            "arguments = $(fileName) $(isData) '{}' $(xSec) $(outName) $(job) $(skimDir)".format(" ".join(channel)),
            "log = {}/log/$(outName)_$(job).log".format(skimDir),
            "error = {}/err/$(outName)_$(job).err".format(skimDir),
            "getenv = True",
            "output = {}/out/$(outName)_$(job).out".format(skimDir),
    ]

    with open("{}/condor.sub".format(skimDir), "w") as condFile:
        for line in condorSub:
            condFile.write(line + "\n")
            
        condFile.write("queue fileName isData xSec outName job skimDir from (\n")
        
        for job in jobInfo:
            condFile.write("{} {} {} {} {} {}\n".format(job[0], str(job[1]), str(job[2]), job[3], str(job[4]), skimDir))

        condFile.write(")")

    subprocess.call(["condor_submit", "{}/condor.sub".format(skimDir)])


def merge(dirName):
    skimDir = "{}/{}".format(os.environ["CHDIR"], dirName)
    procDir = [d for d in os.listdir(skimDir) if d not in ['condor.sub', 'log', 'err', 'out']]

    for dir in procDir:
        subprocess.call(["mkdir", "-p", "{}/{}/merged".format(skimDir, dir)])
        subprocess.call(["hadd", "-f", "{}/{}/merged/{}.root".format(skimDir, dir, dir)] + ["{}/{}/output/{}".format(skimDir, dir, fileName) for fileName in os.listdir("{}/{}/output".format(skimDir, dir))])

def main():
    ##Parser arguments
    args = parser()

    ##Txt with dataset names
    filePath = "{}/src/ChargedSkimming/Skimming/data/filelists".format(os.environ["CMSSW_BASE"])

    if args.process == "all":
        fileLists = [
                    filePath + "/filelist_bkg_2017_NANO.txt",
                    filePath + "/filelist_data_2017_NANO.txt",
        ]

    elif args.process == "data":
        fileLists = [
                    filePath + "/filelist_data_2017_NANO.txt",
        ]

    elif args.process == "bkg":
        fileLists = [
                    filePath + "/filelist_bkg_2017_NANO.txt",
        ]

    if args.submit:
        ##Get xSec
        xSecFile = yaml.load(file("{}/src/ChargedSkimming/Skimming/data/xsec.yaml".format(os.environ["CMSSW_BASE"]), "r"))

        ##Create with each dataset a crab config
        dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
        jobInfo = []

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

                    xSec = 1.

                    for key in xSecFile.keys():
                        if key in name:
                            xSec = xSecFile[key]

                    ##Check if file is true data file
                    isData = True in [n in name for n in ["Electron", "Muon", "MET"]]

                    fileNames = ["root://cms-xrd-global.cern.ch/{}".format(fileName["logical_file_name"]) for fileName in dbs.listFiles(dataset=dataset)]

                    for i, fileName in enumerate(fileNames):
                        jobInfo.append((fileName, isData, xSec, name, i))

        submit(jobInfo, args.channel, args.skim_dir)

    if args.merge:
        merge(args.skim_dir)

if __name__ == "__main__":
    main()
