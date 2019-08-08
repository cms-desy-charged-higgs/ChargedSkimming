from ChargedHiggs.Workflow.task import Task

import os
import sys
import time

sys.path.append("/usr/lib64/python2.6/site-packages/")
import htcondor
import numpy as np

from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config

from ROOT import TFile

class MiniSkim(Task):    
    def __init__(self, config = {}):
        Task.__init__(self, config)
        self._allowParallel = False

        self.badStatus = ["failed"]

    def __crabConfig(self):
        ##Crab config
        self.crabConf = config()

        self.crabConf.General.requestName = "MiniSkim_{}".format(self["datasetName"])
        self.crabConf.General.workArea = self["dir"]
        self.crabConf.General.transferOutputs = True
        self.crabConf.General.transferLogs = False

        self.crabConf.JobType.pluginName = "Analysis"
        self.crabConf.JobType.psetName = "ChargedHiggs/Skimming/python/miniskimmer.py"
        self.crabConf.JobType.pyCfgParams = ["outname={}.root".format(self["datasetName"])]
        self.crabConf.JobType.outputFiles = ["{}.root".format(self["datasetName"])]

        self.crabConf.Data.inputDataset = self["dataset"] 
        self.crabConf.Data.inputDBS = "global"
        self.crabConf.Data.splitting = "EventAwareLumiBased"
        self.crabConf.Data.unitsPerJob = 400000
        self.crabConf.Data.outLFNDirBase = "/store/user/dbrunner/skim"
        self.crabConf.Site.storageSite = "T2_DE_DESY"

    def status(self):
        ##If no job submitted yet, no status will be checked
        if not os.path.exists(self["crab-dir"]):
            return None

        ##Get status
        crabStatus = crabCommand("status", dir=self["crab-dir"])
        nFinished = 0
    
        ##If submission failed, resubmit and leave function
        if crabStatus["status"] == "SUBMITFAILED":
            self.__crabConfig()
            os.system("command rm -r {}".format(self["crab-dir"]))
            crabCommand('submit', config = self.crabConf)

            return None

        ##Set output here, because crab is splitting jobs and #jobs is not know beforehand
        if self["output"] == []:
            for jobNr in crabStatus["jobs"].keys():
                self["output"].append("{}/results/{}_{}.root".format(self["crab-dir"], self["name"], jobNr))

        for status in crabStatus["jobsPerStatus"].keys():
            ##If one jobs failed, resubmit and leave loop
            if status in self.badStatus and (crabStatus["dbStatus"] == "SUBMITTED" or crabStatus["dbStatus"] == "RESUBMITFAILED"):
                if crabStatus["jobsPerStatus"][status] != 0:
                    ##Check if you have to increase memory/run time
                    exitCode = [crabStatus["jobs"][key]["Error"][0] for key in crabStatus["jobs"] if "Error" in crabStatus["jobs"][key]]

                    runTime = "1315" if not 50664 in exitCode else "1600"
                    memory = "2000" if not 50660 in exitCode else "2500"

                    crabCommand("resubmit", dir=self["crab-dir"], maxmemory=memory, maxjobruntime=runTime)
                    break

        if "finished" in crabStatus["jobsPerStatus"]:
            if crabStatus["jobsPerStatus"]["finished"] == len(crabStatus["jobs"].keys()):
                jobsNrs = crabStatus["jobs"].keys()

                ##Skip if output already retrieved
                if(len(os.listdir("{}/results/".format(self["crab-dir"]))) == len(jobsNrs)):
                    self["status"] = "FINISHED"
                    return None

                ##Retrieve output
                for jobList in [jobsNrs[i:i + 500] for i in range(0, len(jobsNrs), 500)]:
                    while(True):
                        try:
                            crabCommand("getoutput", dir=self["crab-dir"], jobids=",".join(jobList))
                            break
                        except:
                            pass       
     
                ##If all jobs finished, retrieve output
                self["status"] = "FINISHED"
            
        time.sleep(60)

    def run(self):
        self.__crabConfig()
        if not os.path.exists(self["crab-dir"]):
            crabCommand('submit', config = self.crabConf)
      
    def output(self):
        self["output"] = []

    @staticmethod
    def configure(fileList):
        datasets = []
        tasks = []

        with open(fileList) as f:
            datasets = [dataset for dataset in f.read().splitlines() if dataset != ""]

        for dataset in datasets:
            if "SIM" in dataset:
                if "ext" in dataset:
                    name = dataset.split("/")[1] + "_ext"
                else: 
                    name = dataset.split("/")[1]
            
            else: 
                name = dataset.split("/")[1] + "_" + dataset.split("/")[2]

            config = {
                        "name": "MiniSkim_{}".format(name),
                        "display-name": "MINI Skim: {}".format(name.split("_")[0]),
                        "dir": "{}/Skimdir/masterSkim/{}".format(os.environ["HOME2"], name),
                        "crab-dir": "{}/Skimdir/masterSkim/{}/crab_MiniSkim_{}".format(os.environ["HOME2"], name, name),
                        "dataset": dataset,                        
                        "datasetName": name, 
            }

            tasks.append(MiniSkim(config))

        return tasks
