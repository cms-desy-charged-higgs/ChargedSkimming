from ChargedAnalysis.Workflow.task import Task

import os
import sys

from ROOT import TFile

class SkimHadd(Task):    
    def __init__(self, config = {}):
        Task.__init__(self, config)

    def status(self):
        if os.path.isfile(self["output"]):
            skimmedFile = TFile.Open(self["output"])
            
            if(skimmedFile.GetSize() == -1.):
                self["status"] = "FAILED"

            else:
                self["status"] = "FINISHED"

    def run(self):
        os.system("hadd -f -j 16 {} {}".format(self["output"], " ".join(self["dependent_files"])))

    def output(self):
        self["output"] = "{}/{}.root".format(self["dir"], self["name"][5:])

    @staticmethod
    def configure(skimTasks, isNANO):
        tasks = []

        if not isNANO:
            for skim in skimTasks:
                config = {
                        "name": "{}".format(skim["name"].replace("MiniSkim", "Hadd")),
                        "display-name": "Hadd: {}".format(skim["name"].split("_")[1]),
                        "dir": " {}/merged".format(skim["dir"]),
                        "dependencies": [skim["name"]],
                }

                tasks.append(SkimHadd(config))

        return tasks
