from ChargedHiggs.Workflow.task import Task

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
        os.system("hadd -j 16 {} {}".format(self["output"], " ".join(self["dependent_files"])))

    def output(self):
        self["output"] = "{}/{}.root".format(self["dir"], self["name"])



