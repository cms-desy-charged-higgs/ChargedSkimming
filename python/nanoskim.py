from ChargedHiggs.Workflow.task import Task

import os
import sys
sys.path.append("/usr/lib64/python2.6/site-packages/")
import htcondor

from ROOT import TFile

class NanoSkim(Task):    
    def __init__(self, config = {}):
        Task.__init__(self, config)

        ##Condor submission config
        self["condor"] = {}

        ##General configuration
        self["CMSSW_DIR"] = os.environ["CMSSW_BASE"] + "/src/"
        
    def __write_condor(self):
        self["condor"]["executable"] = self["CMSSW_DIR"] + "ChargedHiggs/nano_skimming/batch/produceSkim.sh"
        self["condor"]["arguments"] = " ".join([self["input_file"], self["name"] + ".root"] + self["channel"])
        self["condor"]["universe"] = "vanilla"

        self["condor"]["should_transfer_files"] ="YES"
        self["condor"]["transfer_input_files"] = ",".join([self["CMSSW_DIR"] + "ChargedHiggs", self["CMSSW_DIR"] + "x509"])

        self["condor"]["on_exit_hold"] = "(ExitBySignal == True) || (ExitCode != 0)"  
        self["condor"]["periodic_release"] = "(NumJobStarts < 100) && ((CurrentTime - EnteredCurrentStatus) > 60)"

        self["condor"]["log"] = "{}/{}.log".format(self["dir"], self["name"])
        self["condor"]["output"] = "{}/{}.out".format(self["dir"], self["name"])
        self["condor"]["error"] = "{}/{}.err".format(self["dir"], self["name"])

        self["condor"]["when_to_transfer_output"] = "ON_EXIT"
        self["condor"]["transfer_output_remaps"] = '"{}.root = {}/{}.root"'.format(self["dir"].split("/")[-1], self["dir"], self["name"])

    def status(self):
        if os.path.isfile("{}/{}.err".format(self["dir"], self["name"])):
            if os.path.isfile(self["output"]):
                self["status"] = "FINISHED"

            else:
                self["status"] = "FAILED"

    def run(self):
        self.__write_condor()

        ##Create proxy    
        os.system("chmod 755 {}".format(os.environ["X509_USER_PROXY"])) 
        os.system("cp -u {} {}".format(os.environ["X509_USER_PROXY"], self["CMSSW_DIR"])) 

        ##Agressively submit your jobs    
        job = htcondor.Submit()
        job.update(self["condor"])
        schedd = htcondor.Schedd()
      
        while(True):
            try: 
                with schedd.transaction() as txn:
                    job.queue(txn)
                    print "Submit job for file {}".format(self["input_file"])

                break    

            except:
                pass

    def output(self):
        self["output"] = "{}/{}.root".format(self["dir"], self["name"])
