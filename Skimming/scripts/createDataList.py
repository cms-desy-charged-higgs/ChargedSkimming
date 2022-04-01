#!/usr/bin/env python3

import yaml
import os
import subprocess

def getDataSets(listOfProcess, wildcards):
    listOfDasesets = []

    for proc in listOfProcess:
        dataSet = ""

        print("Process: {}".format(proc))

        for card in wildcards:        
            searchString = "/{}*/{}/NANOAOD*".format(proc, card)
            dataSets = subprocess.check_output("dasgoclient --query='dataset={}'".format(searchString), shell = True).decode('ascii').split("\n")[:-1]
            dataSets.sort(key = lambda s: len(s))

            if len(dataSets) != 0:
                dataSet = dataSets[0]
                break

        listOfDasesets.append(dataSet if dataSet else "Not found: {}".format(proc))

    return listOfDasesets

def writeFile(listOfDasesets, wildcards, outDir, type):
    with open("{}/filelist_{}_NANO.yaml".format(outDir, type), "w") as f:
        f.write("##Created with createDataList.py with following wildcards: {}\n\n".format(wildcards))

        listOfDasesets.sort()
        yaml.dump(listOfDasesets, f, default_flow_style = False)


def main():
    processFile = yaml.load(open("{}/src/ChargedSkimming/Skimming/data/xsec.yaml".format(os.environ["CMSSW_BASE"])), Loader=yaml.Loader)
    wildCards = yaml.load(open("{}/src/ChargedSkimming/Skimming/data/wildCards.yaml".format(os.environ["CMSSW_BASE"])), Loader=yaml.Loader)

    for era in wildCards.keys():
        outDir = "{}/src/ChargedSkimming/Skimming/data/filelists/{}/UL".format(os.environ["CMSSW_BASE"], era)

        if not os.path.exists(outDir):
            os.makedirs(outDir)

        sigProcesses = [proc for proc in processFile if processFile[proc]["isUsed"] and "HPlus" in proc]
        bkgProcesses = [proc for proc in processFile if processFile[proc]["isUsed"] and not "HPlus" in proc]
        data = ["SingleMuon", "SingleElectron" if era != 2018 else "EGamma"]
    
        sigSets = getDataSets(sigProcesses, wildCards[era]["MC"])
        writeFile(sigSets, wildCards[era]["MC"], outDir, "sig")

        bkgSets = getDataSets(bkgProcesses, wildCards[era]["MC"])
        writeFile(bkgSets, wildCards[era]["MC"], outDir, "bkg")

        dataSets = []
    
        for e in wildCards[era]["data"]["eras"]:
            wildcards = [w.replace("*", "*{}*".format(e), 1) for w in wildCards[era]["data"]["wildcards"]]
            dataSets.extend(getDataSets(data, wildcards))

        writeFile(dataSets, wildCards[era]["data"]["wildcards"], outDir, "data")

if __name__ == "__main__":
    main()
