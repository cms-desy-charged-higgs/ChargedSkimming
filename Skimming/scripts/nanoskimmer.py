#!/usr/bin/env python

import argparse
import os
import yaml

from ROOT import NanoSkimmer
from ROOT.std import vector

def parser():
    parser = argparse.ArgumentParser(description = "Script to skim NANOAOD root files", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--filename", action = "store", help = "Name of the input NANOAOD root file")
    parser.add_argument("--channel", nargs="+", choices = ["mu4j", "e4j", "mu2j1f", "e2j1f", "mu2f", "e2f"], default = ["mu4j", "e4j", "mu2j1f", "e2j1f", "mu2f", "e2f"], help = "Final state which is of interest")

    parser.add_argument("--out-dir", type = str, default = "{}/src".format(os.environ["CMSSW_BASE"]), help = "Name of output directory")    
    parser.add_argument("--out-name", type = str, default = "outputSkim.root", help = "Output name of skimmed file")

    return parser.parse_args()

def main():
    args = parser()

    xSecFile = yaml.load(file("{}/src/ChargedHiggs/Skimming/data/xsec.yaml".format(os.environ["CMSSW_BASE"]), "r"))

    xSec = 1.

    for key in xSecFile.keys():
        if key in args.out_name:
            xSec = xSecFile[key]["xsec"]

    isData = True in [name in args.out_name for name in ["Electron", "Muon", "MET"]]

    channels = vector("string")()
    [channels.push_back(channel) for channel in args.channel]

    skimmer = NanoSkimmer(args.filename, isData)
    skimmer.EventLoop(channels, xSec)
    skimmer.WriteOutput(args.out_name)

if __name__ == "__main__":
    main()
