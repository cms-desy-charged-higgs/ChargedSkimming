#!/usr/bin/env python

import argparse
import yaml
import os

from dbs.apis.dbsClient import DbsApi

def parser():
    parser = argparse.ArgumentParser(description = "Create list of datasets", formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("out_dir", metavar='out-dir', type=str, action="store", help="Working directory")
    parser.add_argument("--wildcard", nargs='+', required = True, help = "Wildcard used with dasgoclient")

    return parser.parse_args()

def main():
    ##Parser arguments
    args = parser()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    processFile = yaml.load(open("{}/src/ChargedSkimming/Skimming/data/xsec.yaml".format(os.environ["CMSSW_BASE"])), Loader=yaml.Loader)
    listOfProcess = [proc for proc in processFile if processFile[proc]["isUsed"] and not "HPlus" in proc]

    url="https://cmsweb.cern.ch/dbs/prod/{}/DBSReader".format("global")
    api=DbsApi(url=url)

    listOfDasesets = []

    for proc in listOfProcess:
        dataSet = ""

        for card in args.wildcard:
            searchString = "/{}*/{}/MINIAOD*".format(proc, card)
            dataSets = [d["dataset"] for d in api.listDatasets(dataset = searchString)]
            dataSets.sort(key = lambda s: len(s))

            if len(dataSets) != 0:
                dataSet = dataSets[0]
                break

        listOfDasesets.append(dataSet if dataSet else "Not found: {}".format(proc))

    with open("{}/filelist_bkg_MINI.yaml".format(args.out_dir), "w") as f:
        f.write("##Created with createDataList.py with following wildcards: {}\n\n".format(args.wildcard))

        listOfDasesets.sort()
        yaml.dump(listOfDasesets, f, default_flow_style = False)

if __name__ == "__main__":
    main()
