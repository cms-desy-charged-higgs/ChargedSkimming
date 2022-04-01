#!/usr/bin/env python

import glob
import os

def main():
    inDir = "{}/src/jsonpog-integration/POG/*/*/*gz".format(os.environ["CMSSW_BASE"])

    for d in glob.glob(inDir):
        print(d)
        os.system("gunzip {}".format(d))

if __name__ == "__main__":
    main()
