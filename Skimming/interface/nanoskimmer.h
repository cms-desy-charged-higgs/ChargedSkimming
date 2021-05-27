#ifndef NANOSKIMMER_H
#define NANOSKIMMER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <vector>
#include <string>
#include <chrono>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

class NanoSkimmer{
    private:
        //Measure execution time
        std::chrono::steady_clock::time_point start;
        std::chrono::steady_clock::time_point end;

        //Input
        std::string inFile;
        std::string outFile;
        std::vector<std::string> channels;
        float xSec;
        std::string era;
        bool isData;    

        int nEvents = 0; 

        //Output file
        std::vector<TFile*> outputFiles;

        //Trees for each systematic analysis
        std::vector<std::vector<TTree*>> outputTrees;
        std::vector<std::vector<CutFlow>> cutflows; 

        //Vector with wished analyzers for each systematic
        std::vector<std::vector<std::shared_ptr<BaseAnalyzer>>> analyzers;

        //Progress bar function
        void ProgressBar(const int &progress);

        //Configure analysis modules
        void Configure(TTreeReader& reader);
        

    public:
        NanoSkimmer();
        NanoSkimmer(const std::string &inFile, const std::string &outFile, const std::vector<std::string>& channels, const float& xSec, const std::string& era, const bool &isData);
        void EventLoop();
        void WriteOutput(); 
};

#endif
