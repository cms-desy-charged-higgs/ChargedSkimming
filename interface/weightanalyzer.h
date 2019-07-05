#ifndef WEIGHTANALYZER_H
#define WEIGHTANALYZER_H

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>

class WeightAnalyzer : public BaseAnalyzer {
    private:
        //Bool for data
        bool isData;
        float era;    

        //xSec information
        float xSec;
    
        //Values for branch
        float lumi = 1.;
        float genWeight = 1.;
        float puWeight = 1.;
        float eventNumber = 1.;

        //Lumi information
        std::map<int, float> lumis;

        //Information for pile Up reweighting
        std::map<int, std::string> pileUpFiles;
        TH1F* pileUpWeights; 

        TH1F* nGenHist;

        //TTreeReader Values
        std::unique_ptr<TTreeReaderValue<float>> nPU;
        std::unique_ptr<TTreeReaderValue<float>> genWeightValue;
        std::unique_ptr<TTreeReaderValue<ULong64_t>> evtNumber;


    public:
        WeightAnalyzer(const float era, const float xSec, TTreeReader &reader);
        void BeginJob(TTree *tree, bool &isData);
        bool Analyze(std::pair<TH1F*, float> &cutflow);
        void EndJob(TFile* file);
};

#endif
