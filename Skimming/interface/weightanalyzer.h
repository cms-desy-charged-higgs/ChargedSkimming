#ifndef WEIGHTANALYZER_H
#define WEIGHTANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

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
        float eventNumber = 1.;
        float nTrueInt = 1.;
        std::vector<float> prefireWeights;
        std::vector<float> pdfWeights;
        std::vector<float> scaleWeights;

        //Lumi information
        std::map<int, float> lumis;

        //Token for MINIAOD
        puToken pileupToken;
        genToken geninfoToken;
        std::vector<edm::EDGetTokenT<double>> prefireTokens;
        wgtToken scaleToken;
        wgtToken pdfToken;

        //Histograms
        TH1F* puMC; 
        TH1F* nGenHist;
        TH1F* nGenWeightedHist;

        //TTreeReader Values
        std::unique_ptr<TTreeReaderValue<float>> nPU;
        std::unique_ptr<TTreeReaderValue<float>> genWeightValue;
        std::unique_ptr<TTreeReaderValue<ULong64_t>> evtNumber;

    public:
        WeightAnalyzer(const float era, const float xSec, TTreeReader &reader);
        WeightAnalyzer(const float era, const float xSec, puToken &pileupToken, genToken &geninfoToken, const std::vector<edm::EDGetTokenT<double>> prefireTokens, wgtToken& pdfToken, wgtToken& scaleToken);
        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
