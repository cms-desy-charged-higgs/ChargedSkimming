#ifndef WEIGHTANALYZER_H
#define WEIGHTANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

typedef edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken;
typedef edm::EDGetTokenT<GenEventInfoProduct> genToken;
typedef edm::EDGetTokenT<std::vector<float>> wgtToken;

class WeightAnalyzer : public BaseAnalyzer {
    private:
        //Bool for data
        bool isData;
        float era;    

        //xSec information
        float xSec;
    
        //Values for branch
        float genWeight = 1., eventNumber = 1., nTrueInt = 1.;
        std::vector<float> prefireWeights, pdfWeights, scaleWeights;

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
        std::unique_ptr<TTreeReaderValue<float>> nPU, genWeightValue;
        std::unique_ptr<TTreeReaderValue<ULong64_t>> evtNumber;

    public:
        WeightAnalyzer(const float& era, const float& xSec, TTreeReader &reader);
        WeightAnalyzer(const float& era, const float& xSec, puToken &pileupToken, genToken &geninfoToken, const std::vector<edm::EDGetTokenT<double>>& prefireTokens, wgtToken& pdfToken, wgtToken& scaleToken);
        void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
