#ifndef WEIGHTANALYZER_H
#define WEIGHTANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <TParameter.h>

#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

class WeightAnalyzer : public BaseAnalyzer {
    private:
        //Bool for data
        bool isData;
        std::string era;    

        std::shared_ptr<Token> tokens;

        //xSec information
        float xSec, lumi, lumiUp, lumiDown, nGen = 0., nGenWeighted = 0.;
        
        std::string pileUpFile;
    
        //Values for branch
        float genWeight = 1., eventNumber = 1., nTrueInt = 1.;
        float pdfWeights[102], scaleWeights[8];
        std::vector<float> prefireWeights;

        //Histograms
        std::shared_ptr<TH1F> puMC;

        //TTreeReader Values
        std::unique_ptr<TTreeReaderValue<float>> nPU, genWeightValue;
        std::unique_ptr<TTreeReaderValue<ULong64_t>> evtNumber;

    public:
        WeightAnalyzer(const std::string& era, const float& xSec, TTreeReader &reader);
        WeightAnalyzer(const std::string& era, const float& xSec, const std::shared_ptr<Token>& tokens);

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
