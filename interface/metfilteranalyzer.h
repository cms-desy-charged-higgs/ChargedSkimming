#ifndef METFILTERANALYZER_H
#define METFILTERANALYZER_H

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>

#include <numeric>

class MetFilterAnalyzer : public BaseAnalyzer {
    private:
        //Era
        int era;

        //Map with TTreeReaderValues for each era
        std::map<int, std::vector<std::string>> filterNames;

        //Vector with TTreeReaderValues
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> filterValues;

    public:
        MetFilterAnalyzer(const int &era, TTreeReader &reader);
        void BeginJob(TTree *tree, bool &isData);
        bool Analyze(std::pair<TH1F*, float> &cutflow);
        void EndJob(TFile* file);
};

#endif
