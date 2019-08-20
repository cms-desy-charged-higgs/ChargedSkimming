#ifndef METFILTERANALYZER_H
#define METFILTERANALYZER_H

#include <ChargedAnalysis/Skimming/interface/baseanalyzer.h>

#include <numeric>

class MetFilterAnalyzer : public BaseAnalyzer {
    private:
        //Era
        int era;

        //Map with TTreeReaderValues for each era
        std::map<int, std::vector<std::string>> filterNames;

        //EDM Token for MINIAOD analysis
        trigToken triggerToken;

        //Vector with TTreeReaderValues
        std::vector<std::unique_ptr<TTreeReaderValue<bool>>> filterValues;

    public:
        MetFilterAnalyzer(const int &era, TTreeReader &reader);
        MetFilterAnalyzer(const int &era, trigToken& triggerToken);
        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
