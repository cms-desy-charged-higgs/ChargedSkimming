#ifndef GENPARTANALYZER_H
#define GENPARTANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

class GenPartAnalyzer: public BaseAnalyzer{
    private:
        //Bool for checking if data file
        bool isData;

        //EMD Token for MINIAOD
        genPartToken genParticleToken;

        //Set output names
        std::vector<std::string> floatNames;
        std::vector<std::vector<float>> leptonVariables;
        std::vector<std::vector<float>> h1Variables;
        std::vector<std::vector<float>> h2Variables;

    public:
        GenPartAnalyzer(genPartToken& genParticleToken);
        GenPartAnalyzer(TTreeReader &reader);
        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
