#ifndef GENPARTANALYZER_H
#define GENPARTANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

class GenPartAnalyzer: public BaseAnalyzer{
    private:
        //Bool for checking if data file
        bool isData;

        //EMD Token for MINIAOD
        genPartToken genParticleToken;
        std::vector<std::string> particleNames;
        std::vector<int> particleIDs;

        //Set output names
        std::map<std::string, float*> floatVar;
        std::map<std::string, short*> intVar;

        float Pt[100], Eta[100], Phi[100], Mass[100];
        short partID[100], mothID[100];

        short nGenPart;

    public:
        GenPartAnalyzer(genPartToken& genParticleToken, const std::vector<std::string>& particleNames);
        GenPartAnalyzer(TTreeReader &reader, const std::vector<std::string>& particleNames);
        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
