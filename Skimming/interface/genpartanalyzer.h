#ifndef GENPARTANALYZER_H
#define GENPARTANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/tokens.h>
#include <ChargedSkimming/Skimming/interface/util.h>

class GenPartAnalyzer: public BaseAnalyzer{
    private:
        //Bool for checking if data file
        bool isData;
        std::shared_ptr<Token> tokens;

        //EMD std::shared_ptr<Token> for MINIAOD
        std::vector<std::string> particleNames;
        std::vector<int> particleIDs;

        //Set output names
        std::map<std::string, float*> floatVar;
        std::map<std::string, short*> intVar;

        float Pt[100], Eta[100], Phi[100], Mass[100];
        short partID[100], mothID[100];

        short nGenPart;

    public:
        GenPartAnalyzer(const std::shared_ptr<Token>& tokens);
        GenPartAnalyzer(TTreeReader& reader);

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
