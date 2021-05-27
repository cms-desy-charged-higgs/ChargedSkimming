#ifndef MISCANALYZER_H
#define MISCANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <DataFormats/PatCandidates/interface/IsolatedTrack.h>
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>

class MiscAnalyzer: public BaseAnalyzer{
    private:
        //Check if data
        bool isData = false;
        bool isSig = false;

        std::string era;
        std::shared_ptr<Token> tokens;
        float etaCut;

        //Vector with output varirables of the output tree
        std::map<std::string, float*> floatVar;
        std::map<std::string, short*> intVar;

        float Pt[20], Eta[20], Phi[20], Isolation[20], Dxy[20], Dz[20];
        short Charge[20], PDG[20], FromPV[20];

        short nTracks, nParton;

        int GetNParton(const edm::Event* event);

    public:
        MiscAnalyzer(const std::string& era, const std::shared_ptr<Token>& tokens);

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
