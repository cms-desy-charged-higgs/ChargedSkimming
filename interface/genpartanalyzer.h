#ifndef GENPARTANALYZER_H
#define GENPARTANALYZER_H

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>

//Jet class to be safed in tree
struct GenPart {
    TLorentzVector Hc;
    TLorentzVector h1;
    TLorentzVector h2;
    TLorentzVector W;
};


class GenPartAnalyzer: public BaseAnalyzer{
    private:
        //Bool for checking if data file
        bool isData;

        //Gen particle vector
        GenPart genColl;

        //EMD Token for MINIAOD
        genPartToken genParticleToken;

    public:
        GenPartAnalyzer(genPartToken& genParticleToken);
        GenPartAnalyzer(TTreeReader &reader);
        void BeginJob(TTree* tree, bool &isData);
        bool Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
