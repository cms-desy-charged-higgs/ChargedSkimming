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
        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
