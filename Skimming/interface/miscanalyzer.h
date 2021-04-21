#ifndef MISCANALYZER_H
#define MISCANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <DataFormats/PatCandidates/interface/IsolatedTrack.h>
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>

typedef edm::EDGetTokenT<std::vector<pat::IsolatedTrack>> isoToken;
typedef edm::EDGetTokenT<LHEEventProduct> lheToken;

class MiscAnalyzer: public BaseAnalyzer{
    private:
        //Check if data
        bool isData;
        bool isSig = false;

        int era;
        float etaCut;

        isoToken isoTrackToken;
        lheToken lheTok;

        //Vector with output varirables of the output tree
        std::map<std::string, float*> floatVar;
        std::map<std::string, short*> intVar;

        float Pt[20], Eta[20], Phi[20], Isolation[20], Dxy[20], Dz[20];
        short Charge[20], PDG[20], FromPV[20];

        short nTracks, nParton;

        int GetNParton(const edm::Event* event);

    public:
        MiscAnalyzer(const int& era, const float& etaCut, isoToken& isoTrackToken, lheToken& lheTok);

        void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
