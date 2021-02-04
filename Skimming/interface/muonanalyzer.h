#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <DataFormats/PatCandidates/interface/Muon.h>

typedef edm::EDGetTokenT<std::vector<pat::Muon>> muToken;

class MuonAnalyzer: public BaseAnalyzer{
    private:
        //Check if data
        bool isData;
        bool isSig = false;

        //Hist with scale factors
        TH2F* triggerSFhist;
        std::vector<TH2F*> IsoHist, IDHist;

        //Kinematic cut criteria
        int era;
        float ptCut, etaCut;

        //Vector with output varirables of the output tree
        std::map<std::string, float*> floatVar;
        std::map<std::string, short*> intVar;

        float Pt[20], Eta[20], Phi[20], Isolation[20], triggerSF[20], triggerSFUp[20], triggerSFDown[20], looseSF[20], looseSFUp[20], looseSFDown[20], mediumSF[20], mediumSFUp[20], mediumSFDown[20], tightSF[20], tightSFUp[20], tightSFDown[20], looseIsoLooseSF[20], looseIsoLooseSFDown[20], looseIsoLooseSFUp[20], looseIsoMediumSF[20], looseIsoMediumSFDown[20], looseIsoMediumSFUp[20], looseIsoTightSF[20], looseIsoTightSFDown[20], looseIsoTightSFUp[20], tightIsoMediumSF[20], tightIsoMediumSFDown[20], tightIsoMediumSFUp[20], tightIsoTightSF[20], tightIsoTightSFDown[20], tightIsoTightSFUp[20];

        short ID[20], isoID[20], Charge[20], partID[20], mothID[20], grandID[20];

        short nMuons;

        //EDM Token for MINIAOD analysis
        muToken muonToken;
        genPartToken genParticleToken;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> muonPt, muonEta, muonPhi, muonIso;
        std::unique_ptr<TTreeReaderArray<int>> muonCharge;
        std::unique_ptr<TTreeReaderArray<bool>> muonLooseID, muonMediumID, muonTightID;
        std::unique_ptr<TTreeReaderArray<int>> muonGenIdx;

    public:
        MuonAnalyzer(const int& era, const float& ptCut, const float& etaCut, TTreeReader& reader);
        MuonAnalyzer(const int& era, const float& ptCut, const float& etaCut, muToken& muonToken, genPartToken& genParticleToken);

        void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
