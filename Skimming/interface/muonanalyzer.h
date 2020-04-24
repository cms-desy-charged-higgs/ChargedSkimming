#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <DataFormats/PatCandidates/interface/Muon.h>

typedef edm::EDGetTokenT<std::vector<pat::Muon>> muToken;

class MuonAnalyzer: public BaseAnalyzer{
    private:
        //Check if data
        bool isData;

        //Map for SF files
        std::map<int, std::string> isoSFfiles, triggerSFfiles, IDSFfiles;

        //Hist with scale factors
        TH2F* triggerSFhist;
        std::vector<TH2F*> IsoHist, IDHist;

        //Kinematic cut criteria
        int era;
        float ptCut, etaCut;

        //Vector with output varirables of the output tree
        std::map<std::string, std::vector<float>&> floatVar;
        std::map<std::string, std::vector<char>&> intVar;

        std::vector<float> Pt, Eta, Phi, Isolation, triggerSF, triggerSFUp, triggerSFDown, looseSF, looseSFUp, looseSFDown, mediumSF, mediumSFUp, mediumSFDown, tightSF, tightSFUp, tightSFDown, looseIsoLooseSF, looseIsoLooseSFDown, looseIsoLooseSFUp, looseIsoMediumSF, looseIsoMediumSFDown, looseIsoMediumSFUp, looseIsoTightSF, looseIsoTightSFDown, looseIsoTightSFUp, tightIsoMediumSF, tightIsoMediumSFDown, tightIsoMediumSFUp, tightIsoTightSF, tightIsoTightSFDown, tightIsoTightSFUp;

        std::vector<char> ID, isoID, Charge, isFromHPlus;

        char nMuons;

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
