#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

//Muon class to be safed in tree
struct Muon {
    TLorentzVector fourVec;
    Bool_t isMedium;
    Bool_t isTight;
    Bool_t isLooseIso;
    Bool_t isTightIso;

    Float_t looseIsoMediumSF = 1.;
    Float_t tightIsoMediumSF = 1.;
    Float_t looseIsoTightSF = 1.;
    Float_t tightIsoTightSF = 1.;
    Float_t mediumSF = 1.;
    Float_t tightSF = 1.;
    Float_t triggerSF = 1.;

    Bool_t isTriggerMatched;

    TLorentzVector genVec;
    Bool_t isgenMatched = false;
    Bool_t isFromHc = false;

    Float_t charge;
};

class MuonAnalyzer: public BaseAnalyzer{
    private:
        //Bool for checking if data file
        bool isData;

        //Map for SF files
        std::map<int, std::string> isoSFfiles;
        std::map<int, std::string> triggerSFfiles;
        std::map<int, std::string> IDSFfiles;

        //Hist with scale factors
        TH2F* triggerSFhist;

        TH2F* looseIsoMediumSFhist;
        TH2F* tightIsoMediumSFhist;
        TH2F* looseIsoTightSFhist;
        TH2F* tightIsoTightSFhist;

        TH2F* mediumIDHist;
        TH2F* tightIDHist;

        //Input for selecting muon
        int era;
        float ptCut;
        float etaCut;
        unsigned int minNMuon;

        //Pat vector for Mini AOD analysis
        edm::Handle<std::vector<pat::Muon>> muons;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> muonPt;
        std::unique_ptr<TTreeReaderArray<float>> muonEta;
        std::unique_ptr<TTreeReaderArray<float>> muonPhi;
        std::unique_ptr<TTreeReaderArray<float>> muonIso;
        std::unique_ptr<TTreeReaderArray<int>> muonCharge;
        std::unique_ptr<TTreeReaderArray<bool>> muonMediumID;
        std::unique_ptr<TTreeReaderArray<bool>> muonTightID;
        std::unique_ptr<TTreeReaderArray<int>> muonGenIdx;

        //Valid electron collection
        std::vector<Muon> validMuons;

        //Set Gen particle information
        void SetGenParticles(Muon &validMuon, const int &i);

    public:
        MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNMuon, TTreeReader& reader);
        MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNMuon, edm::Handle<std::vector<pat::Muon>> muons);

        void BeginJob(TTree* tree, bool &isData);
        bool Analyze(std::pair<TH1F*, float> &cutflow);
        void EndJob(TFile* file);
};

#endif
