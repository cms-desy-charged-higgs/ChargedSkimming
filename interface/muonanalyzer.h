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
        //Check if data
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

        //Kinematic cut criteria
        int era;
        float ptCut;
        float etaCut;

        //Vector with output varirables of the output tree
        std::vector<std::string> floatNames;
        std::vector<std::string> boolNames;

        std::vector<std::vector<float>> floatVariables;
        std::vector<std::vector<bool>> boolVariables;

        //EDM Token for MINIAOD analysis
        muToken muonToken;
        trigObjToken triggerObjToken; 
        genPartToken genParticleToken;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> muonPt;
        std::unique_ptr<TTreeReaderArray<float>> muonEta;
        std::unique_ptr<TTreeReaderArray<float>> muonPhi;
        std::unique_ptr<TTreeReaderArray<float>> muonIso;
        std::unique_ptr<TTreeReaderArray<int>> muonCharge;
        std::unique_ptr<TTreeReaderArray<bool>> muonMediumID;
        std::unique_ptr<TTreeReaderArray<bool>> muonTightID;
        std::unique_ptr<TTreeReaderArray<int>> muonGenIdx;

    public:
        MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader);
        MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, muToken &muonToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken);

        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
