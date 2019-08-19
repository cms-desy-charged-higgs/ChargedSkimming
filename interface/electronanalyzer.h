#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <ChargedHiggs/Skimming/interface/baseanalyzer.h>

//Electron class to be safed in tree
struct Electron {
    TLorentzVector fourVec;
    Bool_t isMedium;
    Bool_t isTight;

    Bool_t isTriggerMatched;

    Float_t recoSF = 1.;
    Float_t mediumMvaSF = 1.;
    Float_t tightMvaSF = 1.;

    Float_t charge;
    Float_t isolation; 

    TLorentzVector genVec;
    Bool_t isgenMatched = false;
    Bool_t isFromHc = false;
};


class ElectronAnalyzer: public BaseAnalyzer {
    private:
        //Check if data
        bool isData;

        //Map for SF files
        std::map<int, std::string> mediumSFfiles;
        std::map<int, std::string> tightSFfiles;
        std::map<int, std::string> recoSFfiles;

        //Hist with scale factors
        TH2F* mediumSFhist;
        TH2F* tightSFhist;
        TH2F* recoSFhist;

        //Kinematic cut criteria
        int era;
        float ptCut;
        float etaCut;

        //EDM Token for MINIAOD analysis
        eToken eleToken;
        trigObjToken triggerObjToken;
        genPartToken genParticleToken;

        //Vector with output varirables of the output tree
        std::vector<std::string> floatNames;
        std::vector<std::string> boolNames;

        std::vector<std::vector<float>> floatVariables;
        std::vector<std::vector<bool>> boolVariables;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> elePt;
        std::unique_ptr<TTreeReaderArray<float>> eleEta;
        std::unique_ptr<TTreeReaderArray<float>> elePhi;
        std::unique_ptr<TTreeReaderArray<float>> eleIso;
        std::unique_ptr<TTreeReaderArray<int>> eleCharge;
        std::unique_ptr<TTreeReaderArray<bool>> eleMediumMVA;
        std::unique_ptr<TTreeReaderArray<bool>> eleTightMVA;

    public:
        ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, eToken& eleToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken);
        ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader);

        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
