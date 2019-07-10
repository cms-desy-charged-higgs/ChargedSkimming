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
        //Bool for checking if data file
        bool isData;

        //Map for SF files
        std::map<int, std::string> mediumSFfiles;
        std::map<int, std::string> tightSFfiles;
        std::map<int, std::string> recoSFfiles;

        //Hist with scale factors
        TH2F* mediumSFhist;
        TH2F* tightSFhist;
        TH2F* recoSFhist;

        //Input for selecting electrons
        int era;
        float ptCut;
        float etaCut;
        unsigned int minNEle;

        //EDM Token for MINIAOD analysis
        eToken eleToken;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> elePt;
        std::unique_ptr<TTreeReaderArray<float>> eleEta;
        std::unique_ptr<TTreeReaderArray<float>> elePhi;
        std::unique_ptr<TTreeReaderArray<float>> eleIso;
        std::unique_ptr<TTreeReaderArray<int>> eleCharge;
        std::unique_ptr<TTreeReaderArray<bool>> eleMediumMVA;
        std::unique_ptr<TTreeReaderArray<bool>> eleTightMVA;
        std::unique_ptr<TTreeReaderArray<int>> eleGenIdx;

        //Valid electron collection
        std::vector<Electron> validElectrons;

        //Set Gen particle information
        void SetGenParticles(Electron &validElectron, const int &i);

    public:
        ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, edm::EDGetTokenT<std::vector<pat::Electron>>& eleToken);
        ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, TTreeReader& reader);

        void BeginJob(TTree* tree, bool &isData);
        bool Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
