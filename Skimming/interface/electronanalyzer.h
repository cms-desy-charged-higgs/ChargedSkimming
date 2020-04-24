#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <DataFormats/PatCandidates/interface/Electron.h>

typedef edm::EDGetTokenT<std::vector<pat::Electron>> eToken;

class ElectronAnalyzer: public BaseAnalyzer {
    private:
        //Check if data
        bool isData;

        //Map for SF files
        std::map<int, std::string> looseSFfiles, mediumSFfiles, tightSFfiles, recoSFfiles;

        //Hist with scale factors
        TH2F* looseSFhist;
        TH2F* mediumSFhist;
        TH2F* tightSFhist;
        TH2F* recoSFhist;

        //Kinematic cut criteria
        int era;
        float ptCut, etaCut;

        //EDM Token for MINIAOD analysis
        eToken eleToken;
        genPartToken genParticleToken;

        //Name of the energy correction (dependent on systematic study)
        std::string energyCorrection;

        //Vector with output varirables of the output tree
        std::map<std::string, std::vector<float>&> floatVar;
        std::map<std::string, std::vector<char>&> intVar;

        std::vector<float> Pt, Eta, Phi, recoSF, recoSFUp, recoSFDown, looseSF, looseSFUp, looseSFDown, mediumSF, mediumSFUp, mediumSFDown, tightSF, tightSFUp, tightSFDown, Isolation;

        std::vector<char> ID, Charge, isFromHPlus;

        char nElectrons;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> elePt, eleEta, elePhi, eleIso;
        std::unique_ptr<TTreeReaderArray<int>> eleCharge, eleID;

    public:
        ElectronAnalyzer(const int& era, const float& ptCut, const float& etaCut, eToken& eleToken, genPartToken& genParticleToken, const std::string& systematic="");
        ElectronAnalyzer(const int& era, const float& ptCut, const float& etaCut, TTreeReader& reader);

        void BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
