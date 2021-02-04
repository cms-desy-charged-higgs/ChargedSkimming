#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <DataFormats/PatCandidates/interface/Electron.h>

typedef edm::EDGetTokenT<std::vector<pat::Electron>> eToken;

class ElectronAnalyzer: public BaseAnalyzer {
    private:
        //Check if data
        bool isData;

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
        std::map<std::string, float*> floatVar;
        std::map<std::string, short*> intVar;

        float Pt[20], Eta[20], Phi[20], recoSF[20], recoSFUp[20], recoSFDown[20], looseSF[20], looseSFUp[20], looseSFDown[20], mediumSF[20], mediumSFUp[20], mediumSFDown[20], tightSF[20], tightSFUp[20], tightSFDown[20], Isolation[20];

        short ID[20], Charge[20], partID[20], mothID[20], grandID[20];

        short nElectrons;

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
