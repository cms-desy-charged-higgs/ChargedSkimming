#ifndef ELECTRONANALYZER_H
#define ELECTRONANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/tokens.h>

#include <DataFormats/PatCandidates/interface/Electron.h>

class ElectronAnalyzer: public BaseAnalyzer {
    private:
        //Check if data
        bool isData;

        //Input
        std::string era;
        std::shared_ptr<Token> tokens;

        //Kinematic cut criteria
        float ptCut, etaCut;

        //Hist with scale factors
        std::shared_ptr<TH2F> looseSFhist, mediumSFhist, tightSFhist, recoSFhist;

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
        ElectronAnalyzer(const std::string& era, TTreeReader& reader);
        ElectronAnalyzer(const std::string& era, const std::shared_ptr<Token>& tokens, std::string& systematic);

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
