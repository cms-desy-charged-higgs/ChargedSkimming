#ifndef JETANALYZER_H
#define JETANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/btagcsvreader.h>
#include <ChargedSkimming/Skimming/interface/util.h>

#include <cmath>
#include <random>
#include <experimental/filesystem>

#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/JetReco/interface/GenJet.h>

#include <TH2F.h>

class JetAnalyzer: public BaseAnalyzer{
    enum JetType {SUBAK4, AK4, AK8, PF, VTX};

    private:
        //Bool for checking if data file
        bool isData;

        //Cut values
        float ptCut, etaCut;
        float CSVLoose, CSVMedium, CSVTight, DeepLoose, DeepMedium, DeepTight;

        //Map for SF files
        std::vector<std::string> JECMC, JECDATA;

        //Classes for reading btag SF
        BTagCSVReader CSVReader, DeepReader;

        std::vector<TH2F*> bTagEffBLoose, bTagEffBMedium, bTagEffBTight, bTagEffCLoose, bTagEffCMedium, bTagEffCTight, bTagEffLightLoose, bTagEffLightMedium, bTagEffLightTight;
        TH2F* bTotal;
        TH2F* cTotal; 
        TH2F* lightTotal;

        //Classes for reading jet energy SF 
        JME::JetParameters jetParameter;
        std::map<JetType, JME::JetResolution> resolution;
        std::map<JetType, JME::JetResolutionScaleFactor> resolution_sf;

        //JEC/JER systematics
        std::map<JetType, JetCorrectionUncertainty*> jecUnc;
        std::string jecSyst;
        bool isUp;
        bool isJERsyst = false;

        //Input for selecting jets
        std::string era;
        std::shared_ptr<Token> tokens;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> fatJetPt, fatJetEta, fatJetPhi, fatJetMass, fatJetArea, fatJetCSV;
        std::unique_ptr<TTreeReaderArray<float>> fatJetTau1, fatJetTau2, fatJetTau3;

        std::unique_ptr<TTreeReaderArray<int>> jetGenIdx, jetFlavour;
        std::unique_ptr<TTreeReaderArray<float>> jetMass, jetPt, jetEta, jetPhi, jetArea, jetDeepBValue;

        std::unique_ptr<TTreeReaderArray<float>> genJetPt, genJetEta, genJetPhi, genJetMass;
        std::unique_ptr<TTreeReaderArray<float>> genFatJetPt, genFatJetEta, genFatJetPhi, genFatJetMass;

        std::unique_ptr<TTreeReaderValue<float>> metPhi, metPt, jetRho;

        //Parameter for HT
        float metPT, metPHI, metPTUp, metPTDown, metPHIUp, metPHIDown;
        int runNumber;

        //Vector with output varirables of the output tree
        std::map<std::pair<std::string, std::string>, float*> floatVar;
        std::map<std::pair<std::string, std::string>, short*> intVar;

        std::map<JetType, float[300]> Pt, Eta, Phi, Mass, Vx, Vy, Vz, DeepAK8Higgs, DeepAK8Top, DeepAK8W, DeepAK8DY, DeepAK8QCD, Njettiness1, Njettiness2, Njettiness3, looseCSVbTagSF, looseCSVbTagSFUp, looseCSVbTagSFDown, mediumCSVbTagSF, mediumCSVbTagSFUp, mediumCSVbTagSFDown, tightCSVbTagSF, tightCSVbTagSFUp, tightCSVbTagSFDown, looseDeepbTagSF, looseDeepbTagSFUp, looseDeepbTagSFDown, mediumDeepbTagSF, mediumDeepbTagSFUp, mediumDeepbTagSFDown, tightDeepbTagSF, tightDeepbTagSFUp, tightDeepbTagSFDown, DeepScore, CSVScore, corrJEC, corrJER;

        std::map<JetType, short[300]> TrueFlavour, Charge, FatJetIdx, partID, mothID, grandID, bTagCSVID, bTagDeepID, DeepAK8Class;

        short nJets, nSubJets, nFatJets, nVtx, nPFcands;

        //Get jet energy correction
        std::map<JetType, FactorizedJetCorrector*> jetCorrector;
        void SetCorrector(const JetType &type, const int& runNumber);
        float CorrectEnergy(const float& pt, const float& eta, const float& rho, const float& area, const JetType& type);

        //Get JER smear factor
        float SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, const float& coneSize, const JetType& type, const std::vector<reco::GenJet> &genJets = {});

        //Set Gen particle information
        std::map<JetType, ROOT::Math::PtEtaPhiMVector> genJet; 
        std::vector<int> alreadySeenNANO;
        std::vector<const reco::Candidate*> alreadySeenMINI;

        void SetGenParticles(const int& currentSize, const int& i, const float& pt, const float& eta, const float& phi, const std::vector<int>& pdgID, const JetType &type, const std::vector<reco::GenParticle>& genParticle = {});

    public:
        JetAnalyzer(const std::string& era, TTreeReader& reader);
        JetAnalyzer(const std::string& era, const std::shared_ptr<Token>& tokens, const std::string& systematic = "");

        void BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf);
        void Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
