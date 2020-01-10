#ifndef JETANALYZER_H
#define JETANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>

#include <random>

#include <CondFormats/BTauObjects/interface/BTagCalibration.h>
#include <CondTools/BTau/interface/BTagCalibrationReader.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

class JetAnalyzer: public BaseAnalyzer{
    enum JetType {AK4, AK8};

    private:
        //Bool for checking if data file
        bool isData;

        //Map for SF files
        std::map<int, std::vector<std::string>> JECMC, JECDATA;
        std::map<int, std::string> JECUNC, JMESF, JMEPtReso;
        std::map<JetType, std::map<int, std::string>> bTagSF;

        //Values for bTag cuts
        std::map<JetType, std::map<int, std::vector<float>>> bTagCuts; 

        //Classes for reading btag SF
        std::map<JetType, BTagCalibration> calib;
        std::map<JetType, BTagCalibrationReader> looseReader, mediumReader, tightReader;

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
        int era;
        float ptCut, etaCut;

        //EDM Token for MINIAOD analysis
        std::vector<jToken> jetTokens;
        std::vector<genjToken> genjetTokens;
        mToken metToken;
        edm::EDGetTokenT<double> rhoToken;
        genPartToken genParticleToken;
        secvtxToken vertexToken;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> fatJetPt, fatJetEta, fatJetPhi, fatJetMass, fatJetArea, fatJetCSV;
        std::unique_ptr<TTreeReaderArray<float>> fatJetTau1, fatJetTau2, fatJetTau3;

        std::unique_ptr<TTreeReaderArray<int>> jetGenIdx;
        std::unique_ptr<TTreeReaderArray<float>> jetMass, jetPt, jetEta, jetPhi, jetArea, jetDeepBValue;

        std::unique_ptr<TTreeReaderArray<float>> genJetPt, genJetEta, genJetPhi, genJetMass;
        std::unique_ptr<TTreeReaderArray<float>> genFatJetPt, genFatJetEta, genFatJetPhi, genFatJetMass;

        std::unique_ptr<TTreeReaderValue<float>> metPhi, metPt, valueHT, jetRho;

        //Parameter for HT
        float HT, metPx, metPy;
        int runNumber;

        //Vector with output varirables of the output tree
        std::vector<std::string> JetfloatNames, FatJetfloatNames, JetParticlefloatNames, boolNames;
        std::vector<std::vector<float>> JetfloatVariables, FatJetfloatVariables, JetParticlefloatVariables, VertexfloatVariables;
        std::vector<std::vector<bool>> JetboolVariables, FatJetboolVariables;

        //Get jet energy correction
        std::map<JetType, FactorizedJetCorrector*> jetCorrector;
        void SetCorrector(const JetType &type, const int& runNumber);
        float CorrectEnergy(const TLorentzVector &jet, const float &rho, const float &area, const JetType &type);

        //Get JER smear factor
        float SmearEnergy(const TLorentzVector &jet, const float &rho, const float &coneSize, const JetType &type, const std::vector<reco::GenJet> &genJets = {});

        //Set Gen particle information
        std::map<JetType, TLorentzVector> genJet; 
        int SetGenParticles(TLorentzVector& validJet, const int &i, const int &pdgID, const JetType &type, const std::vector<reco::GenParticle>& genParticle={});

    public:
        JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader);
        JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, std::vector<jToken>& jetTokens, std::vector<genjToken>& genjetTokens, mToken &metToken, edm::EDGetTokenT<double> &rhoToken, genPartToken& genParticleToken, secvtxToken& vertexToken, const std::string& systematic = "");

        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
