#ifndef NANOREADER_H
#define NANOREADER_H

#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <Math/Vector4D.h>

#include <ChargedSkimming/Core/interface/input.h>

class NanoInput : public Input {
    private:
        std::shared_ptr<TFile> inputFile;
        std::shared_ptr<TTree> inputTree;

        //Weight related
        TLeaf* pdfWeightL;
        TLeaf* scaleWeightL;
        TLeaf* nTrueIntL;

        //Trigger related
        std::vector<TLeaf*> triggerL;
        std::vector<TLeaf*> METFilterL;
        
        //Electron related
        TLeaf* elePtL;
        TLeaf* eleECorrL;
        TLeaf* eleScaleUpL;
        TLeaf* eleScaleDownL;
        TLeaf* eleSigmaUpL;
        TLeaf* eleSigmaDownL;
        TLeaf* eleMassL;
        TLeaf* eleEtaL;
        TLeaf* elePhiL;
        TLeaf* eleIsoL;
        TLeaf* eleMiniIsoL;
        TLeaf* eleChargeL;
        TLeaf* eleCutIDL;
        TLeaf* eleMVAIDLooseL;
        TLeaf* eleMVAIDMediumL;
        TLeaf* eleMVAIDTightL;
        TLeaf* eleDxyL;
        TLeaf* eleDzL;
        TLeaf* eleConvVetoL;
        TLeaf* eleRelJetIsoL;

        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> eleP4;

        //Muon related
        TLeaf* muPtL;
        TLeaf* muEtaL;
        TLeaf* muPhiL;
        TLeaf* muChargeL;
        TLeaf* muMiniIsoL;
        TLeaf* muIsoL;
        TLeaf* muCutIDLooseL;
        TLeaf* muCutIDMediumL;
        TLeaf* muCutIDTightL;
        TLeaf* muMVAIDL;
        TLeaf* muDxyL;
        TLeaf* muDzL;
        TLeaf* muRelJetIsoL;

        //Jet related
        TLeaf* rhoL;

        TLeaf* metPtL;
        TLeaf* metPhiL;
        TLeaf* metDeltaUnClustXL;
        TLeaf* metDeltaUnClustYL;

        TLeaf* jetPtL;
        TLeaf* jetEtaL;
        TLeaf* jetPhiL;
        TLeaf* jetMassL;
        TLeaf* jetAreaL;
        TLeaf* jetDeepJetL;
        TLeaf* jetDeepCSVL;
        TLeaf* jetPartFlavL;
        TLeaf* jetRawFacL;

        TLeaf* genJetPtL;
        TLeaf* genJetEtaL;
        TLeaf* genJetPhiL;

        TLeaf* fatJetPtL;
        TLeaf* fatJetEtaL;
        TLeaf* fatJetPhiL;
        TLeaf* fatJetMassL;
        TLeaf* fatJetAreaL;
        TLeaf* fatJetTau1L;
        TLeaf* fatJetTau2L;
        TLeaf* fatJetTau3L;
        TLeaf* fatJetDAK8HiggsL;
        TLeaf* fatJetDAK8QCDL;
        TLeaf* fatJetDAK8TvsQCDL;
        TLeaf* fatJetDAK8ZvsQCDL;
        TLeaf* fatJetDAK8WvsQCDL;
        TLeaf* fatJetRawFacL;

        TLeaf* genFatJetPtL;
        TLeaf* genFatJetEtaL;
        TLeaf* genFatJetPhiL;

        //Iso. track related
        TLeaf* isotrkPtL;
        TLeaf* isotrkPhiL;
        TLeaf* isotrkEtaL;
        TLeaf* isotrkDxyL;
        TLeaf* isotrkDzL;
        TLeaf* isotrkPDGL;
        TLeaf* isotrkIsoL;
        TLeaf* isotrkMiniIsoL;

        //Misc related
        TLeaf* evNrL;
        TLeaf* nPartonL;

        //Gen part related
        TLeaf* genPDGL;
        TLeaf* genMotherIdxL;
        TLeaf* genPtL;
        TLeaf* genPhiL; 
        TLeaf* genEtaL;
        TLeaf* genMassL;

        std::size_t entry;

        //Helper function https://www.wolframalpha.com/input/?i=h%2F%28h%2Bt%29+%3D+s+solve+for+h
        float demangleDK8(const float& AvsB, const float& B){
            if(AvsB != 1. and B != 0){
                return - (AvsB*B)/(AvsB -1);
            }

            else{
                return AvsB;
            }
        }

    public:
        NanoInput(const std::string& fileName, const std::string& treeName);
        void SetEntry(const std::size_t& entry){this->entry = entry;}
        std::size_t GetEntries(){return inputTree->GetEntries();}

        void SetWeight();
        void GetWeightEntry();

        void SetTrigger(const std::vector<std::string>& names, const bool& isMETFilter);
        void ReadTrigger();
        void ReadMETFilter();
        void GetTrigger();
        void GetMETFilter();

        void ReadEleEntry();
        void GetElectron(const std::size_t& idx);

        void ReadMuEntry();
        void GetMuon(const std::size_t& idx);

        void ReadJetEntry(const bool& isData);
        void GetJet(const std::size_t& idx);
        void GetFatJet(const std::size_t& idx);
        void GetGenJet(const std::size_t& idx);
        void GetGenFatJet(const std::size_t& idx);

        void ReadIsotrkEntry();
        void GetIsotrk(const std::size_t& idx);

        void ReadMiscEntry(const bool& isData);
        void GetMisc();

        void ReadGenEntry();
        void GetGenPart(const std::size_t& idx);
};

#endif
