#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <vector>
#include <memory>

#include <TTree.h>

const int eleMax = 10;
const int muMax = 10;
const int jetMax = 20;
const int fatJetMax = 5;
const int isotrkMax = 20;

class Output{
    public:
        //Weights
        short nTrueInt;
        float pdfWeight[102], scaleWeight[8];

        //Trigger
        std::vector<short> triggers; 

        //Electron related stuff
        float elePt[eleMax], eleEta[eleMax], elePhi[eleMax], 
              eleRecoSF[eleMax], eleRecoSFUp[eleMax], eleRecoSFDown[eleMax], 
              eleLooseSF[eleMax], eleLooseSFUp[eleMax], eleLooseSFDown[eleMax], 
              eleMediumSF[eleMax], eleMediumSFUp[eleMax], eleMediumSFDown[eleMax], 
              eleTightSF[eleMax], eleTightSFUp[eleMax], eleTightSFDown[eleMax],
              eleMediumMVASF[eleMax], eleMediumMVASFUp[eleMax], eleMediumMVASFDown[eleMax],
              eleTightMVASF[eleMax], eleTightMVASFUp[eleMax], eleTightMVASFDown[eleMax],
              eleDxy[eleMax], eleDz[eleMax], eleRelJetIso[eleMax],
              eleIso[eleMax], eleMiniIso[eleMax],
              eleGenPt[eleMax], eleGenEta[eleMax], eleGenPhi[eleMax];

        short eleCutID[eleMax], eleMVAID[eleMax], eleCharge[eleMax], 
              eleGenID[eleMax], eleGenMotherID[eleMax], eleGenGrandMotherID[eleMax];

        short nElectrons;

        //Muon related stuff
        float muPt[muMax], muEta[muMax], muPhi[muMax],
              muLooseIsoSF[muMax], muLooseIsoSFUp[muMax], muLooseIsoSFDown[muMax],
              muTightIsoSF[muMax], muTightIsoSFUp[muMax], muTightIsoSFDown[muMax],
              muLooseSF[muMax], muLooseSFUp[muMax], muLooseSFDown[muMax],
              muMediumSF[muMax], muMediumSFUp[muMax], muMediumSFDown[muMax],
              muTightSF[muMax], muTightSFUp[muMax], muTightSFDown[muMax],
              muTriggerSF[muMax], muTriggerSFUp[muMax], muTriggerSFDown[muMax],
              muDxy[muMax], muDz[muMax], muRelJetIso[muMax],
              muIso[muMax], muMiniIso[muMax],
              muGenPt[muMax], muGenEta[muMax], muGenPhi[muMax];

        short muCutID[muMax], muMVAID[muMax], muCharge[muMax], 
              muGenID[muMax], muGenMotherID[muMax], muGenGrandMotherID[muMax];

        short nMuons;

        //Jet related stuff
        float jetPt[jetMax], jetEta[jetMax], jetPhi[jetMax], jetMass[jetMax],
              jetLooseDeepCSVSF[jetMax], jetLooseDeepCSVSFUp[jetMax], jetLooseDeepCSVSFDown[jetMax],
              jetMediumDeepCSVSF[jetMax], jetMediumDeepCSVSFUp[jetMax], jetMediumDeepCSVSFDown[jetMax],
              jetTightDeepCSVSF[jetMax], jetTightDeepCSVSFUp[jetMax], jetTightDeepCSVSFDown[jetMax],
              jetLooseDeepJetSF[jetMax], jetLooseDeepJetSFUp[jetMax], jetLooseDeepJetSFDown[jetMax],
              jetMediumDeepJetSF[jetMax], jetMediumDeepJetSFUp[jetMax], jetMediumDeepJetSFDown[jetMax],
              jetTightDeepJetSF[jetMax], jetTightDeepJetSFUp[jetMax], jetTightDeepJetSFDown[jetMax],
              jetDeepJet[jetMax], jetDeepCSV[jetMax],
              jetJEC[jetMax], jetJME[jetMax],
              jetGenPt[jetMax], jetGenEta[jetMax], jetGenPhi[jetMax];
        float subJetPt[jetMax], subJetEta[jetMax], subJetPhi[jetMax], subJetMass[jetMax],
              subJetLooseDeepCSVSF[jetMax], subJetLooseDeepCSVSFUp[jetMax], subJetLooseDeepCSVSFDown[jetMax],
              subJetMediumDeepCSVSF[jetMax], subJetMediumDeepCSVSFUp[jetMax], subJetMediumDeepCSVSFDown[jetMax],
              subJetTightDeepCSVSF[jetMax], subJetTightDeepCSVSFUp[jetMax], subJetTightDeepCSVSFDown[jetMax],
              subJetLooseDeepJetSF[jetMax], subJetLooseDeepJetSFUp[jetMax], subJetLooseDeepJetSFDown[jetMax],
              subJetMediumDeepJetSF[jetMax], subJetMediumDeepJetSFUp[jetMax], subJetMediumDeepJetSFDown[jetMax],
              subJetTightDeepJetSF[jetMax], subJetTightDeepJetSFUp[jetMax], subJetTightDeepJetSFDown[jetMax],
              subJetDeepJet[jetMax], subJetDeepCSV[jetMax],
              subJetJEC[jetMax], subJetJME[jetMax],
              subJetGenPt[jetMax], subJetGenEta[jetMax], subJetGenPhi[jetMax];
        float fatJetPt[fatJetMax], fatJetEta[fatJetMax], fatJetPhi[fatJetMax], fatJetMass[fatJetMax], 
              fatJetTau1[fatJetMax], fatJetTau2[fatJetMax], fatJetTau3[fatJetMax],
              fatJetJEC[fatJetMax], fatJetJME[fatJetMax];
        float metPt, metPtUp, metPtDown, metPhi, metPhiUp, metPhiDown;

        short fatJetIdx[jetMax], 
              jetDeepJetID[jetMax], jetDeepCSVID[jetMax], 
              subJetDeepJetID[jetMax], subJetDeepCSVID[jetMax],
              jetPartFlav[jetMax], subJetPartFlav[jetMax],
              fatJetDAK8ID[fatJetMax],
              jetGenID[jetMax], jetGenMotherID[jetMax], jetGenGrandMotherID[jetMax],
              subJetGenID[jetMax], subJetGenMotherID[jetMax], subJetGenGrandMotherID[jetMax];

        short nJets, nSubJets, nFatJets;

        //Iso track related
        float isotrkPt[isotrkMax], isotrkEta[isotrkMax], isotrkPhi[isotrkMax], 
              isotrkDxy[isotrkMax], isotrkDz[isotrkMax], isotrkIso[isotrkMax], isotrkMiniIso[isotrkMax];

        short isotrkPDG[isotrkMax], isotrkCharge[isotrkMax];

        short isotrkSize;

        //Misc related
        short nParton;
        long evNr;

        //Function to attach branches to output trees
        void RegisterTrigger(const std::vector<std::string>& triggerNames, const std::vector<std::shared_ptr<TTree>>& trees);
        void Register(const std::string& name, const std::vector<std::shared_ptr<TTree>>& trees, const bool& isData = false, const bool& isSyst = false);
};

#endif
