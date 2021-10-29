#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <vector>

struct Input{
    public:
        //Weight related
        float pdfWeight[102], scaleWeight[8];
        short nTrueInt;

        //Trigger related
        std::vector<short> triggers;
        std::vector<short> METFilter;

        //Electron related
        short eleSize, eleCharge, eleCutID, eleMVAID, eleConvVeto;
        float elePt, elePtSigmaUp, elePtSigmaDown, elePtScaleUp, elePtScaleDown, eleEta, elePhi, eleDxy, eleDz, eleRelJetIso, eleIso, eleMiniIso;

        //Muon related
        short muSize, muCharge, muCutID, muMVAID;
        float muPt, muEta, muPhi, muDxy, muDz, muRelJetIso, muIso, muMiniIso;

        //Jet related
        short jetSize, fatJetSize, genJetSize, genFatJetSize,
              jetGenIdx, fatJetGenIdx, fatJetDAK8ID;
        float metPt, metPhi, metDeltaUnClustX, metDeltaUnClustY,
              jetPt, jetMass, jetPtRaw, jetMassRaw,
              jetEta, jetPhi, jetArea,
              jetDeepJet, jetDeepCSV, jetPartFlav,
              genJetEta, genJetPhi, genJetPt,
              fatJetPt, fatJetMass, fatJetPtRaw, fatJetMassRaw,
              fatJetPhi, fatJetEta, fatJetArea,
              fatJetTau1, fatJetTau2, fatJetTau3,
              genFatJetEta, genFatJetPhi, genFatJetPt,
              rho;

        //Iso. track related
        short isotrkSize, isotrkPDG;
        float isotrkPt, isotrkEta, isotrkPhi, isotrkDxy, isotrkDz, isotrkIso, isotrkMiniIso;

        //Misc related
        short nParton;
        long evNr;

        //Gen part related
        short genSize, genPDG, genMotherIdx;
        float genPt, genPhi, genEta, genMass;
};

#endif
