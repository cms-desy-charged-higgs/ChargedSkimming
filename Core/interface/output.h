#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <vector>
#include <memory>
#include <array>

#include <TTree.h>
#include <TH1F.h>

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

const int eleMax = 10;
const int muMax = 10;
const int jetMax = 20;
const int fatJetMax = 5;
const int isotrkMax = 20;

class Output {
    public:
        //Weights
        short nTrueInt;
        float pdfWeight[102], scaleWeight[8];

        //Trigger
        std::vector<short> triggers; 

        //Electron related stuff
        std::array<float, eleMax> elePt, elePtEnergyScaleUp, 
                                  elePtEnergyScaleDown, elePtEnergySigmaUp, elePtEnergySigmaDown,
                                  eleEta, elePhi, 
                                  eleRecoSF, eleRecoSFUp, eleRecoSFDown, 
                                  eleLooseSF, eleLooseSFUp, eleLooseSFDown, 
                                  eleMediumSF, eleMediumSFUp, eleMediumSFDown, 
                                  eleTightSF, eleTightSFUp, eleTightSFDown,
                                  eleMediumMVASF, eleMediumMVASFUp, eleMediumMVASFDown,
                                  eleTightMVASF, eleTightMVASFUp, eleTightMVASFDown,
                                  eleDxy, eleDz, eleRelJetIso,
                                  eleIso03, eleMiniIso,
                                  eleGenPt, eleGenEta, eleGenPhi;

        std::array<short, eleMax> eleCutID, eleMVAID, eleCharge, 
                                  eleGenID, eleGenMotherID, eleGenGrandMotherID;

        short nElectrons;

        //Muon related stuff
        std::array<float, muMax>  muPt, muPtUp, muPtDown, 
                                  muEta, muPhi,
                                  muLooseIsoSF, muLooseIsoSFUp, muLooseIsoSFDown,
                                  muTightIsoSF, muTightIsoSFUp, muTightIsoSFDown,
                                  muLooseSF, muLooseSFUp, muLooseSFDown,
                                  muMediumSF, muMediumSFUp, muMediumSFDown,
                                  muTightSF, muTightSFUp, muTightSFDown,
                                  muTriggerSF, muTriggerSFUp, muTriggerSFDown,
                                  muDxy, muDz, muRelJetIso,
                                  muIso03, muIso04, muMiniIso,
                                  muGenPt, muGenEta, muGenPhi;
            
        std::array<short, muMax> muCutID, muMVAID, muCharge, 
                                 muGenID, muGenMotherID, muGenGrandMotherID;

        short nMuons;

        //Jet related stuff
        std::array<float, jetMax> jetPt, jetPtJMEUp, jetPtJMEDown,
                                  jetMass, jetMassJMEUp, jetMassJMEDown,
                                  jetEta, jetPhi,
                                  jetLooseDeepCSVSF, jetMediumDeepCSVSF, jetTightDeepCSVSF, 
                                  jetLooseDeepJetSF, jetMediumDeepJetSF, jetTightDeepJetSF, 
                                  jetDeepJet, jetDeepCSV,
                                  jetJEC, jetJME,
                                  jetGenPt, jetGenEta, jetGenPhi;
              
        std::vector<std::array<float, jetMax>>  jetPtJECUp, jetPtJECDown, 
                                                jetMassJECUp, jetMassJECDown,
                                                jetLooseDeepCSVSFDown, jetLooseDeepCSVSFLightDown, 
                                                jetLooseDeepCSVSFUp, jetLooseDeepCSVSFLightUp,
                                                jetMediumDeepCSVSFDown, jetMediumDeepCSVSFLightDown, 
                                                jetMediumDeepCSVSFUp, jetMediumDeepCSVSFLightUp,  
                                                jetTightDeepCSVSFDown, jetTightDeepCSVSFLightDown, 
                                                jetTightDeepCSVSFUp, jetTightDeepCSVSFLightUp,
                                                jetLooseDeepJetSFDown, jetLooseDeepJetSFLightDown, 
                                                jetLooseDeepJetSFUp, jetLooseDeepJetSFLightUp,
                                                jetMediumDeepJetSFDown, jetMediumDeepJetSFLightDown, 
                                                jetMediumDeepJetSFUp, jetMediumDeepJetSFLightUp,  
                                                jetTightDeepJetSFDown, jetTightDeepJetSFLightDown, 
                                                jetTightDeepJetSFUp, jetTightDeepJetSFLightUp; 
              
        std::array<float, jetMax> subJetPt, subJetPtJMEUp, subJetPtJMEDown,
                                  subJetMass, subJetMassJMEUp, subJetMassJMEDown,
                                  subJetEta, subJetPhi,
                                  subJetLooseDeepCSVSF, subJetMediumDeepCSVSF, subJetTightDeepCSVSF,
                                  subJetLooseDeepJetSF, subJetMediumDeepJetSF, subJetTightDeepJetSF,
                                  subJetDeepJet, subJetDeepCSV,
                                  subJetJEC, subJetJME,
                                  subJetGenPt, subJetGenEta, subJetGenPhi;
            
        std::vector<std::array<float, jetMax>>  subJetPtJECUp, subJetPtJECDown, 
                                                subJetMassJECUp, subJetMassJECDown,
                                                subJetLooseDeepCSVSFDown, subJetLooseDeepCSVSFLightDown, 
                                                subJetLooseDeepCSVSFUp, subJetLooseDeepCSVSFLightUp,
                                                subJetMediumDeepCSVSFDown, subJetMediumDeepCSVSFLightDown, 
                                                subJetMediumDeepCSVSFUp, subJetMediumDeepCSVSFLightUp,  
                                                subJetTightDeepCSVSFDown, subJetTightDeepCSVSFLightDown, 
                                                subJetTightDeepCSVSFUp, subJetTightDeepCSVSFLightUp,
                                                subJetLooseDeepJetSFDown, subJetLooseDeepJetSFLightDown, 
                                                subJetLooseDeepJetSFUp, subJetLooseDeepJetSFLightUp,
                                                subJetMediumDeepJetSFDown, subJetMediumDeepJetSFLightDown, 
                                                subJetMediumDeepJetSFUp, subJetMediumDeepJetSFLightUp,  
                                                subJetTightDeepJetSFDown, subJetTightDeepJetSFLightDown, 
                                                subJetTightDeepJetSFUp, subJetTightDeepJetSFLightUp; 
              
        std::array<float, fatJetMax> fatJetPt, fatJetPtJMEUp, fatJetPtJMEDown,
                                     fatJetMass, fatJetMassJMEUp, fatJetMassJMEDown,
                                     fatJetEta, fatJetPhi,
                                     fatJetTau1, fatJetTau2, fatJetTau3,
                                     fatJetJEC, fatJetJME;
              
        std::vector<std::array<float, fatJetMax>> fatJetPtJECUp, fatJetPtJECDown, 
                                                  fatJetMassJECUp, fatJetMassJECDown;
              
        float metPt, metPtUnclusteredUp, metPtUnclusteredDown, metPtJMEUp, metPtJMEDown, 
              metPhi, metPhiUnclusteredUp, metPhiUnclusteredDown, metPhiJMEUp, metPhiJMEDown;

        std::vector<float> metPtJECUp, metPtJECDown, metPhiJECUp, metPhiJECDown;

        std::array<short, jetMax> fatJetIdx,
                                  jetID, jetPUID,
                                  jetDeepJetID, jetDeepCSVID, 
                                  subJetDeepJetID, subJetDeepCSVID,
                                  jetPartFlav, subJetPartFlav,
                                  jetGenID, jetGenMotherID, jetGenGrandMotherID,
                                  subJetGenID, subJetGenMotherID, subJetGenGrandMotherID;
              
        std::array<short, fatJetMax> fatJetDAK8ID;

        short nJets, nSubJets, nFatJets;

        //Iso track related
        std::array<float, isotrkMax> isotrkPt, isotrkEta, isotrkPhi, 
                                     isotrkDxy, isotrkDz, 
                                     isotrkIso03, isotrkMiniIso;

        std::array<short, isotrkMax> isotrkPDG, isotrkCharge;

        short isotrkSize;

        //Misc related
        short nParton;

        long evNr;

        float preFire, preFireUp, preFireDown;


        //Function to attach branches to output trees
        void RegisterTrigger(const std::vector<std::string>& triggerNames, const std::vector<std::shared_ptr<TTree>>& trees);
        void Register(const std::string& name, const std::vector<std::shared_ptr<TTree>>& trees, pt::ptree& skim, const bool& isData = false);
};

#endif
