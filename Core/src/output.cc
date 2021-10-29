#include <ChargedSkimming/Core/interface/output.h>
#include <iostream>

void Output::RegisterTrigger(const std::vector<std::string>& triggerNames, const std::vector<std::shared_ptr<TTree>>& trees){
    triggers = std::vector<short>(triggerNames.size(), 1);

    for(std::size_t i = 0; i < triggerNames.size(); ++i){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch(triggerNames[i].c_str(), &triggers[i], (triggerNames[i] + "/S").c_str());
        }
    }
}

void Output::Register(const std::string& name, const std::vector<std::shared_ptr<TTree>>& trees, const bool& isData, const bool& isSyst){
    if(name == "Weight"){
        for(const std::shared_ptr<TTree>& tree: trees){
            if(!isData){
                tree->Branch("Weight_nTrueInt", &nTrueInt, "Weight_nTrueInt/S");

                if(!isSyst){
                    tree->Branch("Weight_pdf", pdfWeight, "Weight_pdf[102]/F");
                    tree->Branch("Weight_scale", scaleWeight, "Weight_scale[8]/F");
                }
            }
        }
    }

    if(name == "Electron"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Electron_Size", &nElectrons, "Electron_Size/S");

            tree->Branch("Electron_Pt", elePt, "Electron_Pt[Electron_Size]/F");
            tree->Branch("Electron_Eta", eleEta, "Electron_Eta[Electron_Size]/F");
            tree->Branch("Electron_Phi", elePhi, "Electron_Phi[Electron_Size]/F");
            tree->Branch("Electron_Isolation", eleIso, "Electron_Isolation[Electron_Size]/F");
            tree->Branch("Electron_MiniIsolation", eleMiniIso, "Electron_MiniIsolation[Electron_Size]/F");
            tree->Branch("Electron_Dxy", eleDxy, "Electron_Dxy[Electron_Size]/F");
            tree->Branch("Electron_Dz", eleDz, "Electron_Dz[Electron_Size]/F");
            tree->Branch("Electron_JetRelIsolation", eleRelJetIso, "Electron_JetRelIsolation[Electron_Size]/F");

            tree->Branch("Electron_Charge", eleCharge, "Electron_Charge[Electron_Size]/S");
            tree->Branch("Electron_CutID", eleCutID, "Electron_CutID[Electron_Size]/S");
            tree->Branch("Electron_MVAID", eleMVAID, "Electron_MVAID[Electron_Size]/S");

            if(!isData){
                tree->Branch("Electron_GenPt", eleGenPt, "Electron_GenPt[Electron_Size]/F");
                tree->Branch("Electron_GenEta", eleGenEta, "Electron_GenEta[Electron_Size]/F");
                tree->Branch("Electron_GenPhi", eleGenPhi, "Electron_GenPhi[Electron_Size]/F");
                tree->Branch("Electron_GenID", eleGenID, "Electron_GenID[Electron_Size]/S");
                tree->Branch("Electron_GenMotherID", eleGenMotherID, "Electron_GenMotherID[Electron_Size]/S");
                tree->Branch("Electron_GenGrandMotherID", eleGenGrandMotherID, "Electron_GenGrandMotherID[Electron_Size]/S");

                tree->Branch("Electron_RecoSF", eleRecoSF, "Electron_RecoSF[Electron_Size]/F");
                tree->Branch("Electron_looseSF", eleLooseSF, "Electron_looseSF[Electron_Size]/F");
                tree->Branch("Electron_mediumSF", eleMediumSF, "Electron_mediumSF[Electron_Size]/F");
                tree->Branch("Electron_tightSF", eleTightSF, "Electron_tightSF[Electron_Size]/F");
                tree->Branch("Electron_mediumMVASF", eleMediumMVASF, "Electron_mediumMVASF[Electron_Size]/F");
                tree->Branch("Electron_tightMVASF", eleTightMVASF, "Electron_tightMVASF[Electron_Size]/F");

                if(!isSyst){
                    tree->Branch("Electron_RecoSFUp", eleRecoSFUp, "Electron_RecoSFUp[Electron_Size]/F");
                    tree->Branch("Electron_looseSFUp", eleLooseSFUp, "Electron_looseSFUp[Electron_Size]/F");
                    tree->Branch("Electron_mediumSFUp", eleMediumSFUp, "Electron_mediumSFUp[Electron_Size]/F");
                    tree->Branch("Electron_tightSFUp", eleTightSFUp, "Electron_tightSFUp[Electron_Size]/F");
                    tree->Branch("Electron_mediumMVASFUp", eleMediumMVASFUp, "Electron_mediumMVASFUp[Electron_Size]/F");
                    tree->Branch("Electron_tightMVASFUp", eleTightMVASFUp, "Electron_tightMVASFUp[Electron_Size]/F");
                    tree->Branch("Electron_RecoSFDown", eleRecoSFDown, "Electron_RecoSFDown[Electron_Size]/F");
                    tree->Branch("Electron_looseSFDown", eleLooseSFDown, "Electron_looseSFDown[Electron_Size]/F");
                    tree->Branch("Electron_mediumSFDown", eleMediumSFDown, "Electron_mediumSFDown[Electron_Size]/F");
                    tree->Branch("Electron_tightSFDown", eleTightSFDown, "Electron_tightSFDown[Electron_Size]/F");
                    tree->Branch("Electron_mediumMVASFDown", eleMediumMVASFDown, "Electron_mediumMVASFDown[Electron_Size]/F");
                    tree->Branch("Electron_tightMVASFDown", eleTightMVASFDown, "Electron_tightMVASFDown[Electron_Size]/F");
                }
            }
        }
    }

    if(name == "Muon"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Muon_Size", &nMuons, "Muon_Size/S");

            tree->Branch("Muon_Pt", muPt, "Muon_Pt[Muon_Size]/F");
            tree->Branch("Muon_Eta", muEta, "Muon_Eta[Muon_Size]/F");
            tree->Branch("Muon_Phi", muPhi, "Muon_Phi[Muon_Size]/F");
            tree->Branch("Muon_Isolation", muIso, "Muon_Isolation[Muon_Size]/F");
            tree->Branch("Muon_MiniIsolation", muMiniIso, "Muon_MiniIsolation[Muon_Size]/F");
            tree->Branch("Muon_Dxy", muDxy, "Muon_Dxy[Muon_Size]/F");
            tree->Branch("Muon_Dz", muDz, "Muon_Dz[Muon_Size]/F");
            tree->Branch("Muon_JetRelIsolation", muRelJetIso, "Muon_JetRelIsolation[Muon_Size]/F");

            tree->Branch("Muon_Charge", muCharge, "Muon_Charge[Muon_Size]/S");
            tree->Branch("Muon_CutID", muCutID, "Muon_CutID[Muon_Size]/S");
            tree->Branch("Muon_MVAID", muMVAID, "Muon_MVAID[Muon_Size]/S");

            if(!isData){
                tree->Branch("Muon_GenPt", muGenPt, "Muon_GenPt[Muon_Size]/F");
                tree->Branch("Muon_GenEta", muGenEta, "Muon_GenEta[Muon_Size]/F");
                tree->Branch("Muon_GenPhi", muGenPhi, "Muon_GenPhi[Muon_Size]/F");
                tree->Branch("Muon_GenID", muGenID, "Muon_GenID[Muon_Size]/S");
                tree->Branch("Muon_GenMotherID", muGenMotherID, "Muon_GenMotherID[Muon_Size]/S");
                tree->Branch("Muon_GenGrandMotherID", muGenGrandMotherID, "Muon_GenGrandMotherID[Muon_Size]/S");

                tree->Branch("Muon_TriggerSF", muTriggerSF, "Muon_TriggerSF[Muon_Size]/F");
                tree->Branch("Muon_looseIsoSF", muLooseIsoSF, "Muon_looseIsoSF[Muon_Size]/F");
                tree->Branch("Muon_tightIsoSF", muTightIsoSF, "Muon_tightIsoSF[Muon_Size]/F");
                tree->Branch("Muon_looseSF", muLooseSF, "Muon_looseSF[Muon_Size]/F");
                tree->Branch("Muon_mediumSF", muMediumSF, "Muon_mediumSF[Muon_Size]/F");
                tree->Branch("Muon_tightSF", muTightSF, "Muon_tightSF[Muon_Size]/F");

                if(!isSyst){
                    tree->Branch("Muon_TriggerSFUp", muTriggerSFUp, "Muon_TriggerSFUp[Muon_Size]/F");
                    tree->Branch("Muon_looseIsoSFUp", muLooseIsoSFUp, "Muon_looseIsoSFUp[Muon_Size]/F");
                    tree->Branch("Muon_tightIsoSFUp", muTightIsoSFUp, "Muon_tightIsoSFUp[Muon_Size]/F");
                    tree->Branch("Muon_looseSFUp", muLooseSFUp, "Muon_looseSFUp[Muon_Size]/F");
                    tree->Branch("Muon_mediumSFUp", muMediumSFUp, "Muon_mediumSFUp[Muon_Size]/F");
                    tree->Branch("Muon_tightSFUp", muTightSFUp, "Muon_tightSFUp[Muon_Size]/F");
                    tree->Branch("Muon_TriggerSFDown", muTriggerSFDown, "Muon_TriggerSFDown[Muon_Size]/F");
                    tree->Branch("Muon_looseIsoSFDown", muLooseIsoSFDown, "Muon_looseIsoSFDown[Muon_Size]/F");
                    tree->Branch("Muon_tightIsoSFDown", muTightIsoSFDown, "Muon_tightIsoSFDown[Muon_Size]/F");
                    tree->Branch("Muon_looseSFDown", muLooseSFDown, "Muon_looseSFDown[Muon_Size]/F");
                    tree->Branch("Muon_mediumSFDown", muMediumSFDown, "Muon_mediumSFDown[Muon_Size]/F");
                    tree->Branch("Muon_tightSFDown", muTightSFDown, "Muon_tightSFDown[Muon_Size]/F");
                }
            }
        }
    }

    if(name == "Jet"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Jet_Size", &nJets, "Jet_Size/S");
            tree->Branch("SubJet_Size", &nSubJets, "SubJet_Size/S");
            tree->Branch("FatJet_Size", &nFatJets, "FatJet_Size/S");

            tree->Branch("MET_Pt", &metPt);
            tree->Branch("MET_PtUp", &metPtUp);
            tree->Branch("MET_PtDown", &metPtDown);
            tree->Branch("MET_Phi", &metPhi);
            tree->Branch("MET_PhiUp", &metPhiUp);
            tree->Branch("MET_PhiDown", &metPhiDown);

            tree->Branch("Jet_Pt", jetPt, "Jet_Pt[Jet_Size]/F");
            tree->Branch("Jet_Eta", jetEta, "Jet_Eta[Jet_Size]/F");
            tree->Branch("Jet_Phi", jetPhi, "Jet_Phi[Jet_Size]/F");
            tree->Branch("Jet_Mass", jetMass, "Jet_Mass[Jet_Size]/F");
            tree->Branch("Jet_DeepJet", jetDeepJet, "Jet_DeepJet[Jet_Size]/F");
            tree->Branch("Jet_DeepCSV", jetDeepCSV, "Jet_DeepCSV[Jet_Size]/F");
            tree->Branch("Jet_DeepJetID", jetDeepJetID, "Jet_DeepJetID[Jet_Size]/S");
            tree->Branch("Jet_DeepCSVID", jetDeepCSVID, "Jet_DeepCSVID[Jet_Size]/S");
            tree->Branch("Jet_JECFac", jetJEC, "Jet_JECFac[Jet_Size]/F");
            tree->Branch("Jet_JMEFac", jetJME, "Jet_JMEFac[Jet_Size]/F");

            tree->Branch("SubJet_Pt", subJetPt, "SubJet_Pt[SubJet_Size]/F");
            tree->Branch("SubJet_Eta", subJetEta, "SubJet_Eta[SubJet_Size]/F");
            tree->Branch("SubJet_Phi", subJetPhi, "SubJet_Phi[SubJet_Size]/F");
            tree->Branch("SubJet_Mass", subJetMass, "SubJet_Mass[SubJet_Size]/F");
            tree->Branch("SubJet_DeepJet", subJetDeepJet, "SubJet_DeepJet[SubJet_Size]/F");
            tree->Branch("SubJet_DeepCSV", subJetDeepCSV, "SubJet_DeepCSV[SubJet_Size]/F");
            tree->Branch("SubJet_DeepJetID", subJetDeepJetID, "SubJet_DeepJetID[SubJet_Size]/S");
            tree->Branch("SubJet_DeepCSVID", subJetDeepCSVID, "SubJet_DeepCSVID[SubJet_Size]/S");
            tree->Branch("SubJet_PartonFlavour", subJetPartFlav, "SubJet_PartonFlavour[SubJet_Size]/S");
            tree->Branch("SubJet_JECFac", subJetJEC, "SubJet_JECFac[SubJet_Size]/F");
            tree->Branch("SubJet_JMEFac", subJetJME, "SubJet_JMEFac[SubJet_Size]/F");
            tree->Branch("SubJet_FatJetIdx", fatJetIdx, "SubJet_FatJetIdx[SubJet_Size]/S");

            tree->Branch("FatJet_Pt", fatJetPt, "FatJet_Pt[FatJet_Size]/F");
            tree->Branch("FatJet_Eta", fatJetEta, "FatJet_Eta[FatJet_Size]/F");
            tree->Branch("FatJet_Phi", fatJetPhi, "FatJet_Phi[FatJet_Size]/F");
            tree->Branch("FatJet_Mass", fatJetMass, "FatJet_Mass[FatJet_Size]/F");
            tree->Branch("FatJet_Tau1", fatJetTau1, "FatJet_Tau1[FatJet_Size]/F");
            tree->Branch("FatJet_Tau2", fatJetTau2, "FatJet_Tau2[FatJet_Size]/F");
            tree->Branch("FatJet_Tau3", fatJetTau3, "FatJet_Tau3[FatJet_Size]/F");
            tree->Branch("FatJet_DeepAK8ID",fatJetDAK8ID, "FatJet_DeepAK8ID[FatJet_Size]/S");
            tree->Branch("FatJet_JECFac", fatJetJEC, "FatJet_JECFac[FatJet_Size]/F");
            tree->Branch("FatJet_JMEFac", fatJetJME, "FatJet_JMEFac[FatJet_Size]/F");

            if(!isData){
                tree->Branch("Jet_PartonFlavour", jetPartFlav, "Jet_PartonFlavour[Jet_Size]/S");
                tree->Branch("Jet_GenPt", jetGenPt, "Jet_GenPt[Jet_Size]/F");
                tree->Branch("Jet_GenEta", jetGenEta, "Jet_GenEta[Jet_Size]/F");
                tree->Branch("Jet_GenPhi", jetGenPhi, "Jet_GenPhi[Jet_Size]/F");
                tree->Branch("Jet_GenID", jetGenID, "Jet_GenID[Jet_Size]/S");
                tree->Branch("Jet_GenMotherID", jetGenMotherID, "Jet_GenMotherID[Jet_Size]/S");
                tree->Branch("Jet_GenGrandMotherID", jetGenGrandMotherID, "Jet_GenGrandMotherID[Jet_Size]/S");

                tree->Branch("SubJet_GenPt", subJetGenPt, "SubJet_GenPt[SubJet_Size]/F");
                tree->Branch("SubJet_GenEta", subJetGenEta, "SubJet_GenEta[SubJet_Size]/F");
                tree->Branch("SubJet_GenPhi", subJetGenPhi, "SubJet_GenPhi[SubJet_Size]/F");
                tree->Branch("SubJet_GenID", subJetGenID, "SubJet_GenID[SubJet_Size]/S");
                tree->Branch("SubJet_GenMotherID", subJetGenMotherID, "SubJet_GenMotherID[SubJet_Size]/S");
                tree->Branch("SubJet_GenGrandMotherID", subJetGenGrandMotherID, "SubJet_GenGrandMotherID[SubJet_Size]/S");

                tree->Branch("Jet_looseDeepCSVSF", jetLooseDeepCSVSF, "Jet_looseDeepCSVSF[Jet_Size]/F");
                tree->Branch("Jet_mediumDeepCSVSF", jetMediumDeepCSVSF, "Jet_mediumDeepCSVSF[Jet_Size]/F");
                tree->Branch("Jet_tightDeepCSVSF", jetTightDeepCSVSF, "Jet_tightDeepCSVSF[Jet_Size]/F");
                tree->Branch("Jet_looseDeepJetSF", jetLooseDeepJetSF, "Jet_looseDeepJetSF[Jet_Size]/F");
                tree->Branch("Jet_mediumDeepJetSF", jetMediumDeepJetSF, "Jet_mediumDeepJetSF[Jet_Size]/F");
                tree->Branch("Jet_tightDeepJetSF", jetTightDeepJetSF, "Jet_tightDeepJetSF[Jet_Size]/F");

                tree->Branch("SubJet_looseDeepCSVSF", subJetLooseDeepCSVSF, "SubJet_looseDeepCSVSF[Jet_Size]/F");
                tree->Branch("SubJet_mediumDeepCSVSF", subJetMediumDeepCSVSF, "SubJet_mediumDeepCSVSF[Jet_Size]/F");
                tree->Branch("SubJet_tightDeepCSVSF", subJetTightDeepCSVSF, "SubJet_tightDeepCSVSF[Jet_Size]/F");
                tree->Branch("SubJet_looseDeepJetSF", subJetLooseDeepJetSF, "SubJet_looseDeepJetSF[Jet_Size]/F");
                tree->Branch("SubJet_mediumDeepJetSF", subJetMediumDeepJetSF, "SubJet_mediumDeepJetSF[Jet_Size]/F");
                tree->Branch("SubJet_tightDeepJetSF", subJetTightDeepJetSF, "SubJet_tightDeepJetSF[Jet_Size]/F");

                if(!isSyst){
                    tree->Branch("Jet_looseDeepCSVSFUp", jetLooseDeepCSVSFUp, "Jet_looseDeepCSVSFUp[Jet_Size]/F");
                    tree->Branch("Jet_looseDeepCSVSFDown", jetLooseDeepCSVSFDown, "Jet_looseDeepCSVSFDown[Jet_Size]/F");
                    tree->Branch("Jet_mediumDeepCSVSFUp", jetMediumDeepCSVSFUp, "Jet_mediumDeepCSVSFUp[Jet_Size]/F");
                    tree->Branch("Jet_mediumDeepCSVSFDown", jetMediumDeepCSVSFDown, "Jet_mediumDeepCSVSFDown[Jet_Size]/F");
                    tree->Branch("Jet_tightDeepCSVSFUp", jetTightDeepCSVSFUp, "Jet_tightDeepCSVSFUp[Jet_Size]/F");
                    tree->Branch("Jet_tightDeepCSVSFDown", jetTightDeepCSVSFDown, "Jet_tightDeepCSVSFDown[Jet_Size]/F");
                    tree->Branch("Jet_looseDeepJetSFUp", jetLooseDeepJetSFUp, "Jet_looseDeepJetSFUp[Jet_Size]/F");
                    tree->Branch("Jet_looseDeepJetSFDown", jetLooseDeepJetSFDown, "Jet_looseDeepJetSFDown[Jet_Size]/F");
                    tree->Branch("Jet_mediumDeepJetSFUp", jetMediumDeepJetSFUp, "Jet_mediumDeepJetSFUp[Jet_Size]/F");
                    tree->Branch("Jet_mediumDeepJetSFDown", jetMediumDeepJetSFDown, "Jet_mediumDeepJetSFDown[Jet_Size]/F");
                    tree->Branch("Jet_tightDeepJetSFUp", jetTightDeepJetSFUp, "Jet_tightDeepJetSFUp[Jet_Size]/F");
                    tree->Branch("Jet_tightDeepJetSFDown", jetTightDeepJetSFDown, "Jet_tightDeepJetSFDown[Jet_Size]/F");

                    tree->Branch("SubJet_looseDeepCSVSFUp", subJetLooseDeepCSVSFUp, "SubJet_looseDeepCSVSFUp[SubJet_Size]/F");
                    tree->Branch("SubJet_looseDeepCSVSFDown", subJetLooseDeepCSVSFDown, "SubJet_looseDeepCSVSFDown[SubJet_Size]/F");
                    tree->Branch("SubJet_mediumDeepCSVSFUp", subJetMediumDeepCSVSFUp, "SubJet_mediumDeepCSVSFUp[SubJet_Size]/F");
                    tree->Branch("SubJet_mediumDeepCSVSFDown", subJetMediumDeepCSVSFDown, "SubJet_mediumDeepCSVSFDown[SubJet_Size]/F");
                    tree->Branch("SubJet_tightDeepCSVSFUp", subJetTightDeepCSVSFUp, "SubJet_tightDeepCSVSFUp[SubJet_Size]/F");
                    tree->Branch("SubJet_tightDeepCSVSFDown", subJetTightDeepCSVSFDown, "SubJet_tightDeepCSVSFDown[SubJet_Size]/F");
                    tree->Branch("SubJet_looseDeepJetSFUp", subJetLooseDeepJetSFUp, "SubJet_looseDeepJetSFUp[SubJet_Size]/F");
                    tree->Branch("SubJet_looseDeepJetSFDown", subJetLooseDeepJetSFDown, "SubJet_looseDeepJetSFDown[SubJet_Size]/F");
                    tree->Branch("SubJet_mediumDeepJetSFUp", subJetMediumDeepJetSFUp, "SubJet_mediumDeepJetSFUp[SubJet_Size]/F");
                    tree->Branch("SubJet_mediumDeepJetSFDown", subJetMediumDeepJetSFDown, "SubJet_mediumDeepJetSFDown[SubJet_Size]/F");
                    tree->Branch("SubJet_tightDeepJetSFUp", subJetTightDeepJetSFUp, "SubJet_tightDeepJetSFUp[SubJet_Size]/F");
                    tree->Branch("SubJet_tightDeepJetSFDown", subJetTightDeepJetSFDown, "SubJet_tightDeepJetSFDown[SubJet_Size]/F");
                }
            }
        }
    }

    if(name == "Isotrack"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("IsoTrack_Size", &isotrkSize, "IsoTrack_Size/S");

            tree->Branch("IsoTrack_PDG", isotrkPDG, "IsoTrack_PDG[IsoTrack_Size]/S");
            tree->Branch("IsoTrack_Charge", isotrkCharge, "IsoTrack_Charge[IsoTrack_Size]/S");
            tree->Branch("IsoTrack_Pt", isotrkPt, "IsoTrack_Pt[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Eta", isotrkEta, "IsoTrack_Eta[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Phi", isotrkPhi, "IsoTrack_Phi[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Dxy", isotrkDxy, "IsoTrack_Dxy[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Dz", isotrkDz, "IsoTrack_Dz[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Isolation", isotrkIso, "IsoTrack_Isolation[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_MiniIsolation", isotrkMiniIso, "IsoTrack_MiniIsolation[IsoTrack_Size]/F");
        }
    }

    if(name == "Misc"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Misc_eventNumber", &evNr, "Misc_eventNumber/I");
            if(!isData) tree->Branch("Misc_nParton", &nParton, "Misc_nParton/S");
        }
    }
}
