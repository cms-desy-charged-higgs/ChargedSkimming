#include <ChargedSkimming/Core/interface/output.h>

void Output::RegisterTrigger(const std::vector<std::string>& triggerNames, const std::vector<std::shared_ptr<TTree>>& trees){
    triggers = std::vector<short>(triggerNames.size(), 1);

    for(std::size_t i = 0; i < triggerNames.size(); ++i){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch(triggerNames[i].c_str(), &triggers[i], (triggerNames[i] + "/S").c_str());
        }
    }
}

void Output::Register(const std::string& name, const std::vector<std::shared_ptr<TTree>>& trees, pt::ptree& skim, const bool& isData){
    if(name == "Weight"){
        for(const std::shared_ptr<TTree>& tree: trees){
            if(!isData){
                tree->Branch("Weight_nTrueInt", &nTrueInt, "Weight_nTrueInt/S");

                tree->Branch("Weight_pdf", pdfWeight, "Weight_pdf[102]/F");
                tree->Branch("Weight_scale", scaleWeight, "Weight_scale[8]/F");
                
                tree->Branch("Weight_L1PreFire", &preFire, "Weight_L1PreFire/F");
                tree->Branch("Weight_L1PreFireUp", &preFireUp, "Weight_L1PreFireUp/F");
                tree->Branch("Weight_L1PreFireDown", &preFireDown, "Weight_L1PreFireDown/F");
            }
        }
    }

    if(name == "Electron"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Electron_Size", &nElectrons, "Electron_Size/S");

            tree->Branch("Electron_Pt", elePt.data(), "Electron_Pt[Electron_Size]/F");
            tree->Branch("Electron_Pt_eleEnergyScaleUp", elePtEnergyScaleUp.data(), "Electron_Pt_eleEnergyScaleUp[Electron_Size]/F");
            tree->Branch("Electron_Pt_eleEnergyScaleDown", elePtEnergyScaleDown.data(), "Electron_Pt_eleEnergyScaleDown[Electron_Size]/F");
            tree->Branch("Electron_Pt_eleEnergySigmaUp", elePtEnergySigmaUp.data(), "Electron_Pt_eleEnergySigmaUp[Electron_Size]/F");
            tree->Branch("Electron_Pt_eleEnergySigmaDown", elePtEnergySigmaDown.data(), "Electron_Pt_eleEnergySigmaDown[Electron_Size]/F");
            tree->Branch("Electron_Eta", eleEta.data(), "Electron_Eta[Electron_Size]/F");
            tree->Branch("Electron_Phi", elePhi.data(), "Electron_Phi[Electron_Size]/F");
            tree->Branch("Electron_Isolation03", eleIso03.data(), "Electron_Isolation03[Electron_Size]/F");
            tree->Branch("Electron_MiniIsolation", eleMiniIso.data(), "Electron_MiniIsolation[Electron_Size]/F");
            tree->Branch("Electron_Dxy", eleDxy.data(), "Electron_Dxy[Electron_Size]/F");
            tree->Branch("Electron_Dz", eleDz.data(), "Electron_Dz[Electron_Size]/F");
            tree->Branch("Electron_JetRelIsolation", eleRelJetIso.data(), "Electron_JetRelIsolation[Electron_Size]/F");

            tree->Branch("Electron_Charge", eleCharge.data(), "Electron_Charge[Electron_Size]/S");
            tree->Branch("Electron_CutID", eleCutID.data(), "Electron_CutID[Electron_Size]/S");
            tree->Branch("Electron_MVAID", eleMVAID.data(), "Electron_MVAID[Electron_Size]/S");

            if(!isData){
                tree->Branch("Electron_GenPt", eleGenPt.data(), "Electron_GenPt[Electron_Size]/F");
                tree->Branch("Electron_GenEta", eleGenEta.data(), "Electron_GenEta[Electron_Size]/F");
                tree->Branch("Electron_GenPhi", eleGenPhi.data(), "Electron_GenPhi[Electron_Size]/F");
                tree->Branch("Electron_GenID", eleGenID.data(), "Electron_GenID[Electron_Size]/S");
                tree->Branch("Electron_GenMotherID", eleGenMotherID.data(), "Electron_GenMotherID[Electron_Size]/S");
                tree->Branch("Electron_GenGrandMotherID", eleGenGrandMotherID.data(), "Electron_GenGrandMotherID[Electron_Size]/S");

                tree->Branch("Electron_RecoSF", eleRecoSF.data(), "Electron_RecoSF[Electron_Size]/F");
                tree->Branch("Electron_looseSF", eleLooseSF.data(), "Electron_looseSF[Electron_Size]/F");
                tree->Branch("Electron_mediumSF", eleMediumSF.data(), "Electron_mediumSF[Electron_Size]/F");
                tree->Branch("Electron_tightSF", eleTightSF.data(), "Electron_tightSF[Electron_Size]/F");
                tree->Branch("Electron_mediumMVASF", eleMediumMVASF.data(), "Electron_mediumMVASF[Electron_Size]/F");
                tree->Branch("Electron_tightMVASF", eleTightMVASF.data(), "Electron_tightMVASF[Electron_Size]/F");

                tree->Branch("Electron_RecoSFUp", eleRecoSFUp.data(), "Electron_RecoSFUp[Electron_Size]/F");
                tree->Branch("Electron_looseSFUp", eleLooseSFUp.data(), "Electron_looseSFUp[Electron_Size]/F");
                tree->Branch("Electron_mediumSFUp", eleMediumSFUp.data(), "Electron_mediumSFUp[Electron_Size]/F");
                tree->Branch("Electron_tightSFUp", eleTightSFUp.data(), "Electron_tightSFUp[Electron_Size]/F");
                tree->Branch("Electron_mediumMVASFUp", eleMediumMVASFUp.data(), "Electron_mediumMVASFUp[Electron_Size]/F");
                tree->Branch("Electron_tightMVASFUp", eleTightMVASFUp.data(), "Electron_tightMVASFUp[Electron_Size]/F");
                tree->Branch("Electron_RecoSFDown", eleRecoSFDown.data(), "Electron_RecoSFDown[Electron_Size]/F");
                tree->Branch("Electron_looseSFDown", eleLooseSFDown.data(), "Electron_looseSFDown[Electron_Size]/F");
                tree->Branch("Electron_mediumSFDown", eleMediumSFDown.data(), "Electron_mediumSFDown[Electron_Size]/F");
                tree->Branch("Electron_tightSFDown", eleTightSFDown.data(), "Electron_tightSFDown[Electron_Size]/F");
                tree->Branch("Electron_mediumMVASFDown", eleMediumMVASFDown.data(), "Electron_mediumMVASFDown[Electron_Size]/F");
                tree->Branch("Electron_tightMVASFDown", eleTightMVASFDown.data(), "Electron_tightMVASFDown[Electron_Size]/F");
            }
        }
    }

    if(name == "Muon"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Muon_Size", &nMuons, "Muon_Size/S");

            tree->Branch("Muon_Pt", muPt.data(), "Muon_Pt[Muon_Size]/F");
            tree->Branch("Muon_Pt_muMomentumScaleUp", muPtUp.data(), "Muon_Pt_muMomentumScaleUp[Muon_Size]/F");
            tree->Branch("Muon_Pt_muMomentumScaleDown", muPtDown.data(), "Muon_Pt_muMomentumScaleDown[Muon_Size]/F");
            tree->Branch("Muon_Eta", muEta.data(), "Muon_Eta[Muon_Size]/F");
            tree->Branch("Muon_Phi", muPhi.data(), "Muon_Phi[Muon_Size]/F");
            tree->Branch("Muon_Isolation03", muIso03.data(), "Muon_Isolation03[Muon_Size]/F");
            tree->Branch("Muon_Isolation04", muIso04.data(), "Muon_Isolation04[Muon_Size]/F");
            tree->Branch("Muon_MiniIsolation", muMiniIso.data(), "Muon_MiniIsolation[Muon_Size]/F");
            tree->Branch("Muon_Dxy", muDxy.data(), "Muon_Dxy[Muon_Size]/F");
            tree->Branch("Muon_Dz", muDz.data(), "Muon_Dz[Muon_Size]/F");
            tree->Branch("Muon_JetRelIsolation", muRelJetIso.data(), "Muon_JetRelIsolation[Muon_Size]/F");

            tree->Branch("Muon_Charge", muCharge.data(), "Muon_Charge[Muon_Size]/S");
            tree->Branch("Muon_CutID", muCutID.data(), "Muon_CutID[Muon_Size]/S");
            tree->Branch("Muon_MVAID", muMVAID.data(), "Muon_MVAID[Muon_Size]/S");

            if(!isData){
                tree->Branch("Muon_GenPt", muGenPt.data(), "Muon_GenPt[Muon_Size]/F");
                tree->Branch("Muon_GenEta", muGenEta.data(), "Muon_GenEta[Muon_Size]/F");
                tree->Branch("Muon_GenPhi", muGenPhi.data(), "Muon_GenPhi[Muon_Size]/F");
                tree->Branch("Muon_GenID", muGenID.data(), "Muon_GenID[Muon_Size]/S");
                tree->Branch("Muon_GenMotherID", muGenMotherID.data(), "Muon_GenMotherID[Muon_Size]/S");
                tree->Branch("Muon_GenGrandMotherID", muGenGrandMotherID.data(), "Muon_GenGrandMotherID[Muon_Size]/S");

                tree->Branch("Muon_TriggerSF", muTriggerSF.data(), "Muon_TriggerSF[Muon_Size]/F");
                tree->Branch("Muon_looseIsoSF", muLooseIsoSF.data(), "Muon_looseIsoSF[Muon_Size]/F");
                tree->Branch("Muon_tightIsoSF", muTightIsoSF.data(), "Muon_tightIsoSF[Muon_Size]/F");
                tree->Branch("Muon_looseSF", muLooseSF.data(), "Muon_looseSF[Muon_Size]/F");
                tree->Branch("Muon_mediumSF", muMediumSF.data(), "Muon_mediumSF[Muon_Size]/F");
                tree->Branch("Muon_tightSF", muTightSF.data(), "Muon_tightSF[Muon_Size]/F");

                tree->Branch("Muon_TriggerSFUp", muTriggerSFUp.data(), "Muon_TriggerSFUp[Muon_Size]/F");
                tree->Branch("Muon_looseIsoSFUp", muLooseIsoSFUp.data(), "Muon_looseIsoSFUp[Muon_Size]/F");
                tree->Branch("Muon_tightIsoSFUp", muTightIsoSFUp.data(), "Muon_tightIsoSFUp[Muon_Size]/F");
                tree->Branch("Muon_looseSFUp", muLooseSFUp.data(), "Muon_looseSFUp[Muon_Size]/F");
                tree->Branch("Muon_mediumSFUp", muMediumSFUp.data(), "Muon_mediumSFUp[Muon_Size]/F");
                tree->Branch("Muon_tightSFUp", muTightSFUp.data(), "Muon_tightSFUp[Muon_Size]/F");
                tree->Branch("Muon_TriggerSFDown", muTriggerSFDown.data(), "Muon_TriggerSFDown[Muon_Size]/F");
                tree->Branch("Muon_looseIsoSFDown", muLooseIsoSFDown.data(), "Muon_looseIsoSFDown[Muon_Size]/F");
                tree->Branch("Muon_tightIsoSFDown", muTightIsoSFDown.data(), "Muon_tightIsoSFDown[Muon_Size]/F");
                tree->Branch("Muon_looseSFDown", muLooseSFDown.data(), "Muon_looseSFDown[Muon_Size]/F");
                tree->Branch("Muon_mediumSFDown", muMediumSFDown.data(), "Muon_mediumSFDown[Muon_Size]/F");
                tree->Branch("Muon_tightSFDown", muTightSFDown.data(), "Muon_tightSFDown[Muon_Size]/F");
            }
        }
    }

    if(name == "Jet"){
        if(!isData){
            for(const std::pair<std::string, boost::property_tree::ptree> j : skim.get_child("Analyzer.Jet.JECSyst")){
                jetPtJECUp.push_back(std::array<float, jetMax>());
                jetPtJECDown.push_back(std::array<float, jetMax>());
                jetMassJECUp.push_back(std::array<float, jetMax>());
                jetMassJECDown.push_back(std::array<float, jetMax>());
                
                fatJetPtJECUp.push_back(std::array<float, fatJetMax>());
                fatJetPtJECDown.push_back(std::array<float, fatJetMax>());
                fatJetMassJECUp.push_back(std::array<float, fatJetMax>());
                fatJetMassJECDown.push_back(std::array<float, fatJetMax>());
                
                subJetPtJECUp.push_back(std::array<float, jetMax>());
                subJetPtJECDown.push_back(std::array<float, jetMax>());
                subJetMassJECUp.push_back(std::array<float, jetMax>());
                subJetMassJECDown.push_back(std::array<float, jetMax>());

                metPtJECUp.push_back(1);
                metPtJECDown.push_back(1);
                metPhiJECUp.push_back(1);
                metPhiJECDown.push_back(1);
            }

            for(const std::pair<std::string, boost::property_tree::ptree> j : skim.get_child("Analyzer.Jet.BTagSyst")){
                jetLooseDeepCSVSFDown.push_back(std::array<float, jetMax>());
                jetLooseDeepCSVSFUp.push_back(std::array<float, jetMax>());
                jetMediumDeepCSVSFDown.push_back(std::array<float, jetMax>());
                jetMediumDeepCSVSFUp.push_back(std::array<float, jetMax>());
                jetTightDeepCSVSFDown.push_back(std::array<float, jetMax>());
                jetTightDeepCSVSFUp.push_back(std::array<float, jetMax>());
                jetLooseDeepJetSFDown.push_back(std::array<float, jetMax>());
                jetLooseDeepJetSFUp.push_back(std::array<float, jetMax>());
                jetMediumDeepJetSFDown.push_back(std::array<float, jetMax>());
                jetMediumDeepJetSFUp.push_back(std::array<float, jetMax>());
                jetTightDeepJetSFDown.push_back(std::array<float, jetMax>()); 
                jetTightDeepJetSFUp.push_back(std::array<float, jetMax>());

                subJetLooseDeepCSVSFDown.push_back(std::array<float, jetMax>());
                subJetLooseDeepCSVSFUp.push_back(std::array<float, jetMax>());
                subJetMediumDeepCSVSFDown.push_back(std::array<float, jetMax>());
                subJetMediumDeepCSVSFUp.push_back(std::array<float, jetMax>());
                subJetTightDeepCSVSFDown.push_back(std::array<float, jetMax>());
                subJetTightDeepCSVSFUp.push_back(std::array<float, jetMax>());
                subJetLooseDeepJetSFDown.push_back(std::array<float, jetMax>());
                subJetLooseDeepJetSFUp.push_back(std::array<float, jetMax>());
                subJetMediumDeepJetSFDown.push_back(std::array<float, jetMax>());
                subJetMediumDeepJetSFUp.push_back(std::array<float, jetMax>());
                subJetTightDeepJetSFDown.push_back(std::array<float, jetMax>()); 
                subJetTightDeepJetSFUp.push_back(std::array<float, jetMax>());
            }

            for(const std::pair<std::string, boost::property_tree::ptree> j : skim.get_child("Analyzer.Jet.BTagSystLight")){
                jetLooseDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                jetLooseDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                jetMediumDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                jetMediumDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                jetTightDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                jetTightDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                jetLooseDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                jetLooseDeepJetSFLightUp.push_back(std::array<float, jetMax>());
                jetMediumDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                jetMediumDeepJetSFLightUp.push_back(std::array<float, jetMax>());  
                jetTightDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                jetTightDeepJetSFLightUp.push_back(std::array<float, jetMax>());

                subJetLooseDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                subJetLooseDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                subJetMediumDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                subJetMediumDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                subJetTightDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                subJetTightDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                subJetLooseDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                subJetLooseDeepJetSFLightUp.push_back(std::array<float, jetMax>());
                subJetMediumDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                subJetMediumDeepJetSFLightUp.push_back(std::array<float, jetMax>());  
                subJetTightDeepJetSFLightDown.push_back(std::array<float, jetMax>());
                subJetTightDeepJetSFLightUp.push_back(std::array<float, jetMax>()); 
            }
        }

        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Jet_Size", &nJets, "Jet_Size/S");
            tree->Branch("SubJet_Size", &nSubJets, "SubJet_Size/S");
            tree->Branch("FatJet_Size", &nFatJets, "FatJet_Size/S");

            tree->Branch("MET_Pt", &metPt, "MET_Pt/F");
            tree->Branch("MET_Pt_JMEUp", &metPtJMEUp,"MET_Pt_JMEUp/F");
            tree->Branch("MET_Pt_JMEDown", &metPtJMEDown, "MET_Pt_JMEDown/F");
            tree->Branch("MET_Pt_UnclusteredUp", &metPtUnclusteredUp, "MET_Pt_UnclusteredUp/F");
            tree->Branch("MET_Pt_UnclusteredDown", &metPtUnclusteredDown, "MET_Pt_UnclusteredDown/F");
            tree->Branch("MET_Phi", &metPhi, "MET_Phi/F");
            tree->Branch("MET_Phi_JMEUp", &metPhiJMEUp, "MET_Phi_JMEUp/F");
            tree->Branch("MET_Phi_JMEDown", &metPhiJMEDown, "MET_Phi_JMEDown/F");
            tree->Branch("MET_Phi_UnclusteredUp", &metPhiUnclusteredUp, "MET_Phi_UnclusteredUp/F");
            tree->Branch("MET_Phi_UnclusteredDown", &metPhiUnclusteredDown, "MET_Phi_UnclusteredDown/F");

            tree->Branch("Jet_Pt", jetPt.data(), "Jet_Pt[Jet_Size]/F");
            tree->Branch("Jet_Pt_JMEUp", jetPtJMEUp.data(), "Jet_Pt_JMEUp[Jet_Size]/F");
            tree->Branch("Jet_Pt_JMEDown", jetPtJMEDown.data(), "Jet_Pt_JMEDown[Jet_Size]/F");
            tree->Branch("Jet_Eta", jetEta.data(), "Jet_Eta[Jet_Size]/F");
            tree->Branch("Jet_Phi", jetPhi.data(), "Jet_Phi[Jet_Size]/F");
            tree->Branch("Jet_Mass", jetMass.data(), "Jet_Mass[Jet_Size]/F");
            tree->Branch("Jet_Mass_JMEUp", jetMassJMEUp.data(), "Jet_Mass_JMEUp[Jet_Size]/F");
            tree->Branch("Jet_Mass_JMEDown", jetMassJMEDown.data(), "Jet_Mass_JMEDown[Jet_Size]/F");
            tree->Branch("Jet_DeepJet", jetDeepJet.data(), "Jet_DeepJet[Jet_Size]/F");
            tree->Branch("Jet_DeepCSV", jetDeepCSV.data(), "Jet_DeepCSV[Jet_Size]/F");
            tree->Branch("Jet_DeepJetID", jetDeepJetID.data(), "Jet_DeepJetID[Jet_Size]/S");
            tree->Branch("Jet_DeepCSVID", jetDeepCSVID.data(), "Jet_DeepCSVID[Jet_Size]/S");
            tree->Branch("Jet_JECFac", jetJEC.data(), "Jet_JECFac[Jet_Size]/F");
            tree->Branch("Jet_JMEFac", jetJME.data(), "Jet_JMEFac[Jet_Size]/F");
            tree->Branch("Jet_ID", jetID.data(), "Jet_ID[Jet_Size]/S");
            tree->Branch("Jet_PUID", jetPUID.data(), "Jet_PUID[Jet_Size]/S");

            tree->Branch("SubJet_Pt", subJetPt.data(), "SubJet_Pt[SubJet_Size]/F");
            tree->Branch("SubJet_Pt_JMEUp", subJetPtJMEUp.data(), "SubJet_Pt_JMEUp[SubJet_Size]/F");
            tree->Branch("SubJet_Pt_JMEDown", subJetPtJMEDown.data(), "SubJet_Pt_JMEDown[SubJet_Size]/F");
            tree->Branch("SubJet_Eta", subJetEta.data(), "SubJet_Eta[SubJet_Size]/F");
            tree->Branch("SubJet_Phi", subJetPhi.data(), "SubJet_Phi[SubJet_Size]/F");
            tree->Branch("SubJet_Mass", subJetMass.data(), "SubJet_Mass[SubJet_Size]/F");
            tree->Branch("SubJet_Mass_JMEUp", subJetMassJMEUp.data(), "SubJet_Mass_JMEUp[SubJet_Size]/F");
            tree->Branch("SubJet_Mass_JMEDown", subJetMassJMEDown.data(), "SubJet_Mass_JMEDown[SubJet_Size]/F");
            tree->Branch("SubJet_DeepJet", subJetDeepJet.data(), "SubJet_DeepJet[SubJet_Size]/F");
            tree->Branch("SubJet_DeepCSV", subJetDeepCSV.data(), "SubJet_DeepCSV[SubJet_Size]/F");
            tree->Branch("SubJet_DeepJetID", subJetDeepJetID.data(), "SubJet_DeepJetID[SubJet_Size]/S");
            tree->Branch("SubJet_DeepCSVID", subJetDeepCSVID.data(), "SubJet_DeepCSVID[SubJet_Size]/S");
            tree->Branch("SubJet_PartonFlavour", subJetPartFlav.data(), "SubJet_PartonFlavour[SubJet_Size]/S");
            tree->Branch("SubJet_JECFac", subJetJEC.data(), "SubJet_JECFac[SubJet_Size]/F");
            tree->Branch("SubJet_JMEFac", subJetJME.data(), "SubJet_JMEFac[SubJet_Size]/F");
            tree->Branch("SubJet_FatJetIdx", fatJetIdx.data(), "SubJet_FatJetIdx[SubJet_Size]/S");

            tree->Branch("FatJet_Pt", fatJetPt.data(), "FatJet_Pt[FatJet_Size]/F");
            tree->Branch("FatJet_Pt_JMEUp", fatJetPtJMEUp.data(), "FatJet_Pt_JMEUp[FatJet_Size]/F");
            tree->Branch("FatJet_Pt_JMEDown", fatJetPtJMEDown.data(), "FatJet_Pt_JMEDown[FatJet_Size]/F");
            tree->Branch("FatJet_Eta", fatJetEta.data(), "FatJet_Eta[FatJet_Size]/F");
            tree->Branch("FatJet_Phi", fatJetPhi.data(), "FatJet_Phi[FatJet_Size]/F");
            tree->Branch("FatJet_Mass", fatJetMass.data(), "FatJet_Mass[FatJet_Size]/F");
            tree->Branch("FatJet_Mass_JMEUp", fatJetMassJMEUp.data(), "FatJet_Mass_JMEUp[FatJet_Size]/F");
            tree->Branch("FatJet_Mass_JMEDown", fatJetMassJMEDown.data(), "FatJet_Mass_JMEDown[FatJet_Size]/F");
            tree->Branch("FatJet_Tau1", fatJetTau1.data(), "FatJet_Tau1[FatJet_Size]/F");
            tree->Branch("FatJet_Tau2", fatJetTau2.data(), "FatJet_Tau2[FatJet_Size]/F");
            tree->Branch("FatJet_Tau3", fatJetTau3.data(), "FatJet_Tau3[FatJet_Size]/F");
            tree->Branch("FatJet_DeepAK8ID",fatJetDAK8ID.data(), "FatJet_DeepAK8ID[FatJet_Size]/S");
            tree->Branch("FatJet_JECFac", fatJetJEC.data(), "FatJet_JECFac[FatJet_Size]/F");
            tree->Branch("FatJet_JMEFac", fatJetJME.data(), "FatJet_JMEFac[FatJet_Size]/F");
            
            if(!isData){
                std::size_t JECIdx = -1;
                
                for(const std::pair<std::string, boost::property_tree::ptree> j : skim.get_child("Analyzer.Jet.JECSyst")){
                    const std::string JECSyst = j.second.get_value<std::string>();
                    ++JECIdx;            
                
                    tree->Branch(("Jet_Pt_JEC" + JECSyst + "Up").c_str(), jetPtJECUp[JECIdx].data(), ("Jet_Pt_JEC" + JECSyst + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_Pt_JEC" + JECSyst + "Down").c_str(), jetPtJECDown[JECIdx].data(), ("Jet_Pt_JEC" + JECSyst + "Down[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_Mass_JEC" + JECSyst + "Up").c_str(), jetMassJECUp[JECIdx].data(), ("Jet_Mass_JEC" + JECSyst + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_Mass_JEC" + JECSyst + "Down").c_str(), jetMassJECDown[JECIdx].data(), ("Jet_Mass_JEC" + JECSyst + "Down[Jet_Size]/F").c_str());
                    
                    tree->Branch(("FatJet_Pt_JEC" + JECSyst + "Up").c_str(), fatJetPtJECUp[JECIdx].data(), ("FatJet_Pt_JEC" + JECSyst + "Up[FatJet_Size]/F").c_str());
                    tree->Branch(("FatJet_Pt_JEC" + JECSyst + "Down").c_str(), fatJetPtJECDown[JECIdx].data(), ("FatJet_Pt_JEC" + JECSyst + "Down[FatJet_Size]/F").c_str());
                    tree->Branch(("FatJet_Mass_JEC" + JECSyst + "Up").c_str(), fatJetMassJECUp[JECIdx].data(), ("FatJet_Mass_JEC" + JECSyst + "Up[FatJet_Size]/F").c_str());
                    tree->Branch(("FatJet_Mass_JEC" + JECSyst + "Down").c_str(), fatJetMassJECDown[JECIdx].data(), ("FatJet_Mass_JEC" + JECSyst + "Down[FatJet_Size]/F").c_str());
                    
                    tree->Branch(("SubJet_Pt_JEC" + JECSyst + "Up").c_str(), subJetPtJECUp[JECIdx].data(), ("SubJet_Pt_JEC" + JECSyst + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_Pt_JEC" + JECSyst + "Down").c_str(), subJetPtJECDown[JECIdx].data(), ("SubJet_Pt_JEC" + JECSyst + "Down[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_Mass_JEC" + JECSyst + "Up").c_str(), subJetMassJECUp[JECIdx].data(), ("SubJet_Mass_JEC" + JECSyst + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_Mass_JEC" + JECSyst + "Down").c_str(), subJetMassJECDown[JECIdx].data(), ("SubJet_Mass_JEC" + JECSyst + "Down[SubJet_Size]/F").c_str());
                    
                    tree->Branch(("MET_Pt_JEC" + JECSyst + "Up").c_str(), &metPtJECUp[JECIdx], ("MET_Pt_JEC" + JECSyst + "Up/F").c_str());
                    tree->Branch(("MET_Pt_JEC" + JECSyst + "Down").c_str(), &metPtJECDown[JECIdx], ("MET_Pt_JEC" + JECSyst + "Down/F").c_str());
                    tree->Branch(("MET_Phi_JEC" + JECSyst + "Up").c_str(), &metPhiJECUp[JECIdx], ("MET_Phi_JEC" + JECSyst + "Up/F").c_str());
                    tree->Branch(("MET_Phi_JEC" + JECSyst + "Down").c_str(), &metPhiJECDown[JECIdx], ("MET_Phi_JEC" + JECSyst + "Down/F").c_str());  
                }

                tree->Branch("Jet_PartonFlavour", jetPartFlav.data(), "Jet_PartonFlavour[Jet_Size]/S");
                tree->Branch("Jet_GenPt", jetGenPt.data(), "Jet_GenPt[Jet_Size]/F");
                tree->Branch("Jet_GenEta", jetGenEta.data(), "Jet_GenEta[Jet_Size]/F");
                tree->Branch("Jet_GenPhi", jetGenPhi.data(), "Jet_GenPhi[Jet_Size]/F");
                tree->Branch("Jet_GenID", jetGenID.data(), "Jet_GenID[Jet_Size]/S");
                tree->Branch("Jet_GenMotherID", jetGenMotherID.data(), "Jet_GenMotherID[Jet_Size]/S");
                tree->Branch("Jet_GenGrandMotherID", jetGenGrandMotherID.data(), "Jet_GenGrandMotherID[Jet_Size]/S");

                tree->Branch("SubJet_GenPt", subJetGenPt.data(), "SubJet_GenPt[SubJet_Size]/F");
                tree->Branch("SubJet_GenEta", subJetGenEta.data(), "SubJet_GenEta[SubJet_Size]/F");
                tree->Branch("SubJet_GenPhi", subJetGenPhi.data(), "SubJet_GenPhi[SubJet_Size]/F");
                tree->Branch("SubJet_GenID", subJetGenID.data(), "SubJet_GenID[SubJet_Size]/S");
                tree->Branch("SubJet_GenMotherID", subJetGenMotherID.data(), "SubJet_GenMotherID[SubJet_Size]/S");
                tree->Branch("SubJet_GenGrandMotherID", subJetGenGrandMotherID.data(), "SubJet_GenGrandMotherID[SubJet_Size]/S");

                tree->Branch("Jet_looseDeepCSVSF", jetLooseDeepCSVSF.data(), "Jet_looseDeepCSVSF[Jet_Size]/F");
                tree->Branch("Jet_mediumDeepCSVSF", jetMediumDeepCSVSF.data(), "Jet_mediumDeepCSVSF[Jet_Size]/F");
                tree->Branch("Jet_tightDeepCSVSF", jetTightDeepCSVSF.data(), "Jet_tightDeepCSVSF[Jet_Size]/F");
                tree->Branch("Jet_looseDeepJetSF", jetLooseDeepJetSF.data(), "Jet_looseDeepJetSF[Jet_Size]/F");
                tree->Branch("Jet_mediumDeepJetSF", jetMediumDeepJetSF.data(), "Jet_mediumDeepJetSF[Jet_Size]/F");
                tree->Branch("Jet_tightDeepJetSF", jetTightDeepJetSF.data(), "Jet_tightDeepJetSF[Jet_Size]/F");

                tree->Branch("SubJet_looseDeepCSVSF", subJetLooseDeepCSVSF.data(), "SubJet_looseDeepCSVSF[SubJet_Size]/F");
                tree->Branch("SubJet_mediumDeepCSVSF", subJetMediumDeepCSVSF.data(), "SubJet_mediumDeepCSVSF[SubJet_Size]/F");
                tree->Branch("SubJet_tightDeepCSVSF", subJetTightDeepCSVSF.data(), "SubJet_tightDeepCSVSF[SubJet_Size]/F");
                tree->Branch("SubJet_looseDeepJetSF", subJetLooseDeepJetSF.data(), "SubJet_looseDeepJetSF[SubJet_Size]/F");
                tree->Branch("SubJet_mediumDeepJetSF", subJetMediumDeepJetSF.data(), "SubJet_mediumDeepJetSF[SubJet_Size]/F");
                tree->Branch("SubJet_tightDeepJetSF", subJetTightDeepJetSF.data(), "SubJet_tightDeepJetSF[SubJet_Size]/F");

                std::size_t bTagUncIdx = -1;

                for(const std::pair<std::string, boost::property_tree::ptree> j : skim.get_child("Analyzer.Jet.BTagSyst")){
                    std::string bUnc = j.second.get_value<std::string>();
                    bUnc[0] = std::toupper(bUnc[0]);
                    ++bTagUncIdx;

                    tree->Branch(("Jet_looseDeepCSVSF" + bUnc + "Up").c_str(), jetLooseDeepCSVSFUp[bTagUncIdx].data(), ("Jet_looseDeepCSVSF" + bUnc + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_looseDeepCSVSF" + bUnc + "Down").c_str(), jetLooseDeepCSVSFDown[bTagUncIdx].data(), ("Jet_looseDeepCSVSF" + bUnc + "Down[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepCSVSF" + bUnc + "Up").c_str(), jetMediumDeepCSVSFUp[bTagUncIdx].data(), ("Jet_mediumDeepCSVSF" + bUnc + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepCSVSF" + bUnc + "Down").c_str(), jetMediumDeepCSVSFDown[bTagUncIdx].data(), ("Jet_mediumDeepCSVSF" + bUnc + "Down[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepCSVSF" + bUnc + "Up").c_str(), jetTightDeepCSVSFUp[bTagUncIdx].data(), ("Jet_tightDeepCSVSF" + bUnc + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepCSVSF" + bUnc + "Down").c_str(), jetTightDeepCSVSFDown[bTagUncIdx].data(), ("Jet_tightDeepCSVSF" + bUnc + "Down[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_looseDeepJetSF" + bUnc + "Up").c_str(), jetLooseDeepJetSFUp[bTagUncIdx].data(), ("Jet_looseDeepJetSF" + bUnc + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_looseDeepJetSF" + bUnc + "Down").c_str(), jetLooseDeepJetSFDown[bTagUncIdx].data(), ("Jet_looseDeepJetSF" + bUnc + "Down[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepJetSF" + bUnc + "Up").c_str(), jetMediumDeepJetSFUp[bTagUncIdx].data(), ("Jet_mediumDeepJetSF" + bUnc + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepJetSF" + bUnc + "Down").c_str(), jetMediumDeepJetSFDown[bTagUncIdx].data(), ("Jet_mediumDeepJetSF" + bUnc + "Down[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepJetSF" + bUnc + "Up").c_str(), jetTightDeepJetSFUp[bTagUncIdx].data(), ("Jet_tightDeepJetSF" + bUnc + "Up[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepJetSF" + bUnc + "Down").c_str(), jetTightDeepJetSFDown[bTagUncIdx].data(), ("Jet_tightDeepJetSF" + bUnc + "Down[Jet_Size]/F").c_str());

                    tree->Branch(("SubJet_looseDeepCSVSF" + bUnc + "Up").c_str(), subJetLooseDeepCSVSFUp[bTagUncIdx].data(), ("SubJet_looseDeepCSVSF" + bUnc + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_looseDeepCSVSF" + bUnc + "Down").c_str(), subJetLooseDeepCSVSFDown[bTagUncIdx].data(), ("SubJet_looseDeepCSVSF" + bUnc + "Down[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepCSVSF" + bUnc + "Up").c_str(), subJetMediumDeepCSVSFUp[bTagUncIdx].data(), ("SubJet_mediumDeepCSVSF" + bUnc + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepCSVSF" + bUnc + "Down").c_str(), subJetMediumDeepCSVSFDown[bTagUncIdx].data(), ("SubJet_mediumDeepCSVSF" + bUnc + "Down[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepCSVSF" + bUnc + "Up").c_str(), subJetTightDeepCSVSFUp[bTagUncIdx].data(), ("SubJet_tightDeepCSVSF" + bUnc + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepCSVSF" + bUnc + "Down").c_str(), subJetTightDeepCSVSFDown[bTagUncIdx].data(), ("SubJet_tightDeepCSVSF" + bUnc + "Down[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_looseDeepJetSF" + bUnc + "Up").c_str(), subJetLooseDeepJetSFUp[bTagUncIdx].data(), ("SubJet_looseDeepJetSF" + bUnc + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_looseDeepJetSF" + bUnc + "Down").c_str(), subJetLooseDeepJetSFDown[bTagUncIdx].data(), ("SubJet_looseDeepJetSF" + bUnc + "Down[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepJetSF" + bUnc + "Up").c_str(), subJetMediumDeepJetSFUp[bTagUncIdx].data(), ("SubJet_mediumDeepJetSF" + bUnc + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepJetSF" + bUnc + "Down").c_str(), subJetMediumDeepJetSFDown[bTagUncIdx].data(), ("SubJet_mediumDeepJetSF" + bUnc + "Down[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepJetSF" + bUnc + "Up").c_str(), subJetTightDeepJetSFUp[bTagUncIdx].data(), ("SubJet_tightDeepJetSF" + bUnc + "Up[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepJetSF" + bUnc + "Down").c_str(), subJetTightDeepJetSFDown[bTagUncIdx].data(), ("SubJet_tightDeepJetSF" + bUnc + "Down[SubJet_Size]/F").c_str());
                }

                bTagUncIdx = 0;

                for(const std::pair<std::string, boost::property_tree::ptree> j : skim.get_child("Analyzer.Jet.BTagSystLight")){
                    std::string bUnc = j.second.get_value<std::string>();
                    bUnc[0] = std::toupper(bUnc[0]);
                    ++bTagUncIdx;

                    jetLooseDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                    jetLooseDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                    jetMediumDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                    jetMediumDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                    jetTightDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                    jetTightDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                    jetLooseDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                    jetLooseDeepJetSFLightUp.push_back(std::array<float, jetMax>());
                    jetMediumDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                    jetMediumDeepJetSFLightUp.push_back(std::array<float, jetMax>());  
                    jetTightDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                    jetTightDeepJetSFLightUp.push_back(std::array<float, jetMax>());

                    subJetLooseDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                    subJetLooseDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                    subJetMediumDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                    subJetMediumDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                    subJetTightDeepCSVSFLightDown.push_back(std::array<float, jetMax>());
                    subJetTightDeepCSVSFLightUp.push_back(std::array<float, jetMax>());
                    subJetLooseDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                    subJetLooseDeepJetSFLightUp.push_back(std::array<float, jetMax>());
                    subJetMediumDeepJetSFLightDown.push_back(std::array<float, jetMax>()); 
                    subJetMediumDeepJetSFLightUp.push_back(std::array<float, jetMax>());  
                    subJetTightDeepJetSFLightDown.push_back(std::array<float, jetMax>());
                    subJetTightDeepJetSFLightUp.push_back(std::array<float, jetMax>()); 

                    tree->Branch(("Jet_looseDeepCSVSF" + bUnc + "LightUp").c_str(), jetLooseDeepCSVSFLightUp[bTagUncIdx].data(), ("Jet_looseDeepCSVSF" + bUnc + "LightUp[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_looseDeepCSVSF" + bUnc + "LightDown").c_str(), jetLooseDeepCSVSFLightDown[bTagUncIdx].data(), ("Jet_looseDeepCSVSF" + bUnc + "LightDown[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepCSVSF" + bUnc + "LightUp").c_str(), jetMediumDeepCSVSFLightUp[bTagUncIdx].data(), ("Jet_mediumDeepCSVSF" + bUnc + "LightUp[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepCSVSF" + bUnc + "LightDown").c_str(), jetMediumDeepCSVSFLightDown[bTagUncIdx].data(), ("Jet_mediumDeepCSVSF" + bUnc + "LightDown[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepCSVSF" + bUnc + "LightUp").c_str(), jetTightDeepCSVSFLightUp[bTagUncIdx].data(), ("Jet_tightDeepCSVSF" + bUnc + "LightUp[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepCSVSF" + bUnc + "LightDown").c_str(), jetTightDeepCSVSFLightDown[bTagUncIdx].data(), ("Jet_tightDeepCSVSF" + bUnc + "LightDown[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_looseDeepJetSF" + bUnc + "LightUp").c_str(), jetLooseDeepJetSFLightUp[bTagUncIdx].data(), ("Jet_looseDeepJetSF" + bUnc + "LightUp[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_looseDeepJetSF" + bUnc + "LightDown").c_str(), jetLooseDeepJetSFLightDown[bTagUncIdx].data(), ("Jet_looseDeepJetSF" + bUnc + "LightDown[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepJetSF" + bUnc + "LightUp").c_str(), jetMediumDeepJetSFLightUp[bTagUncIdx].data(), ("Jet_mediumDeepJetSF" + bUnc + "LightUp[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_mediumDeepJetSF" + bUnc + "LightDown").c_str(), jetMediumDeepJetSFLightDown[bTagUncIdx].data(), ("Jet_mediumDeepJetSF" + bUnc + "LightDown[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepJetSF" + bUnc + "LightUp").c_str(), jetTightDeepJetSFLightUp[bTagUncIdx].data(), ("Jet_tightDeepJetSF" + bUnc + "LightUp[Jet_Size]/F").c_str());
                    tree->Branch(("Jet_tightDeepJetSF" + bUnc + "LightDown").c_str(), jetTightDeepJetSFLightDown[bTagUncIdx].data(), ("Jet_tightDeepJetSF" + bUnc + "LightDown[Jet_Size]/F").c_str());

                    tree->Branch(("SubJet_looseDeepCSVSF" + bUnc + "LightUp").c_str(), subJetLooseDeepCSVSFLightUp[bTagUncIdx].data(), ("SubJet_looseDeepCSVSF" + bUnc + "LightUp[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_looseDeepCSVSF" + bUnc + "LightDown").c_str(), subJetLooseDeepCSVSFLightDown[bTagUncIdx].data(), ("SubJet_looseDeepCSVSF" + bUnc + "LightDown[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepCSVSF" + bUnc + "LightUp").c_str(), subJetMediumDeepCSVSFLightUp[bTagUncIdx].data(), ("SubJet_mediumDeepCSVSF" + bUnc + "LightUp[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepCSVSF" + bUnc + "LightDown").c_str(), subJetMediumDeepCSVSFLightDown[bTagUncIdx].data(), ("SubJet_mediumDeepCSVSF" + bUnc + "LightDown[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepCSVSF" + bUnc + "LightUp").c_str(), subJetTightDeepCSVSFLightUp[bTagUncIdx].data(), ("SubJet_tightDeepCSVSF" + bUnc + "LightUp[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepCSVSF" + bUnc + "LightDown").c_str(), subJetTightDeepCSVSFLightDown[bTagUncIdx].data(), ("SubJet_tightDeepCSVSF" + bUnc + "LightDown[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_looseDeepJetSF" + bUnc + "LightUp").c_str(), subJetLooseDeepJetSFLightUp[bTagUncIdx].data(), ("SubJet_looseDeepJetSF" + bUnc + "LightUp[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_looseDeepJetSF" + bUnc + "LightDown").c_str(), subJetLooseDeepJetSFLightDown[bTagUncIdx].data(), ("SubJet_looseDeepJetSF" + bUnc + "LightDown[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepJetSF" + bUnc + "LightUp").c_str(), subJetMediumDeepJetSFLightUp[bTagUncIdx].data(), ("SubJet_mediumDeepJetSF" + bUnc + "LightUp[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_mediumDeepJetSF" + bUnc + "LightDown").c_str(), subJetMediumDeepJetSFLightDown[bTagUncIdx].data(), ("SubJet_mediumDeepJetSF" + bUnc + "LightDown[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepJetSF" + bUnc + "LightUp").c_str(), subJetTightDeepJetSFLightUp[bTagUncIdx].data(), ("SubJet_tightDeepJetSF" + bUnc + "LightUp[SubJet_Size]/F").c_str());
                    tree->Branch(("SubJet_tightDeepJetSF" + bUnc + "LightDown").c_str(), subJetTightDeepJetSFLightDown[bTagUncIdx].data(), ("SubJet_tightDeepJetSF" + bUnc + "LightDown[SubJet_Size]/F").c_str());
                }
            }
        }
    }

    if(name == "Isotrack"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("IsoTrack_Size", &isotrkSize, "IsoTrack_Size/S");

            tree->Branch("IsoTrack_PDG", isotrkPDG.data(), "IsoTrack_PDG[IsoTrack_Size]/S");
            tree->Branch("IsoTrack_Charge", isotrkCharge.data(), "IsoTrack_Charge[IsoTrack_Size]/S");
            tree->Branch("IsoTrack_Pt", isotrkPt.data(), "IsoTrack_Pt[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Eta", isotrkEta.data(), "IsoTrack_Eta[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Phi", isotrkPhi.data(), "IsoTrack_Phi[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Dxy", isotrkDxy.data(), "IsoTrack_Dxy[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Dz", isotrkDz.data(), "IsoTrack_Dz[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_Isolation03", isotrkIso03.data(), "IsoTrack_Isolation03[IsoTrack_Size]/F");
            tree->Branch("IsoTrack_MiniIsolation", isotrkMiniIso.data(), "IsoTrack_MiniIsolation[IsoTrack_Size]/F");
        }
    }

    if(name == "Misc"){
        for(const std::shared_ptr<TTree>& tree: trees){
            tree->Branch("Misc_eventNumber", &evNr, "Misc_eventNumber/I");
            if(!isData) tree->Branch("Misc_nParton", &nParton, "Misc_nParton/S");
        }
    }
}
