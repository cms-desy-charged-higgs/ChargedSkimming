#include <ChargedSkimming/Core/interface/nanoinput.h>

NanoInput::NanoInput(const std::string& fileName, const std::string& treeName){
    inputFile = std::shared_ptr<TFile>(TFile::Open(fileName.c_str(), "READ"));
    inputTree.reset(static_cast<TTree*>(inputFile->Get(treeName.c_str())));

    //Weight related
    pdfWeightL = inputTree->GetLeaf("LHEPdfWeight");
    scaleWeightL = inputTree->GetLeaf("LHEScaleWeight");
    nTrueIntL = inputTree->GetLeaf("Pileup_nTrueInt");
    preFireL = inputTree->GetLeaf("L1PreFiringWeight_Nom");
    preFireUpL = inputTree->GetLeaf("L1PreFiringWeight_Up");
    preFireDownL = inputTree->GetLeaf("L1PreFiringWeight_Dn");

    //Electron related stuff
    elePtL = inputTree->GetLeaf("Electron_pt");
    eleECorrL = inputTree->GetLeaf("Electron_eCorr");
    eleScaleUpL = inputTree->GetLeaf("Electron_dEscaleUp");
    eleScaleDownL = inputTree->GetLeaf("Electron_dEscaleDown");
    eleSigmaUpL = inputTree->GetLeaf("Electron_dEsigmaUp");
    eleSigmaDownL = inputTree->GetLeaf("Electron_dEsigmaDown");
    eleMassL = inputTree->GetLeaf("Electron_mass");
    eleEtaL = inputTree->GetLeaf("Electron_eta");
    elePhiL = inputTree->GetLeaf("Electron_phi");
    eleIso03L = inputTree->GetLeaf("Electron_pfRelIso03_all");
    eleDxyL = inputTree->GetLeaf("Electron_dxy");
    eleDzL = inputTree->GetLeaf("Electron_dz");
    eleRelJetIsoL = inputTree->GetLeaf("Electron_jetRelIso");
    eleMiniIsoL = inputTree->GetLeaf("Electron_miniPFRelIso_all");
    eleChargeL = inputTree->GetLeaf("Electron_charge");
    eleCutIDL = inputTree->GetLeaf("Electron_cutBased");
    eleMVAIDLooseL = inputTree->GetLeaf("Electron_mvaFall17V2Iso_WPL");
    eleMVAIDMediumL = inputTree->GetLeaf("Electron_mvaFall17V2Iso_WP90");
    eleMVAIDTightL = inputTree->GetLeaf("Electron_mvaFall17V2Iso_WP80");
    eleConvVetoL = inputTree->GetLeaf("Electron_convVeto");

    //Muon related stuff
    muPtL = inputTree->GetLeaf("Muon_pt");
    muEtaL = inputTree->GetLeaf("Muon_eta");
    muPhiL = inputTree->GetLeaf("Muon_phi");
    muIso03L = inputTree->GetLeaf("Muon_pfRelIso03_all");
    muIso04L = inputTree->GetLeaf("Muon_pfRelIso04_all");
    muMiniIsoL = inputTree->GetLeaf("Muon_miniPFRelIso_all");
    muChargeL = inputTree->GetLeaf("Muon_charge");
    muCutIDLooseL = inputTree->GetLeaf("Muon_looseId");
    muCutIDMediumL = inputTree->GetLeaf("Muon_mediumId");
    muCutIDTightL = inputTree->GetLeaf("Muon_tightId");
    muMVAIDL = inputTree->GetLeaf("Muon_mvaId");
    muDxyL = inputTree->GetLeaf("Muon_dxy");
    muDzL = inputTree->GetLeaf("Muon_dz");
    muRelJetIsoL = inputTree->GetLeaf("Muon_jetRelIso");
    muNTrackerLayersL = inputTree->GetLeaf("Muon_nTrackerLayers");

    //Jet related stuff
    rhoL = inputTree->GetLeaf("fixedGridRhoFastjetAll");

    metPtL = inputTree->GetLeaf("MET_pt");
    metPhiL = inputTree->GetLeaf("MET_phi");
    metDeltaUnClustXL = inputTree->GetLeaf("MET_MetUnclustEnUpDeltaX");
    metDeltaUnClustYL = inputTree->GetLeaf("MET_MetUnclustEnUpDeltaY");

    jetPtL = inputTree->GetLeaf("Jet_pt");
    jetEtaL = inputTree->GetLeaf("Jet_eta");
    jetPhiL = inputTree->GetLeaf("Jet_phi");
    jetMassL = inputTree->GetLeaf("Jet_mass");
    jetAreaL = inputTree->GetLeaf("Jet_area");
    jetDeepCSVL = inputTree->GetLeaf("Jet_btagDeepB");
    jetDeepJetL = inputTree->GetLeaf("Jet_btagDeepFlavB");
    jetPartFlavL = inputTree->GetLeaf("Jet_partonFlavour");
    jetRawFacL = inputTree->GetLeaf("Jet_rawFactor");
    jetIDL = inputTree->GetLeaf("Jet_jetId");
    jetPUIDL = inputTree->GetLeaf("Jet_puId");

    genJetPtL = inputTree->GetLeaf("GenJet_pt");
    genJetEtaL = inputTree->GetLeaf("GenJet_eta");
    genJetPhiL = inputTree->GetLeaf("GenJet_phi");

    fatJetPtL = inputTree->GetLeaf("FatJet_pt");
    fatJetEtaL = inputTree->GetLeaf("FatJet_eta");
    fatJetPhiL = inputTree->GetLeaf("FatJet_phi");
    fatJetMassL = inputTree->GetLeaf("FatJet_mass");
    fatJetAreaL = inputTree->GetLeaf("FatJet_area");
    fatJetTau1L = inputTree->GetLeaf("FatJet_tau1");
    fatJetTau2L = inputTree->GetLeaf("FatJet_tau2");
    fatJetTau3L = inputTree->GetLeaf("FatJet_tau3");
    fatJetDAK8HiggsL = inputTree->GetLeaf("FatJet_deepTag_H");
    fatJetDAK8QCDL = inputTree->GetLeaf("FatJet_deepTag_QCD");
    fatJetDAK8TvsQCDL = inputTree->GetLeaf("FatJet_deepTag_TvsQCD");
    fatJetDAK8ZvsQCDL = inputTree->GetLeaf("FatJet_deepTag_ZvsQCD");
    fatJetDAK8WvsQCDL = inputTree->GetLeaf("FatJet_deepTag_WvsQCD");
    fatJetRawFacL = inputTree->GetLeaf("FatJet_rawFactor");

    genFatJetPtL = inputTree->GetLeaf("GenJetAK8_pt");
    genFatJetEtaL = inputTree->GetLeaf("GenJetAK8_eta");
    genFatJetPhiL = inputTree->GetLeaf("GenJetAK8_phi");

    //Iso. track related
    isotrkPtL = inputTree->GetLeaf("IsoTrack_pt");
    isotrkEtaL = inputTree->GetLeaf("IsoTrack_eta");
    isotrkPhiL = inputTree->GetLeaf("IsoTrack_phi"); 
    isotrkDxyL = inputTree->GetLeaf("IsoTrack_dxy"); 
    isotrkDzL = inputTree->GetLeaf("IsoTrack_dz"); 
    isotrkPDGL = inputTree->GetLeaf("IsoTrack_pdgId"); 
    isotrkIso03L = inputTree->GetLeaf("IsoTrack_pfRelIso03_all");
    isotrkIso04L = inputTree->GetLeaf("IsoTrack_pfRelIso04_all");
    isotrkMiniIsoL = inputTree->GetLeaf("IsoTrack_miniPFRelIso_all");

    //Misc related
    evNrL = inputTree->GetLeaf("event");
    nPartonL = inputTree->GetLeaf("LHE_Njets");

    //Gen part related
    genPDGL = inputTree->GetLeaf("GenPart_pdgId");
    genMotherIdxL = inputTree->GetLeaf("GenPart_genPartIdxMother");
    genPtL = inputTree->GetLeaf("GenPart_pt");
    genPhiL = inputTree->GetLeaf("GenPart_phi");
    genEtaL = inputTree->GetLeaf("GenPart_eta");
    genMassL = inputTree->GetLeaf("GenPart_mass");
}

void NanoInput::SetWeight(){
    if(pdfWeightL != nullptr){
        pdfWeightL->GetBranch()->GetEntry(entry);
        scaleWeightL->GetBranch()->GetEntry(entry);
    }

    nTrueIntL->GetBranch()->GetEntry(entry);
    
    if(preFireL != nullptr){
        preFireL->GetBranch()->GetEntry(entry);
        preFireUpL->GetBranch()->GetEntry(entry);
        preFireDownL->GetBranch()->GetEntry(entry);
    }
}

void NanoInput::GetWeightEntry(){
    if(pdfWeightL != nullptr){
        for(int i = 0; i < 102; ++i){
            pdfWeight[i] = pdfWeightL->GetValue(i);
        }

        for(int i = 0; i < 8; ++i){
            scaleWeight[i] = scaleWeightL->GetValue(i);
        }
    }

    else{
        std::fill_n(pdfWeight, 102, 1);
        std::fill_n(scaleWeight, 8, 1);
    }
    
    if(preFireL != nullptr){
        preFire = preFireL->GetValue();
        preFireUp = preFireUpL->GetValue();
        preFireDown = preFireDownL->GetValue();
    }
    
    else{
        preFire = 1.;
        preFireUp = 1.;
        preFireDown = 1.;
    }

    nTrueInt = nTrueIntL->GetValue();
}

void NanoInput::SetTrigger(const std::vector<std::string>& names, const bool& isMETFilter){
    if(!isMETFilter){
        for(const std::string& name : names){
            triggerL.push_back(inputTree->GetLeaf(name.c_str()));
            triggers.push_back(true);
        }
    }

    else{
        for(const std::string& name : names){
            TLeaf* filter = inputTree->GetLeaf(name.c_str());
            if(filter == nullptr){
                std::cout << "MET filter not found: '" + name + "', continue without it!" << std::endl;
                continue;
            }
        
            METFilterL.push_back(filter);
            METFilter.push_back(true);
        }
    }
}

void NanoInput::ReadTrigger(){
    for(TLeaf* trigger : triggerL) trigger->GetBranch()->GetEntry(entry);
}

void NanoInput::ReadMETFilter(){
    for(TLeaf* trigger : METFilterL) trigger->GetBranch()->GetEntry(entry);
}

void NanoInput::GetTrigger(){
    for(std::size_t i = 0; i < triggerL.size(); ++i){
        triggers[i] = triggerL[i]->GetValue();
    }
}

void NanoInput::GetMETFilter(){
    for(std::size_t i = 0; i < METFilterL.size(); ++i){
        METFilter[i] = METFilterL[i]->GetValue();
    }
}

void NanoInput::ReadEleEntry(){
    if(elePtL->GetBranch()->GetReadEntry() == entry) return;

    elePtL->GetBranch()->GetEntry(entry);
    eleMassL->GetBranch()->GetEntry(entry);
    elePhiL->GetBranch()->GetEntry(entry);
    eleEtaL->GetBranch()->GetEntry(entry);
    eleIso03L->GetBranch()->GetEntry(entry);
    eleMiniIsoL->GetBranch()->GetEntry(entry);
    eleChargeL->GetBranch()->GetEntry(entry);
    eleDxyL->GetBranch()->GetEntry(entry);
    eleDzL->GetBranch()->GetEntry(entry);
    eleRelJetIsoL->GetBranch()->GetEntry(entry);
    eleCutIDL->GetBranch()->GetEntry(entry);
    eleMVAIDLooseL->GetBranch()->GetEntry(entry);
    eleMVAIDMediumL->GetBranch()->GetEntry(entry);
    eleMVAIDTightL->GetBranch()->GetEntry(entry);
    eleConvVetoL->GetBranch()->GetEntry(entry);

    if(eleScaleUpL){
        eleScaleUpL->GetBranch()->GetEntry(entry);
        eleScaleDownL->GetBranch()->GetEntry(entry);
        eleSigmaUpL->GetBranch()->GetEntry(entry);
        eleSigmaDownL->GetBranch()->GetEntry(entry);
    }

    if(eleECorrL) eleECorrL->GetBranch()->GetEntry(entry);

    eleSize = elePtL->GetLen();
}

void NanoInput::GetElectron(const std::size_t& idx){
    elePt = elePtL->GetValue(idx);
    elePhi = elePhiL->GetValue(idx);
    eleEta = eleEtaL->GetValue(idx);
    eleIso03 = eleIso03L->GetValue(idx);
    eleMiniIso = eleMiniIsoL->GetValue(idx);
    eleDxy = eleDxyL->GetValue(idx);
    eleDz = eleDzL->GetValue(idx);
    eleRelJetIso = eleRelJetIsoL->GetValue(idx);
    eleConvVeto = eleConvVetoL->GetValue(idx);

    //Apply shift up/down for energy smearing/scaling
    if(eleScaleUpL and eleECorrL){
        float eCorr = eleECorrL->GetValue(idx);

        eleP4.SetPt(elePt/eCorr);
        eleP4.SetEta(eleEta);
        eleP4.SetPhi(elePhi);
        eleP4.SetM(eleMassL->GetValue(idx));

        elePtScaleUp = (eCorr + eleScaleUpL->GetValue(idx)) * eleP4.Pt();
        elePtScaleDown = (eCorr - eleScaleDownL->GetValue(idx)) * eleP4.Pt();
        elePtSigmaUp = (eCorr + eleSigmaUpL->GetValue(idx))* eleP4.Pt();
        elePtSigmaDown = (eCorr - eleSigmaDownL->GetValue(idx)) * eleP4.Pt();
    }

    eleCharge = eleChargeL->GetValue(idx);
    eleCutID = eleCutIDL->GetValue(idx) - 1;
    eleMVAID = eleMVAIDLooseL->GetValue(idx) + eleMVAIDMediumL->GetValue(idx) + eleMVAIDTightL->GetValue(idx);
}

void NanoInput::ReadMuEntry(){
    if(muPtL->GetBranch()->GetReadEntry() == entry) return;

    muRandomNumber = gRandom->Rndm();

    muPtL->GetBranch()->GetEntry(entry);
    muEtaL->GetBranch()->GetEntry(entry);
    muPhiL->GetBranch()->GetEntry(entry);
    muIso03L->GetBranch()->GetEntry(entry);
    muIso04L->GetBranch()->GetEntry(entry);
    muMiniIsoL->GetBranch()->GetEntry(entry);
    muChargeL->GetBranch()->GetEntry(entry);
    muDxyL->GetBranch()->GetEntry(entry);
    muDzL->GetBranch()->GetEntry(entry);
    muRelJetIsoL->GetBranch()->GetEntry(entry);
    muCutIDLooseL->GetBranch()->GetEntry(entry);
    muCutIDMediumL->GetBranch()->GetEntry(entry);
    muCutIDTightL->GetBranch()->GetEntry(entry);
    muMVAIDL->GetBranch()->GetEntry(entry);
    muNTrackerLayersL->GetBranch()->GetEntry(entry);

    muSize = muPtL->GetLen();
}

void NanoInput::GetMuon(const std::size_t& idx){
    muPt = muPtL->GetValue(idx);
    muPhi = muPhiL->GetValue(idx);
    muEta = muEtaL->GetValue(idx);
    muIso03 = muIso03L->GetValue(idx);
    muIso04 = muIso04L->GetValue(idx);
    muMiniIso = muMiniIsoL->GetValue(idx);
    muDxy = muDxyL->GetValue(idx);
    muDz = muDzL->GetValue(idx);
    muRelJetIso = muRelJetIsoL->GetValue(idx);
    muNTrackerLayers = muNTrackerLayersL->GetValue(idx);

    muCharge = muChargeL->GetValue(idx);
    muMVAID = muMVAIDL->GetValue(idx);

    if(muCutIDTightL->GetValue(idx)) muCutID = 3;
    else if(muCutIDMediumL->GetValue(idx)) muCutID = 2;
    else if(muCutIDLooseL->GetValue(idx)) muCutID = 1;
    else muCutID = 0;
}

void NanoInput::ReadJetEntry(const bool& isData){
    if(jetPtL->GetBranch()->GetReadEntry() == entry) return;

    rhoL->GetBranch()->GetEntry(entry);
    rho = rhoL->GetValue();

    metPhiL->GetBranch()->GetEntry(entry);
    metPtL->GetBranch()->GetEntry(entry);
    metDeltaUnClustXL->GetBranch()->GetEntry(entry);
    metDeltaUnClustYL->GetBranch()->GetEntry(entry);
    
    metPt = metPtL->GetValue();
    metPhi = metPhiL->GetValue();
    metDeltaUnClustX = metDeltaUnClustXL->GetValue();
    metDeltaUnClustY = metDeltaUnClustYL->GetValue();

    jetPtL->GetBranch()->GetEntry(entry);
    jetEtaL->GetBranch()->GetEntry(entry);
    jetPhiL->GetBranch()->GetEntry(entry);
    jetMassL->GetBranch()->GetEntry(entry);
    jetAreaL->GetBranch()->GetEntry(entry);
    jetDeepJetL->GetBranch()->GetEntry(entry);
    jetDeepCSVL->GetBranch()->GetEntry(entry);
    jetRawFacL->GetBranch()->GetEntry(entry);
    jetIDL->GetBranch()->GetEntry(entry);
    jetPUIDL->GetBranch()->GetEntry(entry);

    fatJetPtL->GetBranch()->GetEntry(entry);
    fatJetEtaL->GetBranch()->GetEntry(entry);
    fatJetPhiL->GetBranch()->GetEntry(entry);
    fatJetMassL->GetBranch()->GetEntry(entry);
    fatJetAreaL->GetBranch()->GetEntry(entry);
    fatJetTau1L->GetBranch()->GetEntry(entry);
    fatJetTau2L->GetBranch()->GetEntry(entry);
    fatJetTau3L->GetBranch()->GetEntry(entry);
    fatJetDAK8HiggsL->GetBranch()->GetEntry(entry);
    fatJetDAK8QCDL->GetBranch()->GetEntry(entry);
    fatJetDAK8TvsQCDL->GetBranch()->GetEntry(entry);
    fatJetDAK8ZvsQCDL->GetBranch()->GetEntry(entry);
    fatJetDAK8WvsQCDL->GetBranch()->GetEntry(entry);
    fatJetRawFacL->GetBranch()->GetEntry(entry);

    jetSize = jetPtL->GetLen();
    fatJetSize = fatJetPtL->GetLen();

    if(!isData){
        jetPartFlavL->GetBranch()->GetEntry(entry);

        genJetPtL->GetBranch()->GetEntry(entry);
        genJetEtaL->GetBranch()->GetEntry(entry);
        genJetPhiL->GetBranch()->GetEntry(entry);

        genFatJetPtL->GetBranch()->GetEntry(entry);
        genFatJetEtaL->GetBranch()->GetEntry(entry);
        genFatJetPhiL->GetBranch()->GetEntry(entry);

        genJetSize = genJetPtL->GetLen();
        genFatJetSize = genFatJetPtL->GetLen();
    }
}

void NanoInput::GetJet(const std::size_t& idx){
    jetPt = jetPtL->GetValue(idx);
    jetPtRaw = jetPt*(1. - jetRawFacL->GetValue(idx));
    jetMass = jetMassL->GetValue(idx);
    jetMassRaw = jetMass*(1. - jetRawFacL->GetValue(idx));
    jetPhi = jetPhiL->GetValue(idx);
    jetEta = jetEtaL->GetValue(idx);
    jetDeepJet = jetDeepJetL->GetValue(idx);
    jetDeepCSV = jetDeepCSVL->GetValue(idx);
    jetArea = jetAreaL->GetValue(idx);
    jetID = jetIDL->GetValue(idx);
    jetPUID = jetPUIDL->GetValue(idx);
    
    if(jetPartFlavL != nullptr) jetPartFlav = jetPartFlavL->GetValue(idx);
}

void NanoInput::GetFatJet(const std::size_t& idx){
    fatJetPt = fatJetPtL->GetValue(idx);
    fatJetPtRaw = fatJetPt*(1. - fatJetRawFacL->GetValue(idx));
    fatJetMass = fatJetMassL->GetValue(idx);
    fatJetMassRaw = fatJetMass*(1. - fatJetRawFacL->GetValue(idx));
    fatJetPhi = fatJetPhiL->GetValue(idx);
    fatJetEta = fatJetEtaL->GetValue(idx);
    fatJetArea = fatJetAreaL->GetValue(idx);
    fatJetTau1 = fatJetTau1L->GetValue(idx);
    fatJetTau2 = fatJetTau2L->GetValue(idx);
    fatJetTau3 = fatJetTau3L->GetValue(idx);

    //Demangle binary scores to get raw scores of DeepAK8
    float DAK8Higgs = fatJetDAK8HiggsL->GetValue(idx);
    float DAK8QCD = fatJetDAK8QCDL->GetValue(idx);
    float DAK8TvsQCD = fatJetDAK8TvsQCDL->GetValue(idx);
    float DAK8ZvsQCD = fatJetDAK8ZvsQCDL->GetValue(idx);
    float DAK8WvsQCD = fatJetDAK8WvsQCDL->GetValue(idx);

    std::vector<float> DK8Scores = {DAK8QCD, 
                                    demangleDK8(DAK8WvsQCD, DAK8QCD), 
                                    demangleDK8(DAK8TvsQCD, DAK8QCD),
                                    demangleDK8(DAK8ZvsQCD, DAK8QCD),
                                    DAK8Higgs};

    fatJetDAK8ID = std::max_element(DK8Scores.begin(), DK8Scores.end()) - DK8Scores.begin();
}

void NanoInput::GetGenJet(const std::size_t& idx){
    if(genJetPtL->GetBranch()->GetReadEntry() == entry) return;

    genJetPt = genJetPtL->GetValue(idx);
    genJetEta = genJetEtaL->GetValue(idx);
    genJetPhi = genJetPhiL->GetValue(idx);
}

void NanoInput::GetGenFatJet(const std::size_t& idx){
    genFatJetPt = genFatJetPtL->GetValue(idx);
    genFatJetEta = genFatJetEtaL->GetValue(idx);
    genFatJetPhi = genFatJetPhiL->GetValue(idx);
}

void NanoInput::ReadIsotrkEntry(){
    if(isotrkPtL->GetBranch()->GetReadEntry() == entry) return;

    isotrkPtL->GetBranch()->GetEntry(entry);
    isotrkEtaL->GetBranch()->GetEntry(entry);
    isotrkPhiL->GetBranch()->GetEntry(entry);
    isotrkIso03L->GetBranch()->GetEntry(entry);
    isotrkDxyL->GetBranch()->GetEntry(entry);
    isotrkDzL->GetBranch()->GetEntry(entry);
    isotrkPDGL->GetBranch()->GetEntry(entry);
    isotrkMiniIsoL->GetBranch()->GetEntry(entry);

    isotrkSize = isotrkPtL->GetLen();
}

void NanoInput::GetIsotrk(const std::size_t& idx){
    isotrkPt = isotrkPtL->GetValue(idx);
    isotrkEta = isotrkEtaL->GetValue(idx);
    isotrkPhi = isotrkPhiL->GetValue(idx);
    isotrkIso03 = isotrkIso03L->GetValue(idx);
    isotrkMiniIso = isotrkMiniIsoL->GetValue(idx);
    isotrkDxy = isotrkDxyL->GetValue(idx);
    isotrkDz = isotrkDzL->GetValue(idx);
    isotrkPDG = isotrkPDGL->GetValue(idx);
}

void NanoInput::ReadMiscEntry(const bool& isData){
    evNrL->GetBranch()->GetEntry(entry);
    if(nPartonL != nullptr) nPartonL->GetBranch()->GetEntry(entry);
}

void NanoInput::GetMisc(){
    evNr = evNrL->GetValue();
    if(nPartonL != nullptr) nParton = nPartonL->GetTypedValue<short>();
}

void NanoInput::ReadGenEntry(){
    if(genPDGL->GetBranch()->GetReadEntry() == entry) return;

    alreadyMatchedIdx.clear();

    genPDGL->GetBranch()->GetEntry(entry);
    genMotherIdxL->GetBranch()->GetEntry(entry);
    genPtL->GetBranch()->GetEntry(entry);
    genPhiL->GetBranch()->GetEntry(entry);
    genEtaL->GetBranch()->GetEntry(entry);
    genMassL->GetBranch()->GetEntry(entry);

    genSize = genPtL->GetLen();
}

void NanoInput::GetGenPart(const std::size_t& idx){
    genPDG = genPDGL->GetValue(idx);
    genMotherIdx = genMotherIdxL->GetValue(idx);
    genPt = genPtL->GetValue(idx);
    genPhi = genPhiL->GetValue(idx);
    genEta = genEtaL->GetValue(idx);
    genMass = genMassL->GetValue(idx);
}
