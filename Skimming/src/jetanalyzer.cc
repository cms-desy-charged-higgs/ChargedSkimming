#include <ChargedSkimming/Skimming/interface/jetanalyzer.h>

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, std::vector<jToken>& jetTokens, std::vector<genjToken>& genjetTokens, mToken &metToken, edm::EDGetTokenT<double> &rhoToken, genPartToken& genParticleToken, secvtxToken& vertexToken, const std::string& systematic):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    jetTokens(jetTokens),
    metToken(metToken),
    rhoToken(rhoToken),
    genParticleToken(genParticleToken),
    vertexToken(vertexToken)
    {
        if(systematic.find("JEC") != std::string::npos){
            jecSyst = systematic;
            jecSyst.replace(jecSyst.find("JEC"), 3, "");

            if(jecSyst.find("Up") != std::string::npos){
                isUp = true;
                jecSyst.replace(jecSyst.find("Up"), 2, "");    
            }

            else{
                isUp = false;
                jecSyst.replace(jecSyst.find("Down"), 4, "");    
            }
        }

        if(systematic.find("JER") != std::string::npos){
            isJERsyst = true;

            if(systematic.find("Up") != std::string::npos){
                isUp = true;
            }
        }
    }


void JetAnalyzer::SetCorrector(const JetType &type, const int& runNumber){
    std::vector<JetCorrectorParameters> corrVec;

    for(std::string fileName: isData? JECDATA[era] : JECMC[era]){
        if(fileName.find("@") != std::string::npos){
            for(std::pair<std::string, std::pair<int, int>> eraNames: runEras[era]){
                if(eraNames.second.first <= runNumber and runNumber <= eraNames.second.second){
                    fileName.replace(fileName.find("@"), 1, eraNames.first);
                }
            }
        }

        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        corrVec.push_back(JetCorrectorParameters(fileName));
    }
  
    jetCorrector[type] = new FactorizedJetCorrector(corrVec);
}

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetAnalyzer::CorrectEnergy(const TLorentzVector &jet, const float &rho, const float &area, const JetType &type){
    jetCorrector[type]->setJetPt(jet.Pt());
    jetCorrector[type]->setJetEta(jet.Eta());
    jetCorrector[type]->setRho(rho);
    jetCorrector[type]->setJetA(area);

    double correction = jetCorrector[type]->getCorrection();
    
    return correction;
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
float JetAnalyzer::SmearEnergy(const TLorentzVector &jet, const float &rho, const float &coneSize, const JetType &type, const std::vector<reco::GenJet> &genJets){
    //Configure jet SF reader
    jetParameter.setJetPt(jet.Pt()).setJetEta(jet.Eta()).setRho(rho);

    float reso = resolution[type].getResolution(jetParameter);
    float resoSF; 
    if(!isJERsyst) resoSF = resolution_sf[type].getScaleFactor(jetParameter);
    else resoSF = resolution_sf[type].getScaleFactor(jetParameter, isUp ? Variation::UP : Variation::DOWN);
        
    float smearFac = 1.; 
    float dR;
    float genPt, genPhi, genEta, genMass;
    unsigned int size;

    if(isNANO) size = (type == AK4) ? genJetPt->GetSize(): genFatJetPt->GetSize();
    else size = genJets.size();

    //Loop over all gen jets and find match
    for(unsigned int i = 0; i < size; i++){
        if(isNANO){
            genPt = (type == AK4) ? genJetPt->At(i): genFatJetPt->At(i);
            genPhi = (type == AK4) ? genJetPhi->At(i): genFatJetPhi->At(i);
            genEta = (type == AK4) ? genJetEta->At(i): genFatJetEta->At(i);
            genMass = (type == AK4) ? genJetMass->At(i): genFatJetMass->At(i);
        }

        else{
            genPt = genJets.at(i).pt();
            genPhi = genJets.at(i).phi();
            genEta = genJets.at(i).eta();
            genMass = genJets.at(i).mass();
        }

        dR = std::sqrt(std::pow(jet.Phi()-genPhi, 2) + std::pow(jet.Eta()-genEta, 2));

        //Check if jet and gen jet are matched
        if(dR < coneSize/2. and abs(jet.Pt() - genPt) < 3.*reso*jet.Pt()){
            TLorentzVector gJet;
            gJet.SetPtEtaPhiM(genPt, genEta, genPhi, genMass);
            genJet[type] = gJet;
            break;
        }
    }  

    //If you found gen matched 
    if(genJet[type] != TLorentzVector()){
        smearFac = 1.+(resoSF-1)*(jet.Pt() - genJet[type].Pt())/jet.Pt(); 
    }

    //If no match, smear with gaussian pdf
    else if(resoSF > 1.){
        std::default_random_engine generator;
        std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));

        smearFac = 1. + gaus(generator);
    }


    //Check if direction of jet not changed
    if(jet.Pt()*smearFac < 1e-2){
        smearFac = 1e-2/jet.Pt();
    }

    return smearFac;
}

int JetAnalyzer::SetGenParticles(TLorentzVector& validJet, const int &i, const int &pdgID, const JetType &type, const std::vector<reco::GenParticle>& genParticle){
    int nParton=0;
    bool isFromh1 = true;
    bool isFromh2 = true;

    //Check if gen matched particle exist
    if(genJet[type].Pt() != 0){
        float dR;
        
        //Find Gen particle to gen Jet
        int size=isNANO ? genPt->GetSize() : genParticle.size();

        for(int i=0; i < size; i++){
            const reco::Candidate* parton=NULL;
            int index=0;

            int ID = isNANO ? abs(genID->At(i)) : abs(genParticle.at(i).pdgId());
    
            if(ID == pdgID){
                if(isNANO) index = FirstCopy(i, pdgID);
                else parton = FirstCopy(genParticle.at(i), pdgID);
            }
    
            else continue;

            float eta, phi;
            phi = isNANO ? genPhi->At(index) : parton->phi();
            eta = isNANO ? genEta->At(index) : parton->eta();

            dR = std::sqrt(std::pow(phi - genJet[type].Phi(), 2) + std::pow(eta-genJet[type].Eta(), 2));
            float rMin = type == AK4 ? 0.3 : 0.4;

            if(dR <  rMin){
                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(parton->mother()->pdgId());

                if(motherID == 25){
                    const reco::Candidate* hBoson=NULL;

                    nParton++;
                                
                    if(isNANO) index = FirstCopy(genMotherIdx->At(index), 25);
                    else hBoson = FirstCopy(parton->mother(), 25);

                    int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(hBoson->mother()->pdgId());  

                    if(motherID == 37){
                        isFromh1 = isFromh1 && true;
                        isFromh2 = isFromh2 && true;
                    }

                    else{
                        isFromh1 = isFromh1 && false;
                        isFromh2 = isFromh2 && true;
                    }
                }
            }

            if(nParton==1 and type==AK4){
                return !isFromh1 and !isFromh2 ? -1  : isFromh1 ? 1 : 2;
            }

            if(nParton==2 and type==AK8){
                return !isFromh1 and !isFromh2 ? -1  : isFromh1 ? 1 : 2;
            }
        }
    }

    return -1.;
}

void JetAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    JECMC = {
                {2017, {filePath + "/JEC/Fall17_17Nov2017_V32_MC_L1FastJet_&PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L2Relative_&PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L3Absolute_&PFchs.txt"}},
                            
    };

    JECDATA = {
                {2017, {filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L1FastJet_&PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2Relative_&PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L3Absolute_&PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2L3Residual_&PFchs.txt"}},
             
    };

    JECUNC = {
                {2017, filePath + "/JEC/Fall17_17Nov2017_V32_MC_UncertaintySources_&PFchs.txt"},
    };

    JMESF = {
                {2017, filePath + "/JME/Fall17_V3_MC_SF_&PFchs.txt"},     
    };

    JMEPtReso = {
                {2017, filePath + "/JME/Fall17_V3_MC_PtResolution_&PFchs.txt"},    
    };


    bTagSF = {
            {AK4, {
                {2017, filePath + "/btagSF/DeepFlavour_94XSF_V1_B_F.csv"},
                }
            },

            {AK8, {
                 {2017, filePath + "/btagSF/subjet_DeepCSV_94XSF_V4_B_F.csv"},
                }
            },
    };


    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94
    bTagCuts = {
            {AK4, {
                {2017, {0.051, 0.3033, 0.7489}}, //Loose, Medium, Tight
                }
            },

            {AK8, {
                {2017, {0.1522, 0.4941}}, //Loose, Medium, Tight
                }
            },
    };

    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    if(isNANO){
        //Initiliaze TTreeReaderValues
        fatJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_pt");
        fatJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_eta");
        fatJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_phi");
        fatJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_mass");
        fatJetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_area");
        fatJetCSV = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_btagDeepB");
        fatJetTau1 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau1");
        fatJetTau2 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau2");
        fatJetTau3 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau3");

        jetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_pt");
        jetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_eta");
        jetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_phi");
        jetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_mass");
        jetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_area");
        jetRho = std::make_unique<TTreeReaderValue<float>>(*reader, "fixedGridRhoFastjetAll");
        jetDeepBValue = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_btagDeepFlavB");
        valueHT = std::make_unique<TTreeReaderValue<float>>(*reader, "SoftActivityJetHT");
        
        metPhi = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_phi");
        metPt = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_pt");

        if(!this->isData){
            genJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_pt");
            genJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_eta");
            genJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_phi");
            genJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_mass");

            genFatJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_pt");
            genFatJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_eta");
            genFatJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_phi");
            genFatJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_mass");
            
            jetGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_genJetIdx");
        }

        //Set TTreeReader for genpart and trigger obj from baseanalyzer
        SetCollection(this->isData);
    }
    
    for(JetType type: {AK4, AK8}){
        //Set configuration for bTagSF reader  ##https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
        calib[type] = BTagCalibration(std::to_string(type), bTagSF[type][era]);

        looseReader[type] = BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
        mediumReader[type] = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});  
        tightReader[type] = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central", {"up", "down"});  

        looseReader[type].load(calib[type],  BTagEntry::FLAV_B, (type==AK4) ? "comb": "lt");  
        mediumReader[type].load(calib[type],  BTagEntry::FLAV_B, (type==AK4) ? "comb": "lt");  
        
        if(type!=AK8){
            tightReader[type].load(calib[type],  BTagEntry::FLAV_B, (type==AK4) ? "comb": "lt");  
        }
    
        //Set configuration for JER tools
        std::string fileName = JMEPtReso[era];
        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        resolution[type] = JME::JetResolution(fileName);

        fileName = JMESF[era];
        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        resolution_sf[type] = JME::JetResolutionScaleFactor(fileName);

        //Set object to get JEC uncertainty
        if(jecSyst != ""){
            std::string fileName = JECUNC[era];
            fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");  

            jecUnc[type] = new JetCorrectionUncertainty(JetCorrectorParameters(fileName, "Total"));
        }

        else jecUnc[type] = NULL;
    }

    //Set output names
    JetfloatNames = {"E", "Px", "Py", "Pz", "FatJetIdx", "isFromh", "loosebTagSF",  "mediumbTagSF", "tightbTagSF"};
    FatJetfloatNames = {"E", "Px", "Py", "Pz", "oneSubJettiness", "twoSubJettiness", "threeSubJettiness", "isFromh", "loosebTagSF", "mediumbTagSF"};

    if(!isSyst){
        std::vector<std::string> SFJetvariations = {"loosebTagSFUp", "loosebTagDown", "mediumbTagSFUp", "mediumbTagSFDown", "tightbTagSFUp", "tightbTagSFDown"};
        std::vector<std::string> SFFatJetvariations = {"loosebTagSFUp", "loosebTagDown", "mediumbTagSFUp", "mediumbTagSFDown"};

        JetfloatNames.insert(JetfloatNames.end(), SFJetvariations.begin(), SFJetvariations.end()); 
        FatJetfloatNames.insert(FatJetfloatNames.end(), SFFatJetvariations.begin(), SFFatJetvariations.end()); 
    }

    JetParticlefloatNames = {"E", "Px", "Py", "Pz", "Vx", "Vy", "Vz", "Charge", "FatJetIdx"};
    boolNames = {"isLooseB", "isMediumB", "isTightB"};

    JetfloatVariables = std::vector<std::vector<float>>(JetfloatNames.size(), std::vector<float>());
    FatJetfloatVariables = std::vector<std::vector<float>>(FatJetfloatNames.size(), std::vector<float>());
    JetParticlefloatVariables = std::vector<std::vector<float>>(JetParticlefloatNames.size(), std::vector<float>());
    VertexfloatVariables = std::vector<std::vector<float>>(JetParticlefloatNames.size(), std::vector<float>());

    JetboolVariables = std::vector<std::vector<bool>>(boolNames.size(), std::vector<bool>());
    FatJetboolVariables = std::vector<std::vector<bool>>(boolNames.size()-1, std::vector<bool>());

    //Set Branches of output tree
    for(TTree* tree: trees){
        for(unsigned int i=0; i<JetfloatVariables.size(); i++){
            tree->Branch(("Jet_" + JetfloatNames[i]).c_str(), &JetfloatVariables[i]);
        }

        for(unsigned int i=0; i<FatJetfloatVariables.size(); i++){
            tree->Branch(("FatJet_" + FatJetfloatNames[i]).c_str(), &FatJetfloatVariables[i]);
        }

        for(unsigned int i=0; i<JetParticlefloatVariables.size(); i++){
            tree->Branch(("JetParticle_" + JetParticlefloatNames[i]).c_str(), &JetParticlefloatVariables[i]);
        }

        for(unsigned int i=0; i<VertexfloatVariables.size(); i++){
            tree->Branch(("SecondaryVertex_" + JetParticlefloatNames[i]).c_str(), &VertexfloatVariables[i]);
        }

        for(unsigned int i=0; i<boolNames.size(); i++){
            tree->Branch(("Jet_" + boolNames[i]).c_str(), &JetboolVariables[i]);
            if(i!=0){
                tree->Branch(("FatJet_" + boolNames[i]).c_str(), &JetboolVariables[i-1]);
            }
        }

        tree->Branch("MET_Px", &metPx);
        tree->Branch("MET_Py", &metPy);
        tree->Branch("HT", &HT);
    }
}

void JetAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::vector<float>& variable: JetfloatVariables){
        variable.clear();
    }

    for(std::vector<float>& variable: FatJetfloatVariables){
        variable.clear();
    }

    for(std::vector<float>& variable: JetParticlefloatVariables){
        variable.clear();
    }

    for(std::vector<float>& variable: VertexfloatVariables){
        variable.clear();
    }

    for(std::vector<bool>& variable: JetboolVariables){
        variable.clear();
    }

    for(std::vector<bool>& variable: FatJetboolVariables){
        variable.clear();
    }

    int nSubJets=0;
    HT=0;
    runNumber = isNANO ? *run->Get() : event->eventAuxiliary().id().run(); 

    if(jetCorrector.empty()){
        for(const JetType& type: {AK4, AK8}){
            SetCorrector(type, runNumber);
        }
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Jet>> jets, fatJets;
    edm::Handle<std::vector<reco::GenJet>> genJets, genfatJets;
    edm::Handle<std::vector<pat::MET>> MET;
    edm::Handle<double> rho;
    edm::Handle<std::vector<reco::GenParticle>> genParts;
    edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> secVtx;

    if(!isNANO){
        event->getByToken(jetTokens[0], jets);
        event->getByToken(jetTokens[1], fatJets);
        event->getByToken(metToken, MET);
        event->getByToken(rhoToken, rho);
        event->getByToken(vertexToken, secVtx);

        if(!isData){
            event->getByLabel(edm::InputTag("slimmedGenJets"), genJets);
            event->getByLabel(edm::InputTag("slimmedGenJetsAK8"), genfatJets);
        }
    }

    //JER smearing
    float smearFac = 1.;
    float corrFac = 1.;

    //MET values not correct for JER yet
    metPx = isNANO ? *metPt->Get()*std::cos(*metPhi->Get()) : MET->at(0).uncorPx();
    metPy = isNANO ? *metPt->Get()*std::sin(*metPhi->Get()) : MET->at(0).uncorPy();
    
    float fatJetSize = isNANO ? fatJetPt->GetSize() : fatJets->size();
    float jetSize = isNANO ? jetPt->GetSize() : jets->size();

    //Check if initial jet collection has enough jets
    unsigned int nMin = 0;

    for(CutFlow& cutflow: cutflows){
        if(jetSize < cutflow.nMinJet){
            nMin++;
            cutflow.passed = false;
        }
    }

    if(nMin == cutflows.size()) return;
        
    //Loop over all fat jets
    for(unsigned int i = 0; i < fatJetSize; i++){
        float fatPt = isNANO ? fatJetPt->At(i) : fatJets->at(i).pt();
        float fatEta = isNANO ? fatJetEta->At(i) : fatJets->at(i).eta();
        float fatPhi = isNANO ? fatJetPhi->At(i) : fatJets->at(i).phi();
        float fatMass = isNANO ? fatJetMass->At(i) : fatJets->at(i).mass();

        TLorentzVector lVec;
        lVec.SetPtEtaPhiM(fatPt, fatEta, fatPhi, fatMass);
    
        corrFac = isNANO ? CorrectEnergy(lVec,  *jetRho->Get(), fatJetArea->At(i), AK8) : CorrectEnergy(lVec, *rho, jets->at(i).jetArea(), AK8);

        //Get jet uncertainty
        if(jecUnc[AK8] !=  NULL){
            jecUnc[AK8]->setJetPt(corrFac*lVec.Pt());
            jecUnc[AK8]->setJetEta(lVec.Eta());
            float unc = jecUnc[AK8]->getUncertainty(isUp);
            corrFac *= isUp ? 1 + unc : 1 - unc;
        }

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(lVec*corrFac, *jetRho->Get(), 0.8, AK8) : SmearEnergy(lVec*corrFac, *rho, 0.8, AK8, *genfatJets);
            lVec *= smearFac*corrFac;
        }

        else lVec *= corrFac;

        if(lVec.Pt() > 170. and lVec.M() > 40. and abs(lVec.Eta()) < etaCut){
            //Fatjet four momentum components
            FatJetfloatVariables[0].push_back(lVec.E());   //Energy
            FatJetfloatVariables[1].push_back(lVec.Px());  //Px
            FatJetfloatVariables[2].push_back(lVec.Py());  //Py
            FatJetfloatVariables[3].push_back(lVec.Pz());  //Pz

            //Check for btag
            float DeepCSV = 0;

            if(isNANO){ 
                DeepCSV = jetDeepBValue->At(i);
            } 

            else{
                for(std::string disc: {"pfDeepCSVJetTags:probb", "pfDeepCSVJetTags:probbb"}){
                    DeepCSV +=  jets->at(i).bDiscriminator(disc);
                }
            }

            //Check for btag
            FatJetboolVariables[0].push_back(bTagCuts[AK8][era][0] < DeepCSV);
            FatJetboolVariables[1].push_back(bTagCuts[AK8][era][1] < DeepCSV);

            //Nsubjettiness
            FatJetfloatVariables[4].push_back(isNANO ? fatJetTau1->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1"));
            FatJetfloatVariables[5].push_back(isNANO ? fatJetTau2->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2"));
            FatJetfloatVariables[6].push_back(isNANO ? fatJetTau3->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3"));
            if(!isData){
                //btag SF
                FatJetfloatVariables[8].push_back(looseReader[AK8].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                FatJetfloatVariables[9].push_back(mediumReader[AK8].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));

                if(!isSyst){
                    FatJetfloatVariables[10].push_back(looseReader[AK8].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    FatJetfloatVariables[11].push_back(looseReader[AK8].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    FatJetfloatVariables[12].push_back(mediumReader[AK8].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    FatJetfloatVariables[13].push_back(mediumReader[AK8].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                }

                if(isNANO) FatJetfloatVariables[7].push_back(SetGenParticles(lVec, i, 5, AK8));
                else{
                    event->getByToken(genParticleToken, genParts);
                    FatJetfloatVariables[7].push_back(SetGenParticles(lVec, i, 5, AK8, *genParts));
                }
            }

            //Fill in particle flow candidates
            if(!isNANO){
                for(unsigned int j = 0; j < fatJets->at(i).numberOfDaughters(); j++){
                    reco::Candidate const * cand = fatJets->at(i).daughter(j);
                    
                    if(cand->numberOfDaughters() == 0){
                        //Jet particle four momentum components
                        JetParticlefloatVariables[0].push_back(cand->energy());   //Energy
                        JetParticlefloatVariables[1].push_back(cand->px());  //Px
                        JetParticlefloatVariables[2].push_back(cand->py());  //Py
                        JetParticlefloatVariables[3].push_back(cand->pz());  //Pz

                        //Jet particle vertex
                        JetParticlefloatVariables[4].push_back(cand->vx());      
                        JetParticlefloatVariables[5].push_back(cand->vy()); 
                        JetParticlefloatVariables[6].push_back(cand->vz());
                        JetParticlefloatVariables[7].push_back(cand->charge());
            
                        //Fat Jet Index
                        JetParticlefloatVariables[8].push_back(FatJetfloatVariables[0].size()-1);
                    }
                            
                    else{
                        for(unsigned int k = 0; k < cand->numberOfDaughters(); k++){
                            reco::Candidate const * cand2 = cand->daughter(k);

                            //Jet particle four momentum components
                            JetParticlefloatVariables[0].push_back(cand2->energy());   //Energy
                            JetParticlefloatVariables[1].push_back(cand2->px());  //Px
                            JetParticlefloatVariables[2].push_back(cand2->py());  //Py
                            JetParticlefloatVariables[3].push_back(cand2->pz());  //Pz

                            //Jet particle vertex
                            JetParticlefloatVariables[4].push_back(cand2->vx());
                            JetParticlefloatVariables[5].push_back(cand2->vy()); 
                            JetParticlefloatVariables[6].push_back(cand2->vz());
                            JetParticlefloatVariables[7].push_back(cand2->charge());
                
                            //Fat Jet Index
                            JetParticlefloatVariables[8].push_back(FatJetfloatVariables[0].size()-1);
                        }
                    }
                }

                for(const reco::VertexCompositePtrCandidate &vtx: *secVtx){
                    TLorentzVector vtxP4;
                    vtxP4.SetPtEtaPhiM(vtx.p4().Pt(), vtx.p4().Eta(), vtx.p4().Phi(), vtx.p4().M());

                    if(lVec.DeltaR(vtxP4) < 0.8){
                        //SV four momentum components
                        VertexfloatVariables[0].push_back(vtx.energy());   //Energy
                        VertexfloatVariables[1].push_back(vtx.px());  //Px
                        VertexfloatVariables[2].push_back(vtx.py());  //Py
                        VertexfloatVariables[3].push_back(vtx.pz());  //Pz

                        //SV vertex
                        VertexfloatVariables[4].push_back(vtx.vx());
                        VertexfloatVariables[5].push_back(vtx.vy()); 
                        VertexfloatVariables[6].push_back(vtx.vz());
                        VertexfloatVariables[7].push_back(vtx.charge());
            
                        //Fat Jet Index
                        VertexfloatVariables[8].push_back(FatJetfloatVariables[0].size()-1);
                    }
                }
            }
        }
    }

    //Loop over all jets
    for(unsigned int i = 0; i < jetSize; i++){
        float pt = isNANO ? jetPt->At(i) : jets->at(i).pt();
        float eta = isNANO ? jetEta->At(i) : jets->at(i).eta();
        float phi = isNANO ? jetPhi->At(i) : jets->at(i).phi();
        float mass = isNANO ? jetMass->At(i) : jets->at(i).mass();

        //Define here already jet, because of smearing of 4-vec
        TLorentzVector lVec;
        lVec.SetPtEtaPhiM(pt, eta, phi, mass);

        corrFac = isNANO ? CorrectEnergy(lVec,  *jetRho->Get(), jetArea->At(i), AK4) : CorrectEnergy(lVec, *rho, jets->at(i).jetArea(), AK4);

        //Get jet uncertainty
        if(jecUnc[AK8] !=  NULL){
            jecUnc[AK4]->setJetPt(corrFac*lVec.Pt());
            jecUnc[AK4]->setJetEta(lVec.Eta());
            float unc = jecUnc[AK4]->getUncertainty(isUp);
            corrFac *= isUp ? 1 + unc : 1 - unc;
        }

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(lVec*corrFac,  *jetRho->Get(), jetArea->At(i), AK4) : SmearEnergy(lVec, *rho, 0.4, AK4, *genJets);

            lVec*=smearFac*corrFac;
        }

        else lVec*=corrFac;

        //Correct met
        metPx+= lVec.Px()*(1-smearFac*corrFac);
        metPy+= lVec.Py()*(1-smearFac*corrFac);

        lVec *= smearFac*corrFac;

        //Calculate HT for miniAOD
        HT+=lVec.Pt();

        if(lVec.Pt() > ptCut and abs(lVec.Eta()) < etaCut){
            //Fatjet four momentum components
            JetfloatVariables[0].push_back(lVec.E());   //Energy
            JetfloatVariables[1].push_back(lVec.Px());  //Px
            JetfloatVariables[2].push_back(lVec.Py());  //Py
            JetfloatVariables[3].push_back(lVec.Pz());  //Pz

            //Check for btag
            float DeepBValue = 0;

            if(isNANO){ 
                DeepBValue = jetDeepBValue->At(i);
            } 

            else{
                for(std::string disc: {"pfDeepFlavourJetTags:probb", "pfDeepFlavourJetTags:probbb","pfDeepFlavourJetTags:problepb"}){
                    DeepBValue +=  jets->at(i).bDiscriminator(disc);
                }
            }

            JetboolVariables[0].push_back(bTagCuts[AK4][era][0] < DeepBValue);
            JetboolVariables[1].push_back(bTagCuts[AK4][era][1] < DeepBValue);
            JetboolVariables[2].push_back(bTagCuts[AK4][era][2] < DeepBValue);

            if(!isData){
                //btag SF
                JetfloatVariables[6].push_back(looseReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                JetfloatVariables[7].push_back(looseReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                JetfloatVariables[8].push_back(looseReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));

                if(!isSyst){
                    JetfloatVariables[9].push_back(looseReader[AK4].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    JetfloatVariables[10].push_back(looseReader[AK4].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    JetfloatVariables[11].push_back(mediumReader[AK4].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    JetfloatVariables[12].push_back(mediumReader[AK4].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    JetfloatVariables[13].push_back(tightReader[AK4].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                    JetfloatVariables[14].push_back(tightReader[AK4].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(lVec.Eta()), lVec.Pt()));
                }


                if(isNANO) JetfloatVariables[5].push_back(SetGenParticles(lVec, i, 6, AK4));
                else{
                    event->getByToken(genParticleToken, genParts);
                    JetfloatVariables[5].push_back(SetGenParticles(lVec, i, 6, AK4, *genParts));
                }
            }

            //Check overlap with AK4 valid jets
            for(unsigned int j = 0; j < FatJetfloatVariables[0].size(); j++){
                TLorentzVector FatlVec;
                FatlVec.SetPxPyPzE(FatJetfloatVariables[1][j], FatJetfloatVariables[2][j], FatJetfloatVariables[3][j], FatJetfloatVariables[0][j]);

                if(FatlVec.DeltaR(lVec) < 1.2){
                    JetfloatVariables[4].push_back(j);
                    nSubJets++;
                }

                else{
                    JetfloatVariables[4].push_back(-1.);
                }
            }
        } 
    }

    //Set HT
    if(isNANO){
        HT = *valueHT->Get();
    }

    for(CutFlow& cutflow: cutflows){
        //Check if one combination of jet and fatjet number is fullfilled
        if(JetfloatVariables[0].size() - nSubJets >= cutflow.nMinJet and FatJetfloatVariables[0].size() == cutflow.nMinFatjet){
            if(cutflow.passed){
                std::string cutName("N^{AK4}_{jet} >= " + std::to_string(cutflow.nMinJet) + " && N^{AK8}_{jet} == " + std::to_string(cutflow.nMinFatjet));

                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
            }
        }

        else{
            cutflow.passed = false;
        }
    }
}

void JetAnalyzer::EndJob(TFile* file){
    for(const JetType& type: {AK4, AK8}){
        delete jetCorrector[type];
    }
}
