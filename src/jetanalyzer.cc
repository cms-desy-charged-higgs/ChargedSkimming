#include <ChargedHiggs/Skimming/interface/jetanalyzer.h>

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, const std::vector<std::pair<unsigned int, unsigned int>> minNJet, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNJet(minNJet)
    {}

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, const std::vector<std::pair<unsigned int, unsigned int>> minNJet, std::vector<jToken>& jetTokens, std::vector<genjToken>& genjetTokens, mToken &metToken, edm::EDGetTokenT<double> &rhoToken, genPartToken& genParticleToken):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNJet(minNJet),
    jetTokens(jetTokens),
    metToken(metToken),
    rhoToken(rhoToken),
    genParticleToken(genParticleToken)
    {}


//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetAnalyzer::CorrectEnergy(const TLorentzVector &jet, const float &rho, const float &area, const JetType &type, const int& runNumber){
    std::vector<JetCorrectorParameters> corrVec;

    for(std::string fileName: isData? JECDATA[type][era] : JECMC[type][era]){
        if(fileName.find("@") != std::string::npos){
            for(std::pair<std::string, std::pair<int, int>> eraNames: runEras[era]){
                if(eraNames.second.first <= runNumber and runNumber <= eraNames.second.second){
                    fileName.replace(fileName.find("@"), 1, eraNames.first);
                }
            }
        }

        corrVec.push_back(JetCorrectorParameters(fileName));
    }

    FactorizedJetCorrector jetCorrector(corrVec);
    jetCorrector.setJetPt(jet.Pt());
    jetCorrector.setJetEta(jet.Eta());
    jetCorrector.setRho(rho);
    jetCorrector.setJetA(area);

    double correction = jetCorrector.getCorrection();
    
    return correction;
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
float JetAnalyzer::SmearEnergy(const TLorentzVector &jet, const float &rho, const float &coneSize, const JetType &type, const std::vector<reco::GenJet> &genJets){
    //Configure jet SF reader
    jetParameter.setJetPt(jet.Pt()).setJetEta(jet.Eta()).setRho(rho);

    float dR;
    float reso = resolution[type].getResolution(jetParameter);
    float resoSF = resolution_sf[type].getScaleFactor(jetParameter);
    float smearFac = 1.; 

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

void JetAnalyzer::SetGenParticles(Jet& validJet, const int &i, const int &pdgID, const JetType &type, const std::vector<reco::GenParticle>& genParticle){
    //Check if gen matched particle exist
    if(genJet[type].Pt() != 0){
        validJet.genJet = genJet[type];
        float dR;
        
        //Find Gen particle to gen Jet
        int size=isNANO ? genPt->GetSize() : genParticle.size();

        for(int i=0; i < size; i++){
            const reco::Candidate* parton=NULL;
            int index=0;

            int ID = isNANO ? abs(genID->At(genMotherIdx->At(i))) : abs(genParticle.at(i).pdgId());
    
            if(ID == pdgID){
                if(isNANO) index = LastCopy(i, pdgID);
                else parton = LastCopy(genParticle.at(i), pdgID);
            }
    
            else continue;

            float pt, eta, phi, m;
            pt = isNANO ? genPt->At(index) : parton->pt();
            phi = isNANO ? genPhi->At(index) : parton->phi();
            eta = isNANO ? genEta->At(index) : parton->eta();
            m = isNANO ? genMass->At(index) : parton->mass();

            dR = std::sqrt(std::pow(phi - genJet[type].Phi(), 2) + std::pow(eta-genJet[type].Eta(), 2));
            float rMin = type == AK4 ? 0.3 : 0.4;

            if(dR <  rMin){
                TLorentzVector matchedParton;
                matchedParton.SetPtEtaPhiM(pt, eta, phi, m);
                validJet.genPartons.push_back(matchedParton);

                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(parton->mother()->pdgId());

                if(motherID == 25){
                    const reco::Candidate* hBoson=NULL;
                    int index=0;
                                
                    if(isNANO) index = LastCopy(index, 25);
                    else hBoson = LastCopy(parton->mother(), 25);

                    int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(hBoson->mother()->pdgId());  

                    if(motherID == 37){
                        validJet.isFromh1.push_back(true);
                        validJet.isFromh2.push_back(false);
                    }

                    else{
                        validJet.isFromh1.push_back(false);
                        validJet.isFromh2.push_back(true);
                    }
                }
            }

            if(validJet.genPartons.size()==1 and type==AK4){
                break;
            }

            if(validJet.genPartons.size()==2 and type==AK8){
                break;
            }
        }
    }
}


void JetAnalyzer::BeginJob(TTree* tree, bool &isData){
    JECMC = {
            {AK4, {
                {2017, {filePath + "/JEC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt"}}
                },
            },  

            {AK8, {
                {2017, {filePath + "/JEC/Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt"}}
                },
            },               
    };

    JECDATA = {
            {AK4, {
                {2017, {filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L1FastJet_AK4PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2Relative_AK4PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L3Absolute_AK4PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2L3Residual_AK4PFchs.txt"}}
                },
            }, 

            {AK8, {
                {2017, {filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L1FastJet_AK8PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2Relative_AK8PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L3Absolute_AK8PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2L3Residual_AK8PFchs.txt"}}
                },
            },               
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

    JMESF = {
            {AK4, {
                {2017, filePath + "/JME/Fall17_V3_MC_SF_AK4PFchs.txt"},
                }
            },       

            {AK8, {
                {2017, filePath + "/JME/Fall17_V3_MC_SF_AK8PFchs.txt"},
                }
            },             
    };

    JMEPtReso = {
            {AK4, {
                {2017, filePath + "/JME/Fall17_V3_MC_PtResolution_AK4PFchs.txt"},
                }
            },       

            {AK8, {
                {2017, filePath + "/JME/Fall17_V3_MC_PtResolution_AK8PFchs.txt"},
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
        resolution[type] = JME::JetResolution(JMEPtReso[type][era]);
        resolution_sf[type] = JME::JetResolutionScaleFactor(JMESF[type][era]);
    }

    //Set Branches of output tree
    tree->Branch("jet", &validJets);
    tree->Branch("subjet", &subJets);
    tree->Branch("fatjet", &validFatJets);
    tree->Branch("met", &met);
    tree->Branch("HT", &HT);
}


bool JetAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event){
    //Clear jet vector
    validJets.clear();
    subJets.clear();
    validFatJets.clear();
    HT=0;
    int runNumber = isNANO ? *run->Get() : event->eventAuxiliary().id().run(); 

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Jet>> jets;
    edm::Handle<std::vector<pat::Jet>> fatJets;
    edm::Handle<std::vector<reco::GenJet>> genJets;
    edm::Handle<std::vector<reco::GenJet>> genfatJets;
    edm::Handle<std::vector<pat::MET>> MET;
    edm::Handle<double> rho;
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isNANO){
        event->getByToken(jetTokens[0], jets);
        event->getByToken(jetTokens[1], fatJets);
        event->getByToken(metToken, MET);
        event->getByToken(rhoToken, rho);

        if(!isData){
            event->getByLabel(edm::InputTag("slimmedGenJets"), genJets);
            event->getByLabel(edm::InputTag("slimmedGenJetsAK8"), genfatJets);
        }
    }

    //JER smearing
    float smearFac = 1.;
    float corrFac = 1.;

    //MET values not correct for JER yet
    float metPx = isNANO ? *metPt->Get()*std::cos(*metPhi->Get()) : MET->at(0).uncorPx();
    float metPy = isNANO ? *metPt->Get()*std::sin(*metPhi->Get()) : MET->at(0).uncorPy();
    
    float fatJetSize = isNANO ? fatJetPt->GetSize() : fatJets->size();

    //Loop over all fat jets
    for(unsigned int i = 0; i < fatJetSize; i++){
        float fatPt = isNANO ? fatJetPt->At(i) : fatJets->at(i).pt();
        float fatEta = isNANO ? fatJetEta->At(i) : fatJets->at(i).eta();
        float fatPhi = isNANO ? fatJetPhi->At(i) : fatJets->at(i).phi();
        float fatMass = isNANO ? fatJetMass->At(i) : fatJets->at(i).mass();

        //Define fat jet
        FatJet fatJet;
        fatJet.fourVec.SetPtEtaPhiM(fatPt, fatEta, fatPhi, fatMass);
    
        corrFac = isNANO ? CorrectEnergy(fatJet.fourVec,  *jetRho->Get(), fatJetArea->At(i), AK8, runNumber) : CorrectEnergy(fatJet.fourVec, *rho, jets->at(i).jetArea(), AK8, runNumber);

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(fatJet.fourVec*corrFac, *jetRho->Get(), 0.8, AK8) : SmearEnergy(fatJet.fourVec*corrFac, *rho, 0.8, AK8, *genfatJets);
            fatJet.fourVec *= smearFac*corrFac;

            if(isNANO) SetGenParticles(fatJet, i, 5, AK8);
            else{
                event->getByToken(genParticleToken, genParts);
                SetGenParticles(fatJet, i, 5, AK8, *genParts);
            }
        }

        else fatJet.fourVec *= corrFac;

        if(fatJet.fourVec.M() > 50. and abs(fatJet.fourVec.Eta()) < etaCut){
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
            fatJet.isLooseB = bTagCuts[AK8][era][0] < DeepCSV;
            fatJet.isMediumB =  bTagCuts[AK8][era][1] < DeepCSV;

            //Nsubjettiness
            fatJet.oneSubJettiness = isNANO ? fatJetTau1->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1");
            fatJet.twoSubJettiness = isNANO ? fatJetTau2->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2");
            fatJet.threeSubJettiness = isNANO ? fatJetTau3->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3");
            if(!isData){
                //btag SF
                fatJet.loosebTagSF = looseReader[AK8].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(fatJet.fourVec.Eta()), fatJet.fourVec.Pt());
                fatJet.mediumbTagSF = mediumReader[AK8].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(fatJet.fourVec.Eta()), fatJet.fourVec.Pt());
            }

            validFatJets.push_back(fatJet);
        }
    }

    float jetSize = isNANO ? jetPt->GetSize() : jets->size();

    //Loop over all jets
    for(unsigned int i = 0; i < jetSize; i++){
        float pt = isNANO ? jetPt->At(i) : jets->at(i).pt();
        float eta = isNANO ? jetEta->At(i) : jets->at(i).eta();
        float phi = isNANO ? jetPhi->At(i) : jets->at(i).phi();
        float mass = isNANO ? jetMass->At(i) : jets->at(i).mass();

        //Define here already jet, because of smearing of 4-vec
        Jet jetCand;
        jetCand.fourVec.SetPtEtaPhiM(pt, eta, phi, mass);

        corrFac = isNANO ? CorrectEnergy(jetCand.fourVec,  *jetRho->Get(), jetArea->At(i), AK4, runNumber) : CorrectEnergy(jetCand.fourVec, *rho, jets->at(i).jetArea(), AK4, runNumber);

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(jetCand.fourVec*corrFac,  *jetRho->Get(), jetArea->At(i), AK4) : SmearEnergy(jetCand.fourVec, *rho, 0.4, AK4, *genJets);

            jetCand.fourVec*=smearFac*corrFac;
        }

        else jetCand.fourVec*=corrFac;

        //Correct met
        metPx+= jetCand.fourVec.Px()*(1-smearFac*corrFac);
        metPy+= jetCand.fourVec.Py()*(1-smearFac*corrFac);

        jetCand.fourVec *= smearFac*corrFac;

        //Calculate HT for miniAOD
        HT+=jetCand.fourVec.Pt();

        if(jetCand.fourVec.Pt() > ptCut and abs(jetCand.fourVec.Eta()) < etaCut){
            //bool for cleaning fatjets from AK4 jets
            bool isCleaned = true;

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

            jetCand.isLooseB = bTagCuts[AK4][era][0] < DeepBValue;
            jetCand.isMediumB = bTagCuts[AK4][era][1] < DeepBValue;
            jetCand.isTightB = bTagCuts[AK4][era][2] < DeepBValue;

            if(!isData){
                //btag SF
                jetCand.loosebTagSF = looseReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(jetCand.fourVec.Eta()), jetCand.fourVec.Pt());
                jetCand.mediumbTagSF = mediumReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(jetCand.fourVec.Eta()), jetCand.fourVec.Pt());
                jetCand.tightbTagSF = tightReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(jetCand.fourVec.Eta()), jetCand.fourVec.Pt());


                if(isNANO) SetGenParticles(jetCand, i, 5, AK4);
                else{
                    event->getByToken(genParticleToken, genParts);
                    SetGenParticles(jetCand, i, 5, AK4, *genParts);
                }
            }

            //Check overlap with AK4 valid jets
            for(unsigned int j = 0; j < validFatJets.size(); j++){
                if(validFatJets[j].fourVec.DeltaR(jetCand.fourVec) < 1.2){
                    isCleaned = false;
                    jetCand.fatJetIdx = j;
                    subJets.push_back(jetCand);
                }
            }

            //If not cleaned dont use this jet
            if(!isCleaned){
                continue;
            }

            validJets.push_back(jetCand);
        } 
    }


    //Set HT
    if(isNANO){
        HT = *valueHT->Get();
    }

    //Set met
    met.SetPxPyPzE(metPx, metPy, 0 , 0);

    //Check if one combination of jet and fatjet number is fullfilled

    for(std::pair<unsigned int, unsigned int> minN: minNJet){
        if(validJets.size() >= minN.first && validFatJets.size() == minN.second){
            std::string cutName("N^{AK4}_{jet} >= " + std::to_string(minN.first) + " && N^{AK8}_{jet} == " + std::to_string(minN.second));

            cutflow.first->Fill(cutName.c_str(), cutflow.second);
            return true;
        }
    }
    
    return false;
}


void JetAnalyzer::EndJob(TFile* file){
}
