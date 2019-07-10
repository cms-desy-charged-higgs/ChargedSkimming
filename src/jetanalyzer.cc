#include <ChargedHiggs/Skimming/interface/jetanalyzer.h>

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, const std::vector<std::pair<unsigned int, unsigned int>> minNJet, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNJet(minNJet)
    {}

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, const std::vector<std::pair<unsigned int, unsigned int>> minNJet, std::vector<jToken>& jetTokens, std::vector<genjToken>& genjetTokens, mToken &metToken, edm::EDGetTokenT<double> &rhoToken):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNJet(minNJet),
    jetTokens(jetTokens),
    metToken(metToken),
    rhoToken(rhoToken)
    {}


//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetAnalyzer::CorrectEnergy(TLorentzVector &jet, const float &rho, float &area, const JetType &type){
    std::vector<JetCorrectorParameters> corrVec;

    for(std::string fileName: JEC[type][era]){
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
float JetAnalyzer::SmearEnergy(TLorentzVector &jet, const float &rho, const float &coneSize, const JetType &type, const std::vector<reco::GenJet> &genJet){
    //Configure jet SF reader
    jetParameter.setJetPt(jet.Pt()).setJetEta(jet.Eta()).setRho(rho);

    float dR;
    TLorentzVector matchedGenJet;
    float reso = resolution[type].getResolution(jetParameter);
    float resoSF = resolution_sf[type].getScaleFactor(jetParameter);
    float smearFac = 1.; 

    float genPt, genPhi, genEta, genMass;
    unsigned int size;

    if(isNANO) size = (type == AK4) ? genJetPt->GetSize(): genFatJetPt->GetSize();
    else size = genJet.size();

    //Loop over all gen jets and find match
    for(unsigned int i = 0; i < size; i++){
        if(isNANO){
            genPt = (type == AK4) ? genJetPt->At(i): genFatJetPt->At(i);
            genPhi = (type == AK4) ? genJetPhi->At(i): genFatJetPhi->At(i);
            genEta = (type == AK4) ? genJetEta->At(i): genFatJetEta->At(i);
            genMass = (type == AK4) ? genJetMass->At(i): genFatJetMass->At(i);
        }

        else{
            genPt = genJet.at(i).pt();
            genPhi = genJet.at(i).phi();
            genEta = genJet.at(i).eta();
            genMass = genJet.at(i).mass();
        }

        dR = std::sqrt(std::pow(jet.Phi()-genPhi, 2) + std::pow(jet.Eta()-genEta, 2));

        //Check if jet and gen jet are matched
        if(dR < coneSize/2. and abs(jet.Pt() - genPt) < 3.*reso*jet.Pt()){
            matchedGenJet.SetPtEtaPhiM(genPt, genEta, genPhi, genMass);
            break;
        }
    }  

    //If you found gen matched 
    if(matchedGenJet != TLorentzVector()){
        smearFac = 1.+(resoSF-1)*(jet.Pt() - matchedGenJet.Pt())/jet.Pt(); 
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

void JetAnalyzer::SetGenParticles(Jet &validJet, const int &i){
    //Check if gen matched particle exist
    if(jetGenIdx->At(i) != -1){
        int genPartIdx = -1;
        float dR;
        
        //Find Gen particle to gen Jet
        for(unsigned int index=0; index < genPt->GetSize(); index++){
            std::bitset<32> isFirstCopy(genStatus->At(index));

            if(!isFirstCopy[12]){
                continue;
            }

            dR = std::sqrt(std::pow(genPhi->At(index) -genJetPhi->At(jetGenIdx->At(i)), 2) + std::pow(genEta->At(index)-genJetEta->At(jetGenIdx->At(i)), 2));

            if(dR < 0.4 and abs(genID->At(index)) == 5){
                bool alreadyMatched = std::find(alreadyMatchedJet.begin(), alreadyMatchedJet.end(), genPartIdx) != alreadyMatchedJet.end();

                if(!alreadyMatched){
                    genPartIdx = index;
                    alreadyMatchedJet.push_back(index);
                    break;
                }
            }
        }

        if(genPartIdx != -1){
            int idxMotherB = genMotherIdx->At(genPartIdx);

            while(abs(genID->At(idxMotherB)) == 5){
                idxMotherB = genMotherIdx->At(idxMotherB);
            }

            //Check if gen jet from small higgs
            if(abs(genID->At(idxMotherB)) == 25){
                int idxMotherh = genMotherIdx->At(idxMotherB);

                while(abs(genID->At(idxMotherh)) == 25){
                    idxMotherh = genMotherIdx->At(idxMotherh);
                }

                if(abs(genID->At(idxMotherh)) == 37){
                    validJet.isFromh1 = true;
                    validJet.isFromh2 = false;
                }

                else{
                    validJet.isFromh1 = false;
                    validJet.isFromh2 = true;
                }
            }
        }
    }
}

void JetAnalyzer::BeginJob(TTree* tree, bool &isData){
    JEC = {
            {AK4, {
                {2017, {filePath + "/JEC/YYYY_L1FastJet_AK5PF.txt",   
                        filePath + "/JEC/YYYY_L1FastJet_AK5PF.txt"}
                }
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
        fatJetCSV = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_btagDeepB");
        fatJetTau1 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau1");
        fatJetTau2 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau2");
        fatJetTau3 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau3");

        jetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_pt");
        jetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_eta");
        jetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_phi");
        jetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_mass");
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
    alreadyMatchedJet.clear();

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Jet>> jets;
    edm::Handle<std::vector<pat::Jet>> fatJets;
    edm::Handle<std::vector<reco::GenJet>> genJets;
    edm::Handle<std::vector<reco::GenJet>> genfatJets;
    edm::Handle<std::vector<pat::MET>> MET;
    edm::Handle<double> rho;

    if(!isNANO){
        event->getByToken(jetTokens[0], jets);
        event->getByToken(jetTokens[1], fatJets);
        event->getByToken(metToken, MET);
        event->getByToken(rhoToken, rho);

        if(!isData){
            event->getByToken(genjetTokens[0], genJets);
            event->getByToken(genjetTokens[1], genfatJets);
        }
    }

    //JER smearing
    float smearFac = 1.;

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
    
        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(fatJet.fourVec, *jetRho->Get(), 0.8, AK8) : SmearEnergy(fatJet.fourVec, *rho, 0.8, AK8, *genfatJets);
            fatJet.fourVec *= smearFac;
        }

        if(fatJet.fourVec.M() > 50. and abs(fatJet.fourVec.Eta()) < etaCut){
            //Check for btag
            fatJet.isLooseB = isNANO ? bTagCuts[AK8][era][0] < fatJetCSV->At(i) : bTagCuts[AK8][era][0] < fatJets->at(0).bDiscriminator("pfDeepCSVJetTags:probb + pfDeepCSVJetTags:probbb");
            fatJet.isMediumB = isNANO ? bTagCuts[AK8][era][1] < fatJetCSV->At(i) : bTagCuts[AK8][era][0] < fatJets->at(0).bDiscriminator("pfDeepCSVJetTags:probb + pfDeepCSVJetTags:probbb");

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

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(jetCand.fourVec,  *jetRho->Get(), 0.4, AK4) : SmearEnergy(jetCand.fourVec, *rho, 0.4, AK4, *genJets);
            jetCand.fourVec *= smearFac;
        }

        if(!isData){
            //Correct met
            metPx+= std::cos(phi)*pt - std::cos(jetCand.fourVec.Phi())*jetCand.fourVec.Pt();
            metPy+= std::sin(phi)*pt - std::sin(jetCand.fourVec.Phi())*jetCand.fourVec.Pt();
        }

        if(jetCand.fourVec.Pt() > ptCut and abs(jetCand.fourVec.Eta()) < etaCut){
            //bool for cleaning fatjets from AK4 jets
            bool isCleaned = true;

            float DeepBValue = isNANO ? jetDeepBValue->At(i) : jets->at(i).bDiscriminator("      pfDeepFlavourJetTags:probb + pfDeepFlavourJetTags:probbb + pfDeepFlavourJetTags:problepb");

            std::cout << DeepBValue << std::endl;

            //Check for btag
            jetCand.isLooseB = bTagCuts[AK4][era][0] < DeepBValue;
            jetCand.isMediumB = bTagCuts[AK4][era][1] < DeepBValue;
            jetCand.isTightB = bTagCuts[AK4][era][2] < DeepBValue;

            if(!isData){
                //btag SF
                jetCand.loosebTagSF = looseReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(jetCand.fourVec.Eta()), jetCand.fourVec.Pt());
                jetCand.mediumbTagSF = mediumReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(jetCand.fourVec.Eta()), jetCand.fourVec.Pt());
                jetCand.tightbTagSF = tightReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(jetCand.fourVec.Eta()), jetCand.fourVec.Pt());


                //Save gen particle information
                //SetGenParticles(jetCand, i);
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
    //HT = *valueHT->Get();

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
