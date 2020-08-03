#include <ChargedSkimming/Skimming/interface/jetanalyzer.h>

JetAnalyzer::JetAnalyzer(const int& era, const float& ptCut, const float& etaCut, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

JetAnalyzer::JetAnalyzer(const int& era, const float& ptCut, const float& etaCut, std::vector<jToken>& jetTokens, std::vector<genjToken>& genjetTokens, mToken &metToken, edm::EDGetTokenT<double> &rhoToken, genPartToken& genParticleToken, secvtxToken& vertexToken, const std::string& systematic):
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

            if(systematic.find("Up") != std::string::npos) isUp = true;
            else isUp = false;
        }
    }


void JetAnalyzer::SetCorrector(const JetType& type, const int& runNumber){
    std::vector<JetCorrectorParameters> corrVec;

    for(std::string fileName: isData? JECDATA : JECMC){
        if(fileName.find("@") != std::string::npos){
            for(std::pair<std::string, std::pair<int, int>> eraNames: runEras[era]){
                if(eraNames.second.first <= runNumber and runNumber <= eraNames.second.second){
                    fileName.replace(fileName.find("@"), 1, eraNames.first);
                }
            }
        }

        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        fileName = filePath + fileName;

        corrVec.push_back(JetCorrectorParameters(fileName));
    }

    jetCorrector[type] = new FactorizedJetCorrector(corrVec);
}

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetAnalyzer::CorrectEnergy(const float& pt, const float& eta, const float& rho, const float& area, const JetType &type){
    jetCorrector[type]->setJetPt(pt);
    jetCorrector[type]->setJetEta(eta);
    jetCorrector[type]->setRho(rho);
    jetCorrector[type]->setJetA(area);
    
    return jetCorrector[type]->getCorrection();
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
float JetAnalyzer::SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, const float& coneSize, const JetType& type, const std::vector<reco::GenJet>& genJets){
    jetParameter.setJetPt(pt).setJetEta(eta).setRho(rho);

    float reso = resolution[type].getResolution(jetParameter);
    float resoSF; 
    if(!isJERsyst) resoSF = resolution_sf[type].getScaleFactor(jetParameter);
    else resoSF = resolution_sf[type].getScaleFactor(jetParameter, isUp ? Variation::UP : Variation::DOWN);
        
    float smearFac = 1.; 
    float dR;
    float genPt, genPhi, genEta, genMass;
    unsigned int size;
    bool isMatched = false;

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

        dR = std::sqrt(std::pow(phi-genPhi, 2) + std::pow(eta-genEta, 2));

        //Check if jet and gen jet are matched
        if(dR < coneSize/2. and abs(pt - genPt) < 3.*reso*pt){
            genJet[type] = ROOT::Math::PtEtaPhiMVector(genPt, genEta, genPhi, genMass);
            isMatched = true;
            break;
        }
    }  

    //If you found gen matched 
    if(isMatched){
        smearFac = 1.+(resoSF-1)*(pt - genPt)/pt; 
    }

    //If no match, smear with gaussian pdf
    else if(resoSF > 1.){
        std::default_random_engine generator;
        std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));

        smearFac = 1. + gaus(generator);
    }


    //Check if direction of jet not changed
    if(pt*smearFac < 1e-2){
        smearFac = 1e-2/pt;
    }

    return smearFac;
}

void JetAnalyzer::SetGenParticles(const int& i, const float& pt, const float& eta, const float& phi, const std::vector<int>& pdgID, const JetType &type, const std::vector<reco::GenParticle>& genParticle){
    float dR;

    partID.at(type).push_back(-99.);
    mothID.at(type).push_back(-99.);
    grandID.at(type).push_back(-99.);

    float bestDR = 30.; float bestPt = 100.;
    const reco::Candidate* parton=nullptr;
    int index=-1;

    //Check if gen matched particle exist
    if(genJet[type].Pt() != 0){
        //Find Gen particle to gen Jet
        const char& size = isNANO ? genPt->GetSize() : genParticle.size();

        for(int i=0; i < size; i++){
            const int ID = isNANO ? abs(genID->At(i)) : abs(genParticle.at(i).pdgId());

            if(std::find(pdgID.begin(), pdgID.end(), ID) != pdgID.end()){
                if(isNANO) index = FirstCopy(i, ID);
                else parton = FirstCopy(genParticle.at(i), ID);
            }    

            else continue;

            const float& ptGen = isNANO ? genPt->At(index) : parton->pt();
            const float& phiGen = isNANO ? genPhi->At(index) : parton->phi();
            const float& etaGen = isNANO ? genEta->At(index) : parton->eta();

            dR = BaseAnalyzer::DeltaR(etaGen, phiGen, genJet[type].Eta(), genJet[type].Phi());
            const float& dPt = abs((genJet[type].Pt() - ptGen)/genJet[type].Pt());

            if(dR <  bestDR and dPt <  bestPt){
                if(isNANO){
                    if(std::find(alreadySeenNANO.begin(), alreadySeenNANO.end(), index) != alreadySeenNANO.end()) continue; 
                    else alreadySeenNANO.push_back(index);
                }

                else{
                    if(std::find(alreadySeenMINI.begin(), alreadySeenMINI.end(), parton) != alreadySeenMINI.end()) continue; 
                    else alreadySeenMINI.push_back(parton);
                }

                bestDR = dR;    
                bestPt = dPt;

                const reco::Candidate* motherPart=nullptr;
                int motherIdx=0;

                const int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(parton->mother()->pdgId());

                if(isNANO) motherIdx = FirstCopy(genMotherIdx->At(index), motherID);
                else motherPart = FirstCopy(parton->mother(), motherID);

                const int grandMotherID = isNANO ? abs(genID->At(genMotherIdx->At(motherIdx))) : motherPart->mother() != nullptr ? abs(motherPart->mother()->pdgId()) : -99.; 

                partID.at(type).back() = ID;
                mothID.at(type).back() = motherID;
                grandID.at(type).back() = grandMotherID;
                break;
            }
        }
    }

    if(isNANO and index != -1) alreadySeenNANO.push_back(index);
    if(!isNANO and parton != nullptr) alreadySeenMINI.push_back(parton);
}

void JetAnalyzer::BeginJob(std::vector<TTree*>& trees, bool& isData, const bool& isSyst){
    //Read in json config with sf files
    boost::property_tree::ptree sf; 
    boost::property_tree::read_json(std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/config/sf.json", sf);

    //Get for JEC files names
    JECMC = Util::GetVector<std::string>(sf, "Jet.JEC.MC." + std::to_string(era));
    JECDATA = Util::GetVector<std::string>(sf, "Jet.JEC.DATA." + std::to_string(era));

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

            jetFlavour = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_partonFlavour");    
            jetGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_genJetIdx");
        }

        //Set TTreeReader for genpart and trigger obj from baseanalyzer
        SetCollection(this->isData);
    }

    //Set configuration for bTagSF reader  
    CSVReader=BTagCSVReader(filePath + sf.get<std::string>("Jet.BTag.CSV." + std::to_string(era)));
    DeepReader=BTagCSVReader(filePath +sf.get<std::string>("Jet.BTag.DeepJet." + std::to_string(era)));

    //Cut Values for b tagging
    CSVLoose = {{2016, 0.2217}, {2017, 0.1522}, {2018, 0.1241}};
    CSVMedium = {{2016, 0.6321}, {2017, 0.4941}, {2018, 0.4184}};
    CSVTight = {{2016, 0.8953}, {2017, 0.8001}, {2018, 0.7527}};
    DeepLoose = {{2016, 0.0614}, {2017, 0.0521}, {2018, 0.0494}};
    DeepMedium = {{2016, 0.3093}, {2017, 0.3033}, {2018, 0.2770}};
    DeepTight = {{2016, 0.7221}, {2017, 0.7489}, {2018, 0.7264}};

    if(!isData){
        bTotal = new TH2F("nTrueB", "TotalB", 20, ptCut, 400, 20, -etaCut, etaCut);

        for(const std::string& taggerName :  {"Deep", "CSV"}){
            bTagEffLoose.push_back(new TH2F(("nLoose" + taggerName + "bTag").c_str(), "", 20, ptCut, 400, 20, -etaCut, etaCut));
            bTagEffMedium.push_back(new TH2F(("nMedium" + taggerName + "bTag").c_str(), "", 20, ptCut, 400, 20, -etaCut, etaCut));
            bTagEffTight.push_back(new TH2F(("nTight" + taggerName + "bTag").c_str(), "", 20, ptCut, 400, 20, -etaCut, etaCut));
        }
    }
    
    for(JetType type: {AK4, AK8}){    
        //Set configuration for JER tools
        std::string fileName = filePath + sf.get<std::string>("Jet.JMEPtReso." + std::to_string(era));
        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        resolution[type] = JME::JetResolution(fileName);

        fileName = filePath + sf.get<std::string>("Jet.JME." + std::to_string(era));
        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        resolution_sf[type] = JME::JetResolutionScaleFactor(fileName);

        //Set object to get JEC uncertainty
        if(jecSyst != ""){
            std::string fileName = filePath + sf.get<std::string>("Jet.JECUNC." + std::to_string(era));
            fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");  

            jecUnc[type] = new JetCorrectionUncertainty(JetCorrectorParameters(fileName, "Total"));
        }

        else jecUnc[type] = nullptr;
    }

    //Set output names
    floatVar = {
        {{"Jet", "Pt"}, Pt[AK4]},
        {{"Jet", "Eta"}, Eta[AK4]},
        {{"Jet", "Phi"}, Phi[AK4]},
        {{"Jet", "Mass"}, Mass[AK4]},
        {{"Jet", "corrJEC"}, corrJEC[AK4]},
        {{"Jet", "corrJER"}, corrJER[AK4]},
        {{"Jet", "looseCSVbTagSF"}, looseCSVbTagSF[AK4]},
        {{"Jet", "mediumCSVbTagSF"}, mediumCSVbTagSF[AK4]},
        {{"Jet", "tightCSVbTagSF"}, tightCSVbTagSF[AK4]},
        {{"Jet", "looseDeepbTagSF"}, looseDeepbTagSF[AK4]},
        {{"Jet", "mediumDeepbTagSF"}, mediumDeepbTagSF[AK4]},
        {{"Jet", "tightDeepbTagSF"}, tightDeepbTagSF[AK4]},
        {{"Jet", "DeepScore"}, DeepScore[AK4]},
        {{"Jet", "CSVScore"}, CSVScore[AK4]},
        {{"SubJet", "Pt"}, Pt[SUBAK4]},
        {{"SubJet", "Eta"}, Eta[SUBAK4]},
        {{"SubJet", "Phi"}, Phi[SUBAK4]},
        {{"SubJet", "Mass"}, Mass[SUBAK4]},
        {{"SubJet", "corrJEC"}, corrJEC[SUBAK4]},
        {{"SubJet", "corrJER"}, corrJER[SUBAK4]},
        {{"SubJet", "looseCSVbTagSF"}, looseCSVbTagSF[SUBAK4]},
        {{"SubJet", "mediumCSVbTagSF"}, mediumCSVbTagSF[SUBAK4]},
        {{"SubJet", "tightCSVbTagSF"}, tightCSVbTagSF[SUBAK4]},
        {{"SubJet", "looseDeepbTagSF"}, looseDeepbTagSF[SUBAK4]},
        {{"SubJet", "mediumDeepbTagSF"}, mediumDeepbTagSF[SUBAK4]},
        {{"SubJet", "tightDeepbTagSF"}, tightDeepbTagSF[SUBAK4]},
        {{"SubJet", "DeepScore"}, DeepScore[SUBAK4]},
        {{"SubJet", "CSVScore"}, CSVScore[SUBAK4]},
        {{"FatJet", "Pt"}, Pt[AK8]},
        {{"FatJet", "Eta"}, Eta[AK8]},
        {{"FatJet", "Phi"}, Phi[AK8]},
        {{"FatJet", "Mass"}, Mass[AK8]},
        {{"FatJet", "corrJEC"}, corrJEC[AK8]},
        {{"FatJet", "corrJER"}, corrJER[AK8]},
        {{"FatJet", "Njettiness1"}, Njettiness1[AK8]},
        {{"FatJet", "Njettiness2"}, Njettiness2[AK8]},
        {{"FatJet", "Njettiness3"}, Njettiness3[AK8]},
        {{"FatJet", "topVsHiggs"}, topVsHiggs[AK8]},
        {{"FatJet", "QCDVsHiggs"}, QCDVsHiggs[AK8]},
        {{"FatJet", "WVsHiggs"}, WVsHiggs[AK8]},
        {{"JetParticle", "Pt"}, Pt[PF]},
        {{"JetParticle", "Eta"}, Eta[PF]},
        {{"JetParticle", "Phi"}, Phi[PF]},
        {{"JetParticle", "Mass"}, Mass[PF]},
        {{"JetParticle", "Vx"}, Vx[PF]},
        {{"JetParticle", "Vy"}, Vy[PF]},
        {{"JetParticle", "Vz"}, Vz[PF]},
        {{"SecondaryVertex", "Pt"}, Pt[VTX]},
        {{"SecondaryVertex", "Eta"}, Eta[VTX]},
        {{"SecondaryVertex", "Phi"}, Phi[VTX]},
        {{"SecondaryVertex", "Mass"}, Mass[VTX]},
        {{"SecondaryVertex", "Vx"}, Vx[VTX]},
        {{"SecondaryVertex", "Vy"}, Vy[VTX]},
        {{"SecondaryVertex", "Vz"}, Vz[VTX]},
    };

    intVar = {
        {{"Jet", "TrueFlavour"}, TrueFlavour[AK4]},
        {{"Jet", "ParticleID"}, partID[AK4]},
        {{"Jet", "MotherID"}, mothID[AK4]},
        {{"Jet", "GrandMotherID"}, grandID[AK4]},
        {{"SubJet", "TrueFlavour"}, TrueFlavour[SUBAK4]},
        {{"SubJet", "FatJetIdx"}, FatJetIdx[SUBAK4]},
        {{"SubJet", "ParticleID"}, partID[SUBAK4]},
        {{"SubJet", "MotherID"}, mothID[SUBAK4]},
        {{"SubJet", "GrandMotherID"}, grandID[SUBAK4]},
        {{"JetParticle", "Charge"}, Charge[PF]},
        {{"JetParticle", "FatJetIdx"}, FatJetIdx[PF]},
        {{"SecondaryVertex", "FatJetIdx"}, FatJetIdx[VTX]},
        {{"FatJet", "ParticleID"}, partID[AK8]},
        {{"FatJet", "MotherID"}, mothID[AK8]},
        {{"FatJet", "GrandMotherID"}, grandID[AK8]},
    };

    if(!isSyst){
        std::map<std::pair<std::string, std::string>, std::vector<float>&> SFvariations = {
            {{"Jet", "looseDeepbTagSFUp"}, looseDeepbTagSFUp[AK4]},
            {{"Jet", "looseDeepbTagSFDown"}, looseDeepbTagSFUp[AK4]},
            {{"Jet", "mediumDeepbTagSFUp"}, mediumDeepbTagSFUp[AK4]},
            {{"Jet", "mediumDeepbTagSFDown"}, mediumDeepbTagSFUp[AK4]},
            {{"Jet", "tightDeepbTagSFUp"}, tightDeepbTagSFUp[AK4]},
            {{"Jet", "tightDeepbTagSFDown"}, tightDeepbTagSFUp[AK4]},
            {{"Jet", "looseCSVbTagSFUp"}, looseCSVbTagSFUp[AK4]},
            {{"Jet", "looseCSVbTagSFDown"}, looseCSVbTagSFDown[AK4]},
            {{"Jet", "mediumCSVbTagSFUp"}, mediumCSVbTagSFUp[AK4]},
            {{"Jet", "mediumCSVbTagSFDown"}, mediumCSVbTagSFDown[AK4]},
            {{"Jet", "tightCSVbTagSFUp"}, tightCSVbTagSFUp[AK4]},
            {{"Jet", "tightCSVbTagSFDown"}, tightCSVbTagSFDown[AK4]},
            {{"SubJet", "looseDeepbTagSFUp"}, looseDeepbTagSFUp[SUBAK4]},
            {{"SubJet", "looseDeepbTagSFDown"}, looseDeepbTagSFUp[SUBAK4]},
            {{"SubJet", "mediumDeepbTagSFUp"}, mediumDeepbTagSFUp[SUBAK4]},
            {{"SubJet", "mediumDeepbTagSFDown"}, mediumDeepbTagSFUp[SUBAK4]},
            {{"SubJet", "tightDeepbTagSFUp"}, tightDeepbTagSFUp[SUBAK4]},
            {{"SubJet", "tightDeepbTagSFDown"}, tightDeepbTagSFUp[SUBAK4]},
            {{"SubJet", "looseCSVbTagSFUp"}, looseCSVbTagSFUp[SUBAK4]},
            {{"SubJet", "looseCSVbTagSFDown"}, looseCSVbTagSFDown[SUBAK4]},
            {{"SubJet", "mediumCSVbTagSFUp"}, mediumCSVbTagSFUp[SUBAK4]},
            {{"SubJet", "mediumCSVbTagSFDown"}, mediumCSVbTagSFDown[SUBAK4]},
            {{"SubJet", "tightCSVbTagSFUp"}, tightCSVbTagSFUp[SUBAK4]},
            {{"SubJet", "tightCSVbTagSFDown"}, tightCSVbTagSFDown[SUBAK4]},
        };

        floatVar.insert(SFvariations.begin(), SFvariations.end());   
    }

    //Set Branches of output tree
    for(TTree* tree: trees){
        tree->Branch("Jet_Size", &nJets);
        tree->Branch("SubJet_Size", &nSubJets);
        tree->Branch("FatJet_Size", &nFatJets);

        for(std::pair<const std::pair<std::string, std::string>, std::vector<float>&>& var: floatVar){
            tree->Branch((var.first.first + "_" + var.first.second).c_str(), &var.second);
        }

        for(std::pair<const std::pair<std::string, std::string>, std::vector<char>&>& var: intVar){
            tree->Branch((var.first.first + "_" + var.first.second).c_str(), &var.second);
        }

        tree->Branch("MET_Pt", &metPT);
        tree->Branch("MET_Phi", &metPHI);
    }
}

void JetAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::pair<const std::pair<std::string, std::string>, std::vector<float>&>& var: floatVar){
        var.second.clear();
    }

    for(std::pair<const std::pair<std::string, std::string>, std::vector<char>&>& var: intVar){
        var.second.clear();
    }

    alreadySeenNANO.clear();
    alreadySeenMINI.clear();

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

    //Set necessary parameter
    float CSVBValue = 0, DeepBValue = 0, smearFac = 1., corrFac = 1., hScore = 0, topScore = 0, WScore = 0;
    int FatIdx = -1.;
    float metPx = 0, metPy;

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

    //MET values not correct for JER yet
    metPx = isNANO ? *metPt->Get()*std::cos(*metPhi->Get()) : MET->at(0).uncorPx();
    metPy = isNANO ? *metPt->Get()*std::sin(*metPhi->Get()) : MET->at(0).uncorPy();
    
    const char& fatJetSize = isNANO ? fatJetPt->GetSize() : fatJets->size();
    const char& jetSize = isNANO ? jetPt->GetSize() : jets->size();
        
    //Loop over all fat jets
    for(int i = 0; i < fatJetSize; i++){
        //JER smearing
        smearFac = 1., corrFac = 1.;

        const float fatPt = isNANO ? fatJetPt->At(i) : fatJets->at(i).pt();
        const float fatEta = isNANO ? fatJetEta->At(i) : fatJets->at(i).eta();
        const float fatPhi = isNANO ? fatJetPhi->At(i) : fatJets->at(i).phi();
        const float fatMass = isNANO ? fatJetMass->At(i) : fatJets->at(i).mass();
    
        corrFac = isNANO ? CorrectEnergy(fatPt, fatEta, *jetRho->Get(), fatJetArea->At(i), AK8) : CorrectEnergy(fatPt, fatEta, *rho, fatJets->at(i).jetArea(), AK8);

        //Get jet uncertainty
        if(jecUnc[AK8] !=  NULL){
            jecUnc[AK8]->setJetPt(corrFac*fatPt);
            jecUnc[AK8]->setJetEta(fatEta);
            const float& unc = jecUnc[AK8]->getUncertainty(isUp);
            corrFac *= isUp ? 1 + unc : 1 - unc;
        }

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(corrFac*fatPt, fatEta, fatPhi, *jetRho->Get(), 0.8, AK8) : SmearEnergy(corrFac*fatPt, fatEta, fatPhi, *rho, 0.8, AK8, *genfatJets);
        }

        if(smearFac*corrFac*fatPt > 170. and smearFac*corrFac*fatMass > 40. and abs(fatEta) < etaCut){
            if(!isNANO){
                hScore = fatJets->at(i).bDiscriminator("pfDeepBoostedJetTags:probHbb");
                topScore = 0.;
                WScore = 0.;

                for(const std::string& name : {"pfDeepBoostedJetTags:probTbcq", "pfDeepBoostedJetTags:probTbqq", "pfDeepBoostedJetTags:probTbc", "pfDeepBoostedJetTags:probTbq"}){
                    topScore += fatJets->at(i).bDiscriminator(name);
                }

                for(const std::string& name : {"pfDeepBoostedJetTags:probWcq", "pfDeepBoostedJetTags:probWqq"}){
                    WScore += fatJets->at(i).bDiscriminator(name);
                }

                topVsHiggs[AK8].push_back(hScore/(hScore + topScore));
                WVsHiggs[AK8].push_back(hScore/(hScore + WScore));
                QCDVsHiggs[AK8].push_back(fatJets->at(i).bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"));
            }
    
            //Fatjet four momentum components
            Pt[AK8].push_back(smearFac*corrFac*fatPt);
            Eta[AK8].push_back(fatEta);  
            Phi[AK8].push_back(fatPhi);  
            Mass[AK8].push_back(smearFac*corrFac*fatMass);
            corrJEC[AK8].push_back(corrFac);
            if(!isData) corrJER[AK8].push_back(smearFac);

            //Nsubjettiness
            Njettiness1[AK8].push_back(isNANO ? fatJetTau1->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1"));
            Njettiness2[AK8].push_back(isNANO ? fatJetTau2->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2"));
            Njettiness3[AK8].push_back(isNANO ? fatJetTau3->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3"));

            if(!isData){
                if(isNANO) SetGenParticles(i, Pt[AK8].back(), Eta[AK8].back(), Phi[AK8].back(), {24, 25, 23, 6}, AK8);
                else{
                    event->getByToken(genParticleToken, genParts);
                    SetGenParticles(i, Pt[AK8].back(), Eta[AK8].back(), Phi[AK8].back(), {24, 25, 23, 6}, AK8, *genParts);
                }
            }

            //Fill in particle flow candidates
            if(!isNANO){
                for(unsigned int j = 0; j < fatJets->at(i).numberOfDaughters(); j++){
                    const reco::Candidate* cand = fatJets->at(i).daughter(j);
                    
                    if(cand->numberOfDaughters() == 0){
                        //Jet particle four momentum components
                        Pt[PF].push_back(cand->pt());
                        Eta[PF].push_back(cand->eta());
                        Phi[PF].push_back(cand->phi());
                        Mass[PF].push_back(cand->mass());

                        //Jet particle vertex
                        Vx[PF].push_back(cand->vx());      
                        Vy[PF].push_back(cand->vy()); 
                        Vz[PF].push_back(cand->vz());
                        Charge[PF].push_back(cand->charge());
            
                        //Fat Jet Index
                        FatJetIdx[PF].push_back(Pt[AK8].size()-1);
                    }
                            
                    else{
                        for(unsigned int k = 0; k < cand->numberOfDaughters(); k++){
                            const reco::Candidate* cand2 = cand->daughter(k);

                            //Jet particle four momentum components
                            Pt[PF].push_back(cand2->pt());
                            Eta[PF].push_back(cand2->eta());
                            Phi[PF].push_back(cand2->phi());
                            Mass[PF].push_back(cand2->mass());

                            //Jet particle vertex
                            Vx[PF].push_back(cand2->vx());      
                            Vy[PF].push_back(cand2->vy()); 
                            Vz[PF].push_back(cand2->vz());
                            Charge[PF].push_back(cand2->charge());
                
                            //Fat Jet Index
                            FatJetIdx[PF].push_back(Pt[AK8].size()-1);
                        }
                    }
                }

                for(const reco::VertexCompositePtrCandidate& vtx: *secVtx){
                    if(BaseAnalyzer::DeltaR(vtx.p4().Eta(), vtx.p4().Phi(), Eta[AK8].back(), Phi[AK8].back()) < 0.8){
                        //SV four momentum components
                        Pt[VTX].push_back(vtx.p4().Pt());
                        Eta[VTX].push_back(vtx.p4().Eta());
                        Phi[VTX].push_back(vtx.p4().Phi());
                        Mass[VTX].push_back(vtx.p4().M());

                        //SV vertex
                        Vx[VTX].push_back(vtx.vx());
                        Vy[VTX].push_back(vtx.vy()); 
                        Vz[VTX].push_back(vtx.vz());
            
                        //Fat Jet Index
                        FatJetIdx[VTX].push_back(Pt[AK8].size()-1);
                    }
                }
            }
        }
    }

    //Loop over all jets
    for(int i = 0; i < jetSize; i++){
        //JER smearing
        smearFac = 1., corrFac = 1.;

        const float pt = isNANO ? jetPt->At(i) : jets->at(i).pt();
        const float eta = isNANO ? jetEta->At(i) : jets->at(i).eta();
        const float phi = isNANO ? jetPhi->At(i) : jets->at(i).phi();
        const float mass = isNANO ? jetMass->At(i) : jets->at(i).mass();

        //Define here already jet, because of smearing of 4-vec
        corrFac = isNANO ? CorrectEnergy(pt, eta,  *jetRho->Get(), jetArea->At(i), AK4) : CorrectEnergy(pt, eta, *rho, jets->at(i).jetArea(), AK4);

        //Get jet uncertainty
        if(jecUnc[AK4] !=  nullptr){
            jecUnc[AK4]->setJetPt(corrFac*pt);
            jecUnc[AK4]->setJetEta(eta);
            const float unc = jecUnc[AK4]->getUncertainty(isUp);
            corrFac *= isUp ? 1 + unc : 1 - unc;
        }

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(corrFac*pt, eta, phi, *jetRho->Get(), jetArea->At(i), AK4) : SmearEnergy(corrFac*pt, eta, phi, *rho, 0.4, AK4, *genJets);
        }
   
        if(smearFac*corrFac*pt > ptCut and abs(eta) < etaCut){
            FatIdx = -1.;

            //Check overlap with AK4 valid jets
            for(unsigned int j = 0; j < Pt[AK8].size(); j++){
                if(BaseAnalyzer::DeltaR(Eta[AK8][i], Phi[AK8][i], eta, phi) < 1.2){
                    FatIdx = j;
                }
            }

            JetType type = FatIdx == -1 ? AK4 : SUBAK4;

            if(type == SUBAK4) FatJetIdx[SUBAK4].push_back(FatIdx);

            //Jet four momentum components
            Pt[type].push_back(smearFac*corrFac*pt);
            Eta[type].push_back(eta);
            Phi[type].push_back(phi);
            Mass[type].push_back(smearFac*corrFac*mass);
            corrJEC[type].push_back(corrFac);
            if(!isData) corrJER[type].push_back(smearFac);

            //Correct met
            metPx += pt*std::cos(phi) - Pt[type].back()*std::cos(Phi[type].back());
            metPy += pt*std::sin(phi) - Pt[type].back()*std::sin(Phi[type].back());

            //Check for btag
            CSVBValue = 0, DeepBValue = 0;

            if(isNANO){ 
                DeepBValue = jetDeepBValue->At(i);
            } 

            else{
                for(std::string disc: {"pfDeepFlavourJetTags:probb", "pfDeepFlavourJetTags:probbb","pfDeepFlavourJetTags:problepb"}){
                    DeepBValue +=  jets->at(i).bDiscriminator(disc);
                }

                for(std::string disc: {"pfDeepCSVJetTags:probb", "pfDeepCSVJetTags:probbb"}){
                    CSVBValue +=  jets->at(i).bDiscriminator(disc);
                }
            }

            CSVScore[type].push_back(CSVBValue);
            DeepScore[type].push_back(DeepBValue);

            //True Flavour of Jet
            if(!isData){
                TrueFlavour[type].push_back(isNANO ? jetFlavour->At(i) : jets->at(i).partonFlavour());

                if(abs(isNANO ? jetFlavour->At(i) : jets->at(i).partonFlavour()) == 5){
                    bTotal->Fill(pt, eta);
                
                    if(DeepBValue > DeepLoose[era]) bTagEffLoose[0]->Fill(pt, eta);
                    if(DeepBValue > DeepMedium[era]) bTagEffMedium[0]->Fill(pt, eta);
                    if(DeepBValue > DeepTight[era]) bTagEffTight[0]->Fill(pt, eta);
                    if(CSVBValue > CSVLoose[era]) bTagEffLoose[1]->Fill(pt, eta);
                    if(CSVBValue > CSVMedium[era]) bTagEffMedium[1]->Fill(pt, eta);
                    if(CSVBValue > CSVTight[era]) bTagEffTight[1]->Fill(pt, eta);
                }
            }
    
            if(!isData){
                //btag SF
                looseCSVbTagSF[type].push_back(CSVReader.Get(Pt[type].back(), 0));
                mediumCSVbTagSF[type].push_back(CSVReader.Get(Pt[type].back(), 1));
                tightCSVbTagSF[type].push_back(CSVReader.Get(Pt[type].back(), 2));

                looseDeepbTagSF[type].push_back(DeepReader.Get(Pt[type].back(), 0));
                mediumDeepbTagSF[type].push_back(DeepReader.Get(Pt[type].back(), 1));
                tightDeepbTagSF[type].push_back(DeepReader.Get(Pt[type].back(), 2));

                if(!isSyst){
                    looseCSVbTagSFUp[type].push_back(CSVReader.GetUp(Pt[type].back(), 0));
                    looseCSVbTagSFDown[type].push_back(CSVReader.GetDown(Pt[type].back(), 0));
                    mediumCSVbTagSFUp[type].push_back(CSVReader.GetUp(Pt[type].back(), 1));
                    mediumCSVbTagSFDown[type].push_back(CSVReader.GetDown(Pt[type].back(), 1));
                    tightCSVbTagSFUp[type].push_back(CSVReader.GetUp(Pt[type].back(), 2));
                    tightCSVbTagSFDown[type].push_back(CSVReader.GetDown(Pt[type].back(), 2));

                    looseDeepbTagSFUp[type].push_back(DeepReader.GetUp(Pt[type].back(), 0));
                    looseDeepbTagSFDown[type].push_back(DeepReader.GetDown(Pt[type].back(), 0));
                    mediumDeepbTagSFUp[type].push_back(DeepReader.GetUp(Pt[type].back(), 1));
                    mediumDeepbTagSFDown[type].push_back(DeepReader.GetDown(Pt[type].back(), 1));
                    tightDeepbTagSFUp[type].push_back(DeepReader.GetUp(Pt[type].back(), 2));
                    tightDeepbTagSFDown[type].push_back(DeepReader.GetDown(Pt[type].back(), 2));
                }

                if(isNANO) SetGenParticles(i, Pt[type].back(), Eta[type].back(), Phi[type].back(), {5}, AK4);
                else{
                    event->getByToken(genParticleToken, genParts);
                    SetGenParticles(i, Pt[type].back(), Eta[type].back(), Phi[type].back(), {5}, AK4, *genParts);
                }
            }
        } 
    }

    nJets = Pt[AK4].size(); nSubJets = Pt[SUBAK4].size(); nFatJets = Pt[AK8].size();
    metPT = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
    metPHI = std::atan2(metPy, metPx);

    //Resort because of change PT order due to JEC
    std::map<std::string, JetType> jetNames = {{"Jet", AK4}, {"SubJet", SUBAK4}, {"FatJet", AK8}};

    for(std::pair<const std::string, JetType>& jet : jetNames){
        std::vector<int> idx(Pt[jet.second].size());
        std::iota(idx.begin(), idx.end(), 0);

        std::stable_sort(idx.begin(), idx.end(), [&](int i1, int i2) {return Pt[jet.second][i1] > Pt[jet.second][i2];});

        for(std::pair<const std::pair<std::string, std::string>, std::vector<float>&>& var: floatVar){
            if(var.second.size() == 0) continue;

            if(jet.first == var.first.first){
                std::vector<float> tmp(var.second.size());

                for(unsigned int l = 0; l < idx.size(); l++){
                    tmp[l] = var.second.at(idx[l]);
                }

                var.second = std::move(tmp);
            }
        }

        for(std::pair<const std::pair<std::string, std::string>, std::vector<char>&>& var: intVar){
            if(var.second.size() == 0) continue;

            if(jet.first == var.first.first){
                std::vector<char> tmp(var.second.size());

                for(unsigned int l = 0; l < idx.size(); l++){
                    tmp[l] = var.second.at(idx[l]);
                }

                var.second = std::move(tmp);
            }
        }
    }

    for(CutFlow& cutflow: cutflows){
        //Check if one combination of jet and fatjet number is fullfilled
        if(Pt[AK4].size() >= cutflow.nMinJet and Pt[AK8].size() >= cutflow.nMinFatjet){
            if(cutflow.passed and (cutflow.nMinJet!=0 or cutflow.nMinFatjet!=0)){
                std::string cutName("N^{AK4}_{jet} >= " + std::to_string(cutflow.nMinJet) + " && N^{AK8}_{jet} >= " + std::to_string(cutflow.nMinFatjet));

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

    if(!isData){
        for(int i = 0; i < 2; i++){
            for(const std::vector<TH2F*> eff : {bTagEffLoose, bTagEffMedium, bTagEffTight}){
                eff[i]->Write();
                delete eff[i];
            }
        }

        bTotal->Write();
        delete bTotal;
    }
}
