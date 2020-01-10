#include <ChargedSkimming/Skimming/interface/electronanalyzer.h>

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, eToken& eleToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken, const std::string& systematic):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    eleToken(eleToken),
    triggerObjToken(triggerObjToken),
    genParticleToken(genParticleToken)
    {
        energyCorrection = systematic == "" ? "ecalTrkEnergyPostCorr" : systematic; 
    }

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

void ElectronAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    //SF files
    mediumSFfiles = {
                    {2017, filePath + "eleSF/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp80iso.root"},
    };

    tightSFfiles = {
                    {2017, filePath + "eleSF/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp90iso.root"},
    };

    recoSFfiles = {
                    {2017, filePath + "eleSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"},
    };

    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    //Hist with scale factors
    TFile* recoSFfile = TFile::Open(recoSFfiles[era].c_str());
    recoSFhist = (TH2F*)recoSFfile->Get("EGamma_SF2D");

    TFile* mediumSFfile = TFile::Open(mediumSFfiles[era].c_str());
    mediumSFhist = (TH2F*)mediumSFfile->Get("EGamma_SF2D");

    TFile* tightSFfile = TFile::Open(tightSFfiles[era].c_str());
    tightSFhist = (TH2F*)tightSFfile->Get("EGamma_SF2D");

    //Initiliaze TTreeReaderValues then using NANO AOD
    if(isNANO){
        elePt = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pt");
        eleEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_eta");
        elePhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_phi");
        eleCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_charge");
        eleIso = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pfRelIso03_all");
        eleMediumMVA = std::make_unique<TTreeReaderArray<bool>>(*reader, "Electron_mvaFall17V2Iso_WP80");
        eleTightMVA = std::make_unique<TTreeReaderArray<bool>>(*reader, "Electron_mvaFall17V2Iso_WP90");

        //Set TTreeReader for genpart and trigger obj from baseanalyzer    
        SetCollection(this->isData);
    }

    //Set output names
    floatNames = {"E", "Px", "Py", "Pz", "Isolation", "Charge", "recoSF", "mediumSF", "tightSF"};

    if(!isSyst){
        std::vector<std::string> SFvariations = {"recoSFUp", "recoSFDown", "mediumSFUp", "mediumSFDown", "tightSFUp", "tightSFDown"};

        floatNames.insert(floatNames.end(), SFvariations.begin(), SFvariations.end());
    }

    boolNames = {"isMedium", "isTight", "isTriggerMatched", "isFromHc"};

    floatVariables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>());
    boolVariables = std::vector<std::vector<bool>>(boolNames.size(), std::vector<bool>());

    //Set Branches of output tree
    for(TTree* tree: trees){
        for(unsigned int i=0; i<floatVariables.size(); i++){
            tree->Branch(("Electron_" + floatNames[i]).c_str(), &floatVariables[i]);
        }

        for(unsigned int i=0; i<boolVariables.size(); i++){
            tree->Branch(("Electron_" + boolNames[i]).c_str(), &boolVariables[i]);
        }
    }
}

void ElectronAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::vector<float>& variable: floatVariables){
        variable.clear();
    }

    for(std::vector<bool>& variable: boolVariables){
        variable.clear();
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Electron>> electrons;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigObjects;
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isNANO){
        event->getByToken(eleToken, electrons);
        event->getByToken(triggerObjToken, trigObjects);
    }

    float eleSize = isNANO ? elePt->GetSize() : electrons->size();
    
    //Loop over all electrons
    for(unsigned int i = 0; i < eleSize; i++){
        float pt = isNANO ? elePt->At(i) : (electrons->at(i).p4()*electrons->at(i).userFloat(energyCorrection) / electrons->at(i).energy()).Pt();
        float eta = isNANO ? eleEta->At(i) : (electrons->at(i).p4()*electrons->at(i).userFloat(energyCorrection) / electrons->at(i).energy()).Eta();
        float phi = isNANO ? elePhi->At(i) : (electrons->at(i).p4()*electrons->at(i).userFloat(energyCorrection) / electrons->at(i).energy()).Phi();

        if(pt > ptCut && abs(eta) < etaCut){
            TLorentzVector lVec;
            lVec.SetPtEtaPhiM(pt, eta, phi, 0.510*1e-3);

            //Electron four momentum components
            floatVariables[0].push_back(lVec.E());   //Energy
            floatVariables[1].push_back(lVec.Px());  //Px
            floatVariables[2].push_back(lVec.Py());  //Py
            floatVariables[3].push_back(lVec.Pz());  //Pz

            //PF isolation
            float iso;

            if(!isNANO){
                iso = (electrons->at(i).pfIsolationVariables().sumChargedHadronPt + std::max(electrons->at(i).pfIsolationVariables().sumNeutralHadronEt + electrons->at(i).pfIsolationVariables().sumPhotonEt - 0.5 * electrons->at(i).pfIsolationVariables().sumPUPt, 0.0)) / electrons->at(i).pt();
            }

            floatVariables[4].push_back(isNANO ? eleIso->At(i) : iso); //PF isolation
            floatVariables[5].push_back(isNANO ? eleCharge->At(i) : electrons->at(i).charge());  //charge

            //Electron ID
            boolVariables[0].push_back(isNANO ? eleMediumMVA->At(i) : electrons->at(i).electronID("mvaEleID-Fall17-iso-V2-wp80"));  //Medium MVA ID
            boolVariables[1].push_back(isNANO ? eleMediumMVA->At(i) : electrons->at(i).electronID("mvaEleID-Fall17-iso-V2-wp90"));  //Tight MVA ID
            boolVariables[2].push_back(isNANO ? triggerMatching(lVec) : triggerMatching(lVec, *trigObjects)); //Trigger matching

            if(!isData){
               //Fill scale factors
                Int_t recoBin = recoSFhist->FindBin(eta, pt);
                Int_t mediumBin = mediumSFhist->FindBin(eta, pt);
                Int_t tightBin = tightSFhist->FindBin(eta, pt);

                floatVariables[6].push_back(recoSFhist->GetBinContent(recoBin));
                floatVariables[7].push_back(mediumSFhist->GetBinContent(mediumBin));
                floatVariables[8].push_back(tightSFhist->GetBinContent(tightBin));

                if(!isSyst){
                    floatVariables[9].push_back(recoSFhist->GetBinContent(recoBin) + recoSFhist->GetBinErrorUp(recoBin));
                    floatVariables[10].push_back(recoSFhist->GetBinContent(recoBin) - recoSFhist->GetBinErrorLow(recoBin));

                    floatVariables[11].push_back(mediumSFhist->GetBinContent(mediumBin) + mediumSFhist->GetBinErrorUp(mediumBin));
                    floatVariables[12].push_back(mediumSFhist->GetBinContent(mediumBin) - mediumSFhist->GetBinErrorLow(mediumBin));

                    floatVariables[13].push_back(tightSFhist->GetBinContent(tightBin) + tightSFhist->GetBinErrorUp(tightBin));
                    floatVariables[14].push_back(tightSFhist->GetBinContent(tightBin) - tightSFhist->GetBinErrorLow(tightBin));
                }

                //Save gen particle information
                if(isNANO) boolVariables[3].push_back(SetGenParticles(lVec, i, 11));
                else{
                    event->getByToken(genParticleToken, genParts);
                    boolVariables[3].push_back(SetGenParticles(lVec, i, 11, *genParts));
                }
            }
        } 
    }

    for(CutFlow &cutflow: cutflows){
        if(cutflow.nMinEle <= floatVariables[0].size()){
            if(cutflow.nMinEle!=0 and cutflow.passed){
                std::string cutName("N_{e} >= " + std::to_string(cutflow.nMinEle) + " (no iso/ID req)");
                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
            }
        }

        else{
            cutflow.passed = false;
        }
    }
}

void ElectronAnalyzer::EndJob(TFile* file){
}
