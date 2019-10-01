#include <ChargedSkimming/Skimming/interface/muonanalyzer.h>

MuonAnalyzer::MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

MuonAnalyzer::MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, muToken& muonToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken):
    BaseAnalyzer(), 
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    muonToken(muonToken),
    triggerObjToken(triggerObjToken),
    genParticleToken(genParticleToken)
    {}

void MuonAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData){
    isoSFfiles = {
        {2017, filePath + "/muonSF/RunBCDEF_SF_ISO.root"},
    };

    triggerSFfiles = {
        {2017, filePath + "/muonSF/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root"},
    };

    IDSFfiles = {
        {2017, filePath + "/muonSF/RunBCDEF_SF_ID.root"},
    };

    //Set data bool
    this->isData = isData;

    //Hist with scale factors
    TFile* triggerSFfile = TFile::Open(triggerSFfiles[era].c_str());
    triggerSFhist = (TH2F*)triggerSFfile->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");

    TFile* isoSFfile = TFile::Open(isoSFfiles[era].c_str());
    looseIsoMediumSFhist = (TH2F*)isoSFfile->Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta");
    tightIsoMediumSFhist = (TH2F*)isoSFfile->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta");
    looseIsoTightSFhist = (TH2F*)isoSFfile->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_pt_abseta");
    tightIsoTightSFhist = (TH2F*)isoSFfile->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

    TFile* IDSFfile = TFile::Open(IDSFfiles[era].c_str());
    mediumIDHist = (TH2F*)IDSFfile->Get("NUM_MediumID_DEN_genTracks_pt_abseta");
    tightIDHist = (TH2F*)IDSFfile->Get("NUM_TightID_DEN_genTracks_pt_abseta");

    if(isNANO){
        //Initiliaze TTreeReaderValues
        muonPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_pt");
        muonEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_eta");
        muonPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_phi");
        muonCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_charge");
        muonIso = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_miniPFRelIso_all");
        muonMediumID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_mediumId");
        muonTightID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_tightId");

        if(!this->isData){
            muonGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_genPartIdx");
        }


        //Set TTreeReader for genpart and trigger obj from baseanalyzer
        SetCollection(this->isData);
    }

    //Set output names
    floatNames = {"E", "Px", "Py", "Pz", "Charge", "looseIsoMediumSF", "tightIsoMediumSF", "looseIsoTightSF", "tightIsoTightSF", "mediumSF", "tightSF", "triggerSF"};
    boolNames = {"isMediumIso", "isTightIso", "isMedium", "isTight", "isTriggerMatched", "isFromHc"};

    floatVariables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>());
    boolVariables = std::vector<std::vector<bool>>(boolNames.size(), std::vector<bool>());

    //Set Branches of output tree
    for(TTree* tree: trees){
        for(unsigned int i=0; i<floatVariables.size(); i++){
            tree->Branch(("Muon_" + floatNames[i]).c_str(), &floatVariables[i]);
        }

        for(unsigned int i=0; i<boolVariables.size(); i++){
            tree->Branch(("Muon_" + boolNames[i]).c_str(), &boolVariables[i]);
        }
    }
}

void MuonAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::vector<float>& variable: floatVariables){
        variable.clear();
    }

    for(std::vector<bool>& variable: boolVariables){
        variable.clear();
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Muon>> muons;
    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigObjects;
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isNANO){
        event->getByToken(muonToken, muons);
        event->getByToken(triggerObjToken, trigObjects);
    }

    float muSize = isNANO ? muonPt->GetSize() : muons->size();

    //Loop over all electrons
    for(unsigned int i = 0; i < muSize; i++){
        float pt = isNANO ? muonPt->At(i) : muons->at(i).pt();
        float eta = isNANO ? muonEta->At(i) : muons->at(i).eta();
        float phi = isNANO ? muonPhi->At(i) : muons->at(i).phi();

        if(pt > ptCut && abs(eta) < etaCut){
            TLorentzVector lVec;
            lVec.SetPtEtaPhiM(pt, eta, phi, 105.658*1e-3);

            //Muon four momentum components
            floatVariables[0].push_back(lVec.E());   //Energy
            floatVariables[1].push_back(lVec.Px());  //Px
            floatVariables[2].push_back(lVec.Py());  //Py
            floatVariables[3].push_back(lVec.Pz());  //Pz

            floatVariables[4].push_back(isNANO ? muonCharge->At(i) : muons->at(i).charge());

            //Isolation and ID
            boolVariables[0].push_back(isNANO ? muonIso->At(i) < 0.25 : muons->at(i).PFIsoLoose);
            boolVariables[1].push_back(isNANO ? muonIso->At(i) < 0.15 : muons->at(i).PFIsoTight);
            boolVariables[2].push_back(isNANO ? muonMediumID->At(i) : muons->at(i).CutBasedIdMedium);
            boolVariables[3].push_back(isNANO ? muonTightID->At(i) : muons->at(i).CutBasedIdTight);
            boolVariables[4].push_back(isNANO ? triggerMatching(lVec) : triggerMatching(lVec, *trigObjects));
            
            if(!isData){
                //Scale factors
                floatVariables[5].push_back(triggerSFhist->GetBinContent(triggerSFhist->FindBin(pt, abs(eta))) != 0 ? triggerSFhist->GetBinContent(triggerSFhist->FindBin(pt, abs(eta))) : 1.);
                floatVariables[6].push_back(mediumIDHist->GetBinContent(mediumIDHist->FindBin(pt, abs(eta))) != 0 ? mediumIDHist->GetBinContent(mediumIDHist->FindBin(pt, abs(eta))) : 1.);
                floatVariables[7].push_back(tightIDHist->GetBinContent(tightIDHist->FindBin(pt, abs(eta))) != 0 ? tightIDHist->GetBinContent(tightIDHist->FindBin(pt, abs(eta))) : 1.);
                floatVariables[8].push_back(looseIsoMediumSFhist->GetBinContent(looseIsoMediumSFhist->FindBin(pt, abs(eta))) != 0 ? looseIsoMediumSFhist->GetBinContent(looseIsoMediumSFhist->FindBin(pt, abs(eta))) : 1.);
                floatVariables[9].push_back(tightIsoMediumSFhist->GetBinContent(tightIsoMediumSFhist->FindBin(pt, abs(eta))) != 0 ? tightIsoMediumSFhist->GetBinContent(tightIsoMediumSFhist->FindBin(pt, abs(eta))) : 1.);
                floatVariables[10].push_back(looseIsoTightSFhist->GetBinContent(looseIsoTightSFhist->FindBin(pt, abs(eta))) != 0 ? looseIsoTightSFhist->GetBinContent(looseIsoTightSFhist->FindBin(pt, abs(eta))) : 1.);
                floatVariables[11].push_back(tightIsoTightSFhist->GetBinContent(tightIsoTightSFhist->FindBin(pt, abs(eta))) != 0 ? tightIsoTightSFhist->GetBinContent(tightIsoTightSFhist->FindBin(pt, abs(eta))) : 1.);

                //Save gen particle information
                if(isNANO) boolVariables[5].push_back(SetGenParticles(lVec, i, 13));
                else{
                    event->getByToken(genParticleToken, genParts);
                    boolVariables[5].push_back(SetGenParticles(lVec, i, 13, *genParts));
                }
             }
        } 
    }
    
    //Check if event has enough electrons
    for(CutFlow &cutflow: cutflows){
        if(cutflow.nMinMu <= floatVariables[0].size()){
            if(cutflow.nMinMu!=0 and cutflow.passed){
                std::string cutName("N_{#mu} >= " + std::to_string(cutflow.nMinMu) + " (no iso/ID req)");
                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
            }
        }

        else{
            cutflow.passed = false;
        }
    }
}

void MuonAnalyzer::EndJob(TFile* file){
}
