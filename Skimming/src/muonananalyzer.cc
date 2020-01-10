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

void MuonAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
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
    this->isSyst = isSyst;

    //Hist with scale factors
    TFile* triggerSFfile = TFile::Open(triggerSFfiles[era].c_str());
    triggerSFhist = (TH2F*)triggerSFfile->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");

    TFile* isoSFfile = TFile::Open(isoSFfiles[era].c_str());
    IsoHist.push_back((TH2F*)isoSFfile->Get("NUM_LooseRelIso_DEN_LooseID_pt_abseta"));
    IsoHist.push_back((TH2F*)isoSFfile->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"));

    TFile* IDSFfile = TFile::Open(IDSFfiles[era].c_str());
    IDHist.push_back((TH2F*)IDSFfile->Get("NUM_LooseID_DEN_genTracks_pt_abseta"));
    IDHist.push_back((TH2F*)IDSFfile->Get("NUM_TightID_DEN_genTracks_pt_abseta"));

    if(isNANO){
        //Initiliaze TTreeReaderValues
        muonPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_pt");
        muonEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_eta");
        muonPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_phi");
        muonCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_charge");
        muonIso = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_miniPFRelIso_all");
        muonLooseID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_looseId");
        muonTightID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_tightId");

        if(!this->isData){
            muonGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_genPartIdx");
        }


        //Set TTreeReader for genpart and trigger obj from baseanalyzer
        SetCollection(this->isData);
    }

    //Set output names
    floatNames = {"E", "Px", "Py", "Pz", "Charge", "looseIsoLooseSF", "tightIsoTightSF", "looseSF", "tightSF", "triggerSF"};

    if(!isSyst){
        std::vector<std::string> SFvariations = {"looseIsoLooseSFUp", "looseIsoLooseSFDown", "tightIsoTightSFUp", "tightIsoTightSFDown", "looseSFUp", "looseSDown", "tightSFUp", "tightSFDown", "triggerSFUp", "triggerSFDown"}; 

        floatNames.insert(floatNames.end(), SFvariations.begin(), SFvariations.end());   
    }

    boolNames = {"isLooseIso", "isTightIso", "isLoose", "isTight", "isTriggerMatched", "isFromHc"};

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
            boolVariables[0].push_back(isNANO ? muonIso->At(i) < 0.25 : muons->at(i).passed(reco::Muon::PFIsoLoose));
            boolVariables[1].push_back(isNANO ? muonIso->At(i) < 0.15 : muons->at(i).passed(reco::Muon::PFIsoTight));

            boolVariables[2].push_back(isNANO ? muonLooseID->At(i) : muons->at(i).passed(reco::Muon::CutBasedIdLoose));
            boolVariables[3].push_back(isNANO ? muonTightID->At(i) : muons->at(i).passed(reco::Muon::CutBasedIdTight));
            boolVariables[4].push_back(isNANO ? triggerMatching(lVec) : triggerMatching(lVec, *trigObjects));
            
            if(!isData){
                //Scale factors
                Int_t looseIsoBin = IsoHist[0]->FindBin(pt, abs(eta));
                Int_t tightIsoBin = IsoHist[1]->FindBin(pt, abs(eta));
                Int_t looseIDBin = IDHist[0]->FindBin(pt, abs(eta));
                Int_t tightIDBin = IDHist[1]->FindBin(pt, abs(eta));
                Int_t triggerBin = triggerSFhist->FindBin(pt, abs(eta));

                floatVariables[5].push_back(IsoHist[0]->GetBinContent(looseIsoBin));
                floatVariables[6].push_back(IsoHist[1]->GetBinContent(tightIsoBin));
                floatVariables[7].push_back(IDHist[0]->GetBinContent(looseIDBin));
                floatVariables[8].push_back(IDHist[1]->GetBinContent(tightIDBin));
                floatVariables[9].push_back(triggerSFhist->GetBinContent(triggerBin));

                if(!isSyst){
                    floatVariables[10].push_back(IsoHist[0]->GetBinContent(looseIsoBin) + IsoHist[0]->GetBinErrorUp(looseIsoBin));
                    floatVariables[11].push_back(IsoHist[0]->GetBinContent(looseIsoBin) - IsoHist[0]->GetBinErrorLow(looseIsoBin));

                    floatVariables[12].push_back(IsoHist[1]->GetBinContent(tightIsoBin) + IsoHist[1]->GetBinErrorUp(tightIsoBin));
                    floatVariables[13].push_back(IsoHist[1]->GetBinContent(tightIsoBin) - IsoHist[1]->GetBinErrorLow(tightIsoBin));

                    floatVariables[14].push_back(IDHist[0]->GetBinContent(looseIDBin) + IDHist[0]->GetBinErrorUp(looseIDBin));
                    floatVariables[15].push_back(IDHist[0]->GetBinContent(looseIDBin) - IDHist[0]->GetBinErrorLow(looseIDBin));

                    floatVariables[16].push_back(IDHist[1]->GetBinContent(tightIDBin) + IDHist[1]->GetBinErrorUp(tightIDBin));
                    floatVariables[17].push_back(IDHist[1]->GetBinContent(tightIDBin) - IDHist[1]->GetBinErrorLow(tightIDBin));

                    floatVariables[18].push_back(triggerSFhist->GetBinContent(triggerBin) + triggerSFhist->GetBinErrorUp(triggerBin));
                    floatVariables[19].push_back(triggerSFhist->GetBinContent(triggerBin) - triggerSFhist->GetBinErrorLow(triggerBin));
                }

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
