#include <ChargedHiggs/Skimming/interface/muonanalyzer.h>

MuonAnalyzer::MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNMuon, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNMuon(minNMuon)
    {}

MuonAnalyzer::MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNMuon, muToken& muonToken):
    BaseAnalyzer(), 
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNMuon(minNMuon),
    muonToken(muonToken)
    {}


void MuonAnalyzer::SetGenParticles(Muon &validMuon, const int &i){
    //Check if gen matched particle exist
    if(muonGenIdx->At(i) != -1){
        validMuon.isgenMatched = true;        
        int idxMotherMu = genMotherIdx->At(muonGenIdx->At(i));

        while(abs(genID->At(idxMotherMu)) == 13){
            idxMotherMu = genMotherIdx->At(idxMotherMu);
        }

        validMuon.genVec.SetPtEtaPhiM(genPt->At(muonGenIdx->At(i)), genEta->At(muonGenIdx->At(i)), genPhi->At(muonGenIdx->At(i)), genMass->At(muonGenIdx->At(i)));

        if(abs(genID->At(idxMotherMu)) == 24){
            float idxMotherW = genMotherIdx->At(idxMotherMu);

            while(abs(genID->At(idxMotherW)) == 24){
                idxMotherW = genMotherIdx->At(idxMotherW);
            }

            if(abs(genID->At(idxMotherW)) == 37){
                validMuon.isFromHc = true;
            }
        }
    }
}

void MuonAnalyzer::BeginJob(TTree* tree, bool &isData){
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

    //Set Branches of output tree
    tree->Branch("muon", &validMuons);
}

bool MuonAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event){
    //Clear electron vector3
    validMuons.clear();

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Muon>> muons;

    if(!isNANO){
        event->getByToken(muonToken, muons);
    }

    float muSize = isNANO ? muonPt->GetSize() : muons->size();

    //Loop over all electrons
    for(unsigned int i = 0; i < muSize; i++){
        float pt = isNANO ? muonPt->At(i) : muons->at(i).pt();
        float eta = isNANO ? muonEta->At(i) : muons->at(i).eta();
        float phi = isNANO ? muonPhi->At(i) : muons->at(i).phi();

        if(pt > ptCut && abs(eta) < etaCut){
            Muon validMuon;

            //Set muon information
            validMuon.fourVec.SetPtEtaPhiM(pt, eta, phi, 105.658*1e-3);

            validMuon.isMedium = isNANO ? muonMediumID->At(i) : muons->at(i).isMediumMuon();
            //validMuon.isTight = isNANO ? muonTightID->At(i) : muons->at(i).isTightMuon();
            validMuon.isLooseIso = isNANO ? muonIso->At(i) < 0.25 : muons->at(i).PFIsoLoose;
            validMuon.isTightIso = isNANO ? muonIso->At(i) < 0.15 : muons->at(i).PFIsoTight;
            validMuon.charge = isNANO ? muonCharge->At(i) : muons->at(i).charge();
            //validMuon.isTriggerMatched = triggerMatching(validMuon.fourVec, 13);
            
            if(!isData){
                validMuon.triggerSF = triggerSFhist->GetBinContent(triggerSFhist->FindBin(pt, abs(eta))) != 0 ? triggerSFhist->GetBinContent(triggerSFhist->FindBin(pt, abs(eta)))  : 1.;
                validMuon.mediumSF = mediumIDHist->GetBinContent(mediumIDHist->FindBin(pt, abs(eta))) != 0 ? mediumIDHist->GetBinContent(mediumIDHist->FindBin(pt, abs(eta))) : 1.;
                validMuon.tightSF = tightIDHist->GetBinContent(tightIDHist->FindBin(pt, abs(eta))) != 0 ? tightIDHist->GetBinContent(tightIDHist->FindBin(pt, abs(eta))) : 1.;
                validMuon.looseIsoMediumSF = looseIsoMediumSFhist->GetBinContent(looseIsoMediumSFhist->FindBin(pt, abs(eta))) != 0 ? looseIsoMediumSFhist->GetBinContent(looseIsoMediumSFhist->FindBin(pt, abs(eta))) : 1.;
                validMuon.tightIsoMediumSF = tightIsoMediumSFhist->GetBinContent(tightIsoMediumSFhist->FindBin(pt, abs(eta))) != 0 ? tightIsoMediumSFhist->GetBinContent(tightIsoMediumSFhist->FindBin(pt, abs(eta))) : 1.;
                validMuon.looseIsoTightSF = looseIsoTightSFhist->GetBinContent(looseIsoTightSFhist->FindBin(pt, abs(eta))) != 0 ? looseIsoTightSFhist->GetBinContent(looseIsoTightSFhist->FindBin(pt, abs(eta))) : 1.;
                validMuon.tightIsoTightSF = tightIsoTightSFhist->GetBinContent(tightIsoTightSFhist->FindBin(pt, abs(eta))) != 0 ? tightIsoTightSFhist->GetBinContent(tightIsoTightSFhist->FindBin(pt, abs(eta))) : 1.;

                //Save gen particle information
                //SetGenParticles(validMuon, i);
             }

            //Fill electron in collection
            validMuons.push_back(validMuon);       
        } 
    }
    
    //Check if event has enough electrons
    if(validMuons.size() < minNMuon){
        return false;
    }

    if(minNMuon != 0){
        std::string cutName("N_{#mu} >= " + std::to_string(minNMuon) + " (no iso/ID req.)");
        cutflow.first->Fill(cutName.c_str(), cutflow.second);
    }
    return true;
}

void MuonAnalyzer::EndJob(TFile* file){
}
