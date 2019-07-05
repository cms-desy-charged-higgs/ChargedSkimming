#include <ChargedHiggs/Skimming/interface/muonanalyzer.h>

MuonAnalyzer::MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNMuon, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNMuon(minNMuon)
    {}

MuonAnalyzer::MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNMuon, edm::Handle<std::vector<pat::Muon>> muons):
    BaseAnalyzer(), 
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNMuon(minNMuon),
    muons(muons)
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

    //Set Branches of output tree
    tree->Branch("muon", &validMuons);
}

bool MuonAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow){
    //Clear electron vector3
    validMuons.clear();

    //Loop over all electrons
    for(unsigned int i = 0; i < muonPt->GetSize(); i++){
        if(muonPt->At(i) > ptCut && abs(muonEta->At(i)) < etaCut){
            Muon validMuon;

            //Set muon information
            validMuon.fourVec.SetPtEtaPhiM(muonPt->At(i), muonEta->At(i), muonPhi->At(i), 105.658*1e-3);

            validMuon.isMedium = muonMediumID->At(i);
            validMuon.isTight = muonTightID->At(i);
            validMuon.isLooseIso = muonIso->At(i) < 0.25;
            validMuon.isTightIso = muonIso->At(i) < 0.15;
            validMuon.charge = muonCharge->At(i);
            validMuon.isTriggerMatched = triggerMatching(validMuon.fourVec, 13);
            
            if(!isData){
                validMuon.triggerSF = triggerSFhist->GetBinContent(triggerSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? triggerSFhist->GetBinContent(triggerSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i))))  : 1.;
                validMuon.mediumSF = mediumIDHist->GetBinContent(mediumIDHist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? mediumIDHist->GetBinContent(mediumIDHist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) : 1.;
                validMuon.tightSF = tightIDHist->GetBinContent(tightIDHist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? tightIDHist->GetBinContent(tightIDHist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) : 1.;
                validMuon.looseIsoMediumSF = looseIsoMediumSFhist->GetBinContent(looseIsoMediumSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? looseIsoMediumSFhist->GetBinContent(looseIsoMediumSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) : 1.;
                validMuon.tightIsoMediumSF = tightIsoMediumSFhist->GetBinContent(tightIsoMediumSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? tightIsoMediumSFhist->GetBinContent(tightIsoMediumSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) : 1.;
                validMuon.looseIsoTightSF = looseIsoTightSFhist->GetBinContent(looseIsoTightSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? looseIsoTightSFhist->GetBinContent(looseIsoTightSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) : 1.;
                validMuon.tightIsoTightSF = tightIsoTightSFhist->GetBinContent(tightIsoTightSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) != 0 ? tightIsoTightSFhist->GetBinContent(tightIsoTightSFhist->FindBin(muonPt->At(i), abs(muonEta->At(i)))) : 1.;

                //Save gen particle information
                SetGenParticles(validMuon, i);
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
