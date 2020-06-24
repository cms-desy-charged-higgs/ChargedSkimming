#include <ChargedSkimming/Skimming/interface/muonanalyzer.h>

MuonAnalyzer::MuonAnalyzer(const int& era, const float& ptCut, const float& etaCut, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

MuonAnalyzer::MuonAnalyzer(const int& era, const float& ptCut, const float& etaCut, muToken& muonToken, genPartToken& genParticleToken):
    BaseAnalyzer(), 
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    muonToken(muonToken),
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
    IsoHist.push_back((TH2F*)isoSFfile->Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta"));
    IsoHist.push_back((TH2F*)isoSFfile->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_pt_abseta"));
    IsoHist.push_back((TH2F*)isoSFfile->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta"));
    IsoHist.push_back((TH2F*)isoSFfile->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"));

    TFile* IDSFfile = TFile::Open(IDSFfiles[era].c_str());
    IDHist.push_back((TH2F*)IDSFfile->Get("NUM_LooseID_DEN_genTracks_pt_abseta"));
    IDHist.push_back((TH2F*)IDSFfile->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
    IDHist.push_back((TH2F*)IDSFfile->Get("NUM_TightID_DEN_genTracks_pt_abseta"));

    if(isNANO){
        //Initiliaze TTreeReaderValues
        muonPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_pt");
        muonEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_eta");
        muonPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_phi");
        muonCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_charge");
        muonIso = std::make_unique<TTreeReaderArray<float>>(*reader, "Muon_pfRelIso04_all");
        muonLooseID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_looseId");
        muonMediumID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_mediumId");
        muonTightID = std::make_unique<TTreeReaderArray<bool>>(*reader, "Muon_tightId");

        if(!this->isData){
            muonGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Muon_genPartIdx");
        }

        //Set TTreeReader for genpart and trigger obj from baseanalyzer
        SetCollection(this->isData);
    }

    //Variable name mapping to branch name
    floatVar = {
            {"Pt", Pt},
            {"Eta", Eta},
            {"Phi", Phi},
            {"triggerSF", triggerSF},
            {"looseSF", looseSF},
            {"mediumSF", mediumSF},
            {"tightSF", tightSF},
            {"looseIsoLooseSF", looseIsoLooseSF},
            {"looseIsoMediumSF", looseIsoMediumSF},
            {"looseIsoTightSF", looseIsoTightSF},
            {"tightIsoMediumSF", tightIsoMediumSF},
            {"tightIsoTightSF", tightIsoTightSF},
    };

    intVar = {
            {"Charge", Charge},
            {"ID", ID},
            {"isoID", isoID},
            {"ParticleID", partID},
            {"MotherID", mothID},
            {"GrandMotherID", grandID},
    };

    if(!isSyst){
        std::map<std::string, std::vector<float>&> SFvariations = {
            {"triggerSFUp", triggerSFUp},
            {"triggerSFDown", triggerSFDown},
            {"looseSFUp", looseSFUp},
            {"looseSFDown", looseSFDown},
            {"mediumSFUp", mediumSFUp},
            {"mediumSFDown", mediumSFDown},
            {"tightSFUp", tightSFUp},
            {"tightSFDown", tightSFDown},
            {"looseIsoLooseSFDown", looseIsoLooseSFDown},
            {"looseIsoLooseSFUp", looseIsoLooseSFUp},
            {"looseIsoMediumSFDown", looseIsoMediumSFDown},
            {"looseIsoMediumSFUp", looseIsoMediumSFUp},
            {"looseIsoTightSFDown", looseIsoTightSFDown},
            {"looseIsoTightSFUp", looseIsoTightSFUp},
            {"tightIsoMediumSFDown", tightIsoMediumSFDown},
            {"tightIsoMediumSFUp", tightIsoMediumSFUp},
            {"tightIsoTightSFDown", tightIsoTightSFDown},
            {"tightIsoTightSFUp", tightIsoTightSFUp}
        };

        floatVar.insert(SFvariations.begin(), SFvariations.end());   
    }

    //Set Branches of output tree
    for(TTree* tree: trees){
        tree->Branch("Muon_Size", &nMuons);

        for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
            tree->Branch(("Muon_" + var.first).c_str(), &var.second);
        }

        for(std::pair<const std::string, std::vector<char>&>& var: intVar){
            tree->Branch(("Muon_" + var.first).c_str(), &var.second);
        }
    }
}

void MuonAnalyzer::Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
        var.second.clear();
    }

    for(std::pair<const std::string, std::vector<char>&>& var: intVar){
        var.second.clear();
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Muon>> muons;
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isNANO){
        event->getByToken(muonToken, muons);
    }

    const char& muSize = isNANO ? muonPt->GetSize() : muons->size();

    //Loop over all electrons
    for(int i = 0; i < muSize; i++){
        const float& pt = isNANO ? muonPt->At(i) : muons->at(i).pt();
        const float& eta = isNANO ? muonEta->At(i) : muons->at(i).eta();
        const float& phi = isNANO ? muonPhi->At(i) : muons->at(i).phi();

        if(pt > ptCut && abs(eta) < etaCut){
            //Muon four momentum components
            Pt.push_back(pt);
            Eta.push_back(eta);
            Phi.push_back(phi);

            //ID
            if(isNANO ? muonTightID->At(i) : muons->at(i).passed(reco::Muon::MvaTight)) ID.push_back(3);
            else if(isNANO ? muonMediumID->At(i) : muons->at(i).passed(reco::Muon::MvaMedium)) ID.push_back(2);
            else if(isNANO ? muonLooseID->At(i) : muons->at(i).passed(reco::Muon::MvaLoose)) ID.push_back(1);
            else ID.push_back(0);

            //Isolation
            if(isNANO ? muonIso->At(i) < 0.15 : muons->at(i).passed(reco::Muon::PFIsoTight)) isoID.push_back(3);
            else if(isNANO ? muonIso->At(i) < 0.2 : muons->at(i).passed(reco::Muon::PFIsoMedium)) isoID.push_back(2);
            else if(isNANO ? muonIso->At(i) < 0.25 : muons->at(i).passed(reco::Muon::PFIsoLoose)) isoID.push_back(1);
            else isoID.push_back(0);

            //Charge
            Charge.push_back(isNANO ? muonCharge->At(i) : muons->at(i).charge());
            
            if(!isData){
                //Scale factors
                const Int_t& looseIsoLooseIDBin = IsoHist[0]->FindBin(pt, abs(eta));
                const Int_t& looseIsoMediumIDBin = IsoHist[1]->FindBin(pt, abs(eta));
                const Int_t& looseIsoTightIDBin = IsoHist[2]->FindBin(pt, abs(eta));
                const Int_t& tightIsoMediumIDBin = IsoHist[3]->FindBin(pt, abs(eta));
                const Int_t& tightIsoTightIDBin = IsoHist[4]->FindBin(pt, abs(eta));

                const Int_t& looseIDBin = IDHist[0]->FindBin(pt, abs(eta));
                const Int_t& mediumIDBin = IDHist[1]->FindBin(pt, abs(eta));
                const Int_t& tightIDBin = IDHist[2]->FindBin(pt, abs(eta));
                const Int_t& triggerBin = triggerSFhist->FindBin(pt, abs(eta));

                looseIsoLooseSF.push_back(IsoHist[0]->GetBinContent(looseIsoLooseIDBin));
                looseIsoMediumSF.push_back(IsoHist[1]->GetBinContent(looseIsoMediumIDBin));
                looseIsoTightSF.push_back(IsoHist[2]->GetBinContent(looseIsoTightIDBin));
                tightIsoMediumSF.push_back(IsoHist[3]->GetBinContent(tightIsoMediumIDBin));
                tightIsoTightSF.push_back(IsoHist[4]->GetBinContent(tightIsoTightIDBin));
                looseSF.push_back(IDHist[0]->GetBinContent(looseIDBin));
                mediumSF.push_back(IDHist[1]->GetBinContent(mediumIDBin));
                tightSF.push_back(IDHist[2]->GetBinContent(triggerBin));
                triggerSF.push_back(triggerSFhist->GetBinContent(triggerBin));

                if(!isSyst){
                    looseIsoLooseSFUp.push_back(IsoHist[0]->GetBinContent(looseIsoLooseIDBin) + IsoHist[0]->GetBinErrorUp(looseIsoLooseIDBin));
                    looseIsoLooseSFDown.push_back(IsoHist[0]->GetBinContent(looseIsoLooseIDBin) - IsoHist[0]->GetBinErrorLow(looseIsoLooseIDBin));
                    looseIsoMediumSFUp.push_back(IsoHist[1]->GetBinContent(looseIsoMediumIDBin) + IsoHist[1]->GetBinErrorUp(looseIsoMediumIDBin));
                    looseIsoMediumSFDown.push_back(IsoHist[1]->GetBinContent(looseIsoMediumIDBin) - IsoHist[1]->GetBinErrorLow(looseIsoMediumIDBin));
                    looseIsoTightSFUp.push_back(IsoHist[2]->GetBinContent(looseIsoTightIDBin) + IsoHist[2]->GetBinErrorUp(looseIsoTightIDBin));
                    looseIsoTightSFDown.push_back(IsoHist[2]->GetBinContent(looseIsoTightIDBin) - IsoHist[2]->GetBinErrorLow(looseIsoTightIDBin));
                    tightIsoMediumSFUp.push_back(IsoHist[3]->GetBinContent(tightIsoMediumIDBin) + IsoHist[3]->GetBinErrorUp(tightIsoMediumIDBin));
                    tightIsoMediumSFDown.push_back(IsoHist[3]->GetBinContent(tightIsoMediumIDBin) - IsoHist[3]->GetBinErrorLow(tightIsoMediumIDBin));
                    tightIsoTightSFUp.push_back(IsoHist[4]->GetBinContent(looseIsoTightIDBin) + IsoHist[4]->GetBinErrorUp(looseIsoTightIDBin));
                    tightIsoTightSFDown.push_back(IsoHist[4]->GetBinContent(tightIsoTightIDBin) - IsoHist[4]->GetBinErrorLow(tightIsoTightIDBin));

                    looseSFUp.push_back(IDHist[0]->GetBinContent(looseIDBin) + IDHist[0]->GetBinErrorUp(looseIDBin));
                    looseSFDown.push_back(IDHist[0]->GetBinContent(looseIDBin) - IDHist[0]->GetBinErrorLow(looseIDBin));
                    mediumSFUp.push_back(IDHist[1]->GetBinContent(mediumIDBin) + IDHist[1]->GetBinErrorUp(mediumIDBin));
                    mediumSFDown.push_back(IDHist[1]->GetBinContent(mediumIDBin) - IDHist[1]->GetBinErrorLow(mediumIDBin));
                    tightSFUp.push_back(IDHist[2]->GetBinContent(tightIDBin) + IDHist[2]->GetBinErrorUp(tightIDBin));
                    tightSFDown.push_back(IDHist[2]->GetBinContent(tightIDBin) - IDHist[2]->GetBinErrorLow(tightIDBin));

                    triggerSFUp.push_back(triggerSFhist->GetBinContent(triggerBin) + triggerSFhist->GetBinErrorUp(triggerBin));
                    triggerSFDown.push_back(triggerSFhist->GetBinContent(triggerBin) - triggerSFhist->GetBinErrorLow(triggerBin));
                }

                //Save gen particle information
                std::tuple<int, int, int> IDs;

                if(isNANO){
                    IDs = SetGenParticles(pt, eta, phi, i, 13);
                    partID.push_back(std::get<0>(IDs));
                    mothID.push_back(std::get<1>(IDs));
                    grandID.push_back(std::get<2>(IDs));
                }
    
                else{
                    event->getByToken(genParticleToken, genParts);
                    IDs = SetGenParticles(pt, eta, phi, i, 13, *genParts);
                    partID.push_back(std::get<0>(IDs));
                    mothID.push_back(std::get<1>(IDs));
                    grandID.push_back(std::get<2>(IDs));
                }
             }
        }
    }

    nMuons = Pt.size();

    //Check if event has enough electrons
    for(CutFlow &cutflow: cutflows){
        if(cutflow.nMinMu <= Pt.size()){
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
