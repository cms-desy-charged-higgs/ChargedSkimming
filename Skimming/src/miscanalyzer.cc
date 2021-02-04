#include <ChargedSkimming/Skimming/interface/miscanalyzer.h>

MiscAnalyzer::MiscAnalyzer(const int& era, const float& etaCut, isoToken& isoTrackToken):
    BaseAnalyzer(), 
    era(era),
    etaCut(etaCut),
    isoTrackToken(isoTrackToken)
    {}

void MiscAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    //Variable name mapping to branch name
    floatVar = {
            {"Pt", Pt},
            {"Eta", Eta},
            {"Phi", Phi},
            {"Dxy", Dxy},
            {"Dz", Dz},
            {"Isolation", Isolation},
    };

    intVar = {
            {"Charge", Charge},
            {"PDG", PDG},
            {"FromPV", FromPV},
    };


    //Set Branches of output tree
    for(TTree* tree: trees){
        tree->Branch("IsoTrack_Size", &nTracks, "IsoTrack_Size/S");

        for(std::pair<const std::string, float*>& var: floatVar){
            tree->Branch(("IsoTrack_" + var.first).c_str(), var.second, ("IsoTrack_" + var.first + "[IsoTrack_Size]/F").c_str());
        }

        for(std::pair<const std::string, short*>& var: intVar){
            tree->Branch(("IsoTrack_" + var.first).c_str(), var.second, ("IsoTrack_" + var.first + "[IsoTrack_Size]/S").c_str());
        }
    }
}

void MiscAnalyzer::Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event){
    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::IsolatedTrack>> isoTracks;

    if(!isNANO){
        event->getByToken(isoTrackToken, isoTracks);
    }

    const short& isoTrackSize = isNANO ? 0 : isoTracks->size();
    nTracks = 0;

    //Loop over all electrons
    for(int i = 0; i < isoTrackSize; i++){
        const float& pt = isNANO ? 0 : isoTracks->at(i).pt();
        const float& eta = isNANO ? 0 : isoTracks->at(i).eta();
        const float& phi = isNANO ? 0 : isoTracks->at(i).phi();

        if(abs(eta) < etaCut and pt > 20 and abs(isoTracks->at(i).dz()) < 0.1 and isoTracks->at(i).fromPV() > 1){
            float iso = (isoTracks->at(i).pfIsolationDR03().chargedHadronIso() + std::max({isoTracks->at(i).pfIsolationDR03().neutralHadronIso() + isoTracks->at(i).pfIsolationDR03().photonIso() - isoTracks->at(i).pfIsolationDR03().puChargedHadronIso()/2., 0.}))/pt;

            Pt[nTracks] = pt;
            Eta[nTracks] = eta;
            Phi[nTracks] = phi;
            Isolation[nTracks] = iso;
            Dxy[nTracks] = isoTracks->at(i).dxy();
            Dz[nTracks] = isoTracks->at(i).dz();

            //Charge
            Charge[nTracks] = isNANO ? 0 : isoTracks->at(i).charge();
            PDG[nTracks] = isNANO ? 0 : isoTracks->at(i).pdgId();
            FromPV[nTracks] = isNANO ? 0 : isoTracks->at(i).fromPV();

            ++nTracks;
        }
    }
}

void MiscAnalyzer::EndJob(TFile* file){
}
