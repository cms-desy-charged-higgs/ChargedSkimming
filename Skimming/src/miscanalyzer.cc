#include <ChargedSkimming/Skimming/interface/miscanalyzer.h>

MiscAnalyzer::MiscAnalyzer(const int& era, const float& etaCut, isoToken& isoTrackToken, lheToken& lheTok):
    BaseAnalyzer(), 
    era(era),
    etaCut(etaCut),
    isoTrackToken(isoTrackToken),
    lheTok(lheTok)
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
        tree->Branch("Misc_NParton", &nParton, "Misc_NParton/S");

        for(std::pair<const std::string, float*>& var: floatVar){
            tree->Branch(("IsoTrack_" + var.first).c_str(), var.second, ("IsoTrack_" + var.first + "[IsoTrack_Size]/F").c_str());
        }

        for(std::pair<const std::string, short*>& var: intVar){
            tree->Branch(("IsoTrack_" + var.first).c_str(), var.second, ("IsoTrack_" + var.first + "[IsoTrack_Size]/S").c_str());
        }
    }
}

int MiscAnalyzer::GetNParton(const edm::Event* event){
    edm::Handle<LHEEventProduct> lheEventProduct;  
    event->getByToken(lheTok, lheEventProduct);
    const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
    
    size_t nOutgoing = 0;

    for(size_t idxParticle = 0; idxParticle < lheEvent.PUP.size(); ++idxParticle){
        int absPdgId = std::abs(lheEvent.IDUP[idxParticle]);
        int status = lheEvent.ISTUP[idxParticle];
        
        if (status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21)){  // quarks and gluons
            ++nOutgoing;
        } 
    }

    return nOutgoing;
}

void MiscAnalyzer::Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event){
    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::IsolatedTrack>> isoTracks;

    if(!isNANO){
        event->getByToken(isoTrackToken, isoTracks);
    }

    const short& isoTrackSize = isNANO ? 0 : isoTracks->size();
    nTracks = 0, nParton = 0;

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

    //Get Number of true partons
    if(!isData and event){
        nParton = GetNParton(event);
    }
}

void MiscAnalyzer::EndJob(TFile* file){
}
