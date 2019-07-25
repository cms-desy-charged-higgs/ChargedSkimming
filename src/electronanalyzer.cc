#include <ChargedHiggs/Skimming/interface/electronanalyzer.h>

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, eToken& eleToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNEle(minNEle),
    eleToken(eleToken),
    triggerObjToken(triggerObjToken),
    genParticleToken(genParticleToken)
    {}

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNEle(minNEle)
    {}

void ElectronAnalyzer::BeginJob(TTree* tree, bool &isData){
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

    //Hist with scale factors
    TFile* recoSFfile = TFile::Open(recoSFfiles[era].c_str());
    recoSFhist = (TH2F*)recoSFfile->Get("EGamma_SF2D");

    TFile* mediumSFfile = TFile::Open(mediumSFfiles[era].c_str());
    mediumSFhist = (TH2F*)mediumSFfile->Get("EGamma_SF2D");

    TFile* tightSFfile = TFile::Open(tightSFfiles[era].c_str());
    tightSFhist = (TH2F*)tightSFfile->Get("EGamma_SF2D");

    //Initiliaze TTreeReaderValues then using NANO AOD
    if(isNANO){
        TTree* eventTree = reader->GetTree();

        elePt = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pt");
        eleEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_eta");
        elePhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_phi");
        eleCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_charge");
        eleIso = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pfRelIso03_all");

        if(eventTree->GetBranchStatus("Electron_mvaFall17V1Iso_WP80")){
            eleMediumMVA = std::make_unique<TTreeReaderArray<bool>>(*reader, "Electron_mvaFall17V1Iso_WP80");
            eleTightMVA = std::make_unique<TTreeReaderArray<bool>>(*reader, "Electron_mvaFall17V1Iso_WP80");
        }

        else{
            eleMediumMVA = std::make_unique<TTreeReaderArray<bool>>(*reader, "Electron_mvaFall17Iso_WP80");
            eleTightMVA = std::make_unique<TTreeReaderArray<bool>>(*reader, "Electron_mvaFall17Iso_WP80");
        }

        //Set TTreeReader for genpart and trigger obj from baseanalyzer    
        SetCollection(this->isData);
    }

    //Set Branches of output tree
    tree->Branch("electron", &validElectrons);
}

bool ElectronAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow, const edm::Event* event){
    //Clear electron vector
    validElectrons.clear();

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
        float pt = isNANO ? elePt->At(i) : electrons->at(i).pt();
        float eta = isNANO ? eleEta->At(i) : electrons->at(i).eta();
        float phi = isNANO ? elePhi->At(i) : electrons->at(i).phi();

        if(pt > ptCut && abs(eta) < etaCut){
            Electron validElectron;

            float iso;

            if(!isNANO){
                iso = (electrons->at(i).pfIsolationVariables().sumChargedHadronPt + std::max(electrons->at(i).pfIsolationVariables().sumNeutralHadronEt + electrons->at(i).pfIsolationVariables().sumPhotonEt - 0.5 * electrons->at(i).pfIsolationVariables().sumPUPt, 0.0)) / electrons->at(i).pt();
            }

            //Set electron information
            validElectron.fourVec.SetPtEtaPhiM(pt, eta, phi, 0.510*1e-3);
            validElectron.isolation = isNANO ? eleIso->At(i) : iso;
            validElectron.isMedium = isNANO ? eleMediumMVA->At(i) : electrons->at(i).electronID("mvaEleID-Fall17-iso-V2-wp80");
            validElectron.isTight = isNANO ? eleMediumMVA->At(i) : electrons->at(i).electronID("mvaEleID-Fall17-iso-V2-wp90");
            validElectron.charge = isNANO ? eleCharge->At(i) : electrons->at(i).charge();
            validElectron.isTriggerMatched = triggerMatching(validElectron.fourVec, *trigObjects);

            if(!isData){
               //Fill scale factors
                validElectron.recoSF = recoSFhist->GetBinContent(recoSFhist->FindBin(eta, pt));
                validElectron.mediumMvaSF = mediumSFhist->GetBinContent(mediumSFhist->FindBin(eta, pt));
                validElectron.tightMvaSF = tightSFhist->GetBinContent(tightSFhist->FindBin(eta, pt));

                //Save gen particle information
                if(isNANO) SetGenParticles<Electron>(validElectron, i, 11);
                else{
                    event->getByToken(genParticleToken, genParts);
                    SetGenParticles<Electron>(validElectron, i, 11, *genParts);
                }
            }

            //Fill electron in collection
            validElectrons.push_back(validElectron);
        } 

    }
    
    //Check if event has enough electrons
    if(validElectrons.size() < minNEle){
        return false;
    }

    if(minNEle != 0){
        std::string cutName("N_{e} >= " + std::to_string(minNEle) + " (no iso/ID req)");
        cutflow.first->Fill(cutName.c_str(), cutflow.second);
    }

    return true;
}


void ElectronAnalyzer::EndJob(TFile* file){
}
