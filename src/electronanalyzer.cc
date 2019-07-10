#include <ChargedHiggs/Skimming/interface/electronanalyzer.h>

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, eToken& eleToken):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNEle(minNEle),
    eleToken(eleToken)
    {}

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNEle(minNEle)
    {}

void ElectronAnalyzer::SetGenParticles(Electron &validElectron, const int &i){
    //Check if gen matched particle exist
    if(eleGenIdx->At(i) != -1){
        validElectron.isgenMatched = true;
        int idxMotherEle = genMotherIdx->At(eleGenIdx->At(i));

        while(abs(genID->At(idxMotherEle)) == 11){
            idxMotherEle = genMotherIdx->At(idxMotherEle);
        }

        validElectron.genVec.SetPtEtaPhiM(genPt->At(eleGenIdx->At(i)), genEta->At(eleGenIdx->At(i)), genPhi->At(eleGenIdx->At(i)), genMass->At(eleGenIdx->At(i)));

        if(abs(genID->At(idxMotherEle)) == 24){
            float idxMotherW = genMotherIdx->At(idxMotherEle);

            while(abs(genID->At(idxMotherW)) == 24){
                idxMotherW = genMotherIdx->At(idxMotherW);
            }

            if(abs(genID->At(idxMotherW)) == 37){
                validElectron.isFromHc = true;
            }
        }
    }
}

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

        if(!this->isData){
            eleGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_genPartIdx");
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

    if(!isNANO){
        event->getByToken(eleToken, electrons);
    }

    float eleSize = isNANO ? elePt->GetSize() : electrons->size();

    //Loop over all electrons
    for(unsigned int i = 0; i < eleSize; i++){
        float pt = isNANO ? elePt->At(i) : electrons->at(i).pt();
        float eta = isNANO ? eleEta->At(i) : electrons->at(i).eta();
        float phi = isNANO ? elePhi->At(i) : electrons->at(i).phi();

        if(pt > ptCut && abs(eta) < etaCut){
            Electron validElectron;

            //Set electron information
            validElectron.fourVec.SetPtEtaPhiM(pt, eta, phi, 0.510*1e-3);
            //validElectron.isolation = isNANO ? eleIso->At(i) : electrons->at(i).miniPFIsolation();
            validElectron.isMedium = isNANO ? eleMediumMVA->At(i) : electrons->at(i).electronID("mvaEleID-Fall17-iso-V2-wp80");
            validElectron.isTight = isNANO ? eleMediumMVA->At(i) : electrons->at(i).electronID("mvaEleID-Fall17-iso-V2-wp90");
            validElectron.charge = isNANO ? eleCharge->At(i) : electrons->at(i).charge();
           // validElectron.isTriggerMatched = triggerMatching(validElectron.fourVec, 11);

            if(!isData){
               //Fill scale factors
               validElectron.recoSF = recoSFhist->GetBinContent(recoSFhist->FindBin(eta, pt));
               validElectron.mediumMvaSF = mediumSFhist->GetBinContent(mediumSFhist->FindBin(eta, pt));
               validElectron.tightMvaSF = tightSFhist->GetBinContent(tightSFhist->FindBin(eta, pt));

               //Save gen particle information
               //SetGenParticles(validElectron, i);
                    
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
