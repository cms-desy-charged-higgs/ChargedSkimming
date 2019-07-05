#include <ChargedHiggs/Skimming/interface/electronanalyzer.h>

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNEle(minNEle)
    {}

ElectronAnalyzer::ElectronAnalyzer(const int &era, const float &ptCut, const float &etaCut, const int &minNEle, edm::Handle<std::vector<pat::Electron>>& electrons):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    minNEle(minNEle),
    electrons(electrons)
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
    TTree* eventTree = reader->GetTree();

    //Hist with scale factors
    TFile* recoSFfile = TFile::Open(recoSFfiles[era].c_str());
    recoSFhist = (TH2F*)recoSFfile->Get("EGamma_SF2D");

    TFile* mediumSFfile = TFile::Open(mediumSFfiles[era].c_str());
    mediumSFhist = (TH2F*)mediumSFfile->Get("EGamma_SF2D");

    TFile* tightSFfile = TFile::Open(tightSFfiles[era].c_str());
    tightSFhist = (TH2F*)tightSFfile->Get("EGamma_SF2D");

    //Initiliaze TTreeReaderValues
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

    //Set Branches of output tree
    tree->Branch("electron", &validElectrons);
}

bool ElectronAnalyzer::Analyze(std::pair<TH1F*, float> &cutflow){
    //Clear electron vector
    validElectrons.clear();

    //Loop over all electrons
    for(unsigned int i = 0; i < elePt->GetSize(); i++){
        if(elePt->At(i) > ptCut && abs(eleEta->At(i)) < etaCut){
            Electron validElectron;

            //Set electron information
            validElectron.fourVec.SetPtEtaPhiM(elePt->At(i), eleEta->At(i), elePhi->At(i), 0.510*1e-3);
            validElectron.isolation = eleIso->At(i);
            validElectron.isMedium = eleMediumMVA->At(i);
            validElectron.isTight = eleTightMVA->At(i);
            validElectron.charge = eleCharge->At(i);
            validElectron.isTriggerMatched = triggerMatching(validElectron.fourVec, 11);

            if(!isData){
                    //Fill scale factors
                validElectron.recoSF = recoSFhist->GetBinContent(recoSFhist->FindBin(eleEta->At(i), elePt->At(i)));
                validElectron.mediumMvaSF = mediumSFhist->GetBinContent(mediumSFhist->FindBin(eleEta->At(i), elePt->At(i)));
                validElectron.tightMvaSF = tightSFhist->GetBinContent(tightSFhist->FindBin(eleEta->At(i), elePt->At(i)));

                //Save gen particle information
                SetGenParticles(validElectron, i);
                    
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
