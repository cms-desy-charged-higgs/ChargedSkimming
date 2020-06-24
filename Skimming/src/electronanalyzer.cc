#include <ChargedSkimming/Skimming/interface/electronanalyzer.h>

ElectronAnalyzer::ElectronAnalyzer(const int& era, const float& ptCut, const float& etaCut, eToken& eleToken, genPartToken& genParticleToken, const std::string& systematic):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    eleToken(eleToken),
    genParticleToken(genParticleToken)
    {
        energyCorrection = systematic == "" ? "ecalTrkEnergyPostCorr" : systematic; 
    }

ElectronAnalyzer::ElectronAnalyzer(const int& era, const float& ptCut, const float& etaCut, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

void ElectronAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    //SF files
    looseSFfiles = {
                    {2017, filePath + "eleSF/2017_ElectronLoose.root"},
    };

    mediumSFfiles = {
                    {2017, filePath + "eleSF/2017_ElectronMedium.root"},
    };

    tightSFfiles = {
                    {2017, filePath + "eleSF/2017_ElectronTight.root"},
    };

    recoSFfiles = {
                    {2017, filePath + "eleSF/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"},
    };

    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    //Hist with scale factors
    TFile* recoSFfile = TFile::Open(recoSFfiles[era].c_str());
    recoSFhist = (TH2F*)recoSFfile->Get("EGamma_SF2D");

    TFile* looseSFfile = TFile::Open(looseSFfiles[era].c_str());
    looseSFhist = (TH2F*)looseSFfile->Get("EGamma_SF2D");

    TFile* mediumSFfile = TFile::Open(mediumSFfiles[era].c_str());
    mediumSFhist = (TH2F*)mediumSFfile->Get("EGamma_SF2D");

    TFile* tightSFfile = TFile::Open(tightSFfiles[era].c_str());
    tightSFhist = (TH2F*)tightSFfile->Get("EGamma_SF2D");

    //Initiliaze TTreeReaderValues then using NANO AOD
    if(isNANO){
        elePt = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pt");
        eleEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_eta");
        elePhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_phi");
        eleCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_charge");
        eleIso = std::make_unique<TTreeReaderArray<float>>(*reader, "Electron_pfRelIso03_all");
        eleID = std::make_unique<TTreeReaderArray<int>>(*reader, "Electron_cutBased");

        //Set TTreeReader for genpart and trigger obj from baseanalyzer    
        SetCollection(this->isData);
    }

    //Variable name mapping to branch name
    floatVar = {
            {"Pt", Pt},
            {"Eta", Eta},
            {"Phi", Phi},
            {"Isolation", Isolation},
            {"recoSF", recoSF},
            {"looseSF", looseSF},
            {"mediumSF", mediumSF},
            {"tightSF", tightSF},
    };

    intVar = {
            {"Charge", Charge},
            {"ID", ID},
            {"ParticleID", partID},
            {"MotherID", mothID},
            {"GrandMotherID", grandID},
    };

    if(!isSyst){
        std::map<std::string, std::vector<float>&> SFvariations = {
            {"recoSFUp", recoSFUp},
            {"recoSFDown", recoSFDown},
            {"looseSFUp", looseSFUp},
            {"looseSFDown", looseSFDown},
            {"mediumSFUp", mediumSFUp},
            {"mediumSFDown", mediumSFDown},
            {"tightSFUp", tightSFUp},
            {"tightSFDown", tightSFDown},
        };

        floatVar.insert(SFvariations.begin(), SFvariations.end());   
    }

    //Set Branches of output tree
    for(TTree* tree: trees){
        tree->Branch("Electron_Size", &nElectrons);

        for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
            tree->Branch(("Electron_" + var.first).c_str(), &var.second);
        }

        for(std::pair<const std::string, std::vector<char>&>& var: intVar){
            tree->Branch(("Electron_" + var.first).c_str(), &var.second);
        }
    }
}

void ElectronAnalyzer::Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
        var.second.clear();
    }

    for(std::pair<const std::string, std::vector<char>&>& var: intVar){
        var.second.clear();
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Electron>> electrons;
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    if(!isNANO){
        event->getByToken(eleToken, electrons);
    }

    const int& eleSize = isNANO ? elePt->GetSize() : electrons->size();
    
    //Loop over all electrons
    for(int i = 0; i < eleSize; i++){
        const float& pt = isNANO ? elePt->At(i) : (electrons->at(i).p4()*electrons->at(i).userFloat(energyCorrection) / electrons->at(i).energy()).Pt();
        const float& eta = isNANO ? eleEta->At(i) : (electrons->at(i).p4()*electrons->at(i).userFloat(energyCorrection) / electrons->at(i).energy()).Eta();
        const float& phi = isNANO ? elePhi->At(i) : (electrons->at(i).p4()*electrons->at(i).userFloat(energyCorrection) / electrons->at(i).energy()).Phi();

        if(pt > ptCut && abs(eta) < etaCut){
            //Electron four momentum components
            Pt.push_back(pt);
            Eta.push_back(eta);
            Phi.push_back(phi);

            Charge.push_back(isNANO ? eleCharge->At(i) : electrons->at(i).charge());  //charge

            //Isolation
            const float& iso = isNANO ? eleIso->At(i) : (electrons->at(i).pfIsolationVariables().sumChargedHadronPt + std::max(electrons->at(i).pfIsolationVariables().sumNeutralHadronEt + electrons->at(i).pfIsolationVariables().sumPhotonEt - 0.5 * electrons->at(i).pfIsolationVariables().sumPUPt, 0.0)) / electrons->at(i).pt();
   
            Isolation.push_back(iso);

            //Electron ID
            if(isNANO ? eleID->At(i) == 4 : electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V2-tight")) ID.push_back(3);
            else if(isNANO ? eleID->At(i) == 3 : electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V2-medium")) ID.push_back(2);
            else if(isNANO ? eleID->At(i) == 2 : electrons->at(i).electronID("cutBasedElectronID-Fall17-94X-V2-loose")) ID.push_back(1);
            else ID.push_back(0);

            if(!isData){
               //Fill scale factors
                const Int_t& recoBin = recoSFhist->FindBin(eta, pt);
                const Int_t& looseBin = looseSFhist->FindBin(eta, pt);
                const Int_t& mediumBin = mediumSFhist->FindBin(eta, pt);
                const Int_t& tightBin = tightSFhist->FindBin(eta, pt);

                recoSF.push_back(recoSFhist->GetBinContent(recoBin));
                looseSF.push_back(looseSFhist->GetBinContent(looseBin));
                mediumSF.push_back(mediumSFhist->GetBinContent(mediumBin));
                tightSF.push_back(tightSFhist->GetBinContent(tightBin));

                if(!isSyst){
                    recoSFUp.push_back(recoSFhist->GetBinContent(recoBin) + recoSFhist->GetBinErrorUp(recoBin));
                    recoSFDown.push_back(recoSFhist->GetBinContent(recoBin) - recoSFhist->GetBinErrorLow(recoBin));

                    looseSFUp.push_back(looseSFhist->GetBinContent(looseBin) + looseSFhist->GetBinErrorUp(looseBin));
                    looseSFDown.push_back(looseSFhist->GetBinContent(looseBin) - looseSFhist->GetBinErrorLow(looseBin));

                    mediumSFUp.push_back(mediumSFhist->GetBinContent(mediumBin) + mediumSFhist->GetBinErrorUp(mediumBin));
                    mediumSFDown.push_back(mediumSFhist->GetBinContent(mediumBin) - mediumSFhist->GetBinErrorLow(mediumBin));

                    tightSFUp.push_back(tightSFhist->GetBinContent(tightBin) + tightSFhist->GetBinErrorUp(tightBin));
                    tightSFDown.push_back(tightSFhist->GetBinContent(tightBin) - tightSFhist->GetBinErrorLow(tightBin));
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

    nElectrons = Pt.size();

    for(CutFlow &cutflow: cutflows){
        if(cutflow.nMinEle <= Pt.size()){
            if(cutflow.nMinEle!=0 and cutflow.passed){
                std::string cutName("N_{e} >= " + std::to_string(cutflow.nMinEle) + " (no iso/ID req)");
                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
            }
        }

        else{
            cutflow.passed = false;
        }
    }
}

void ElectronAnalyzer::EndJob(TFile* file){
}
