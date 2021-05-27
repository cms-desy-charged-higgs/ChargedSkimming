#include <ChargedSkimming/Skimming/interface/electronanalyzer.h>

ElectronAnalyzer::ElectronAnalyzer(const std::string& era, const std::shared_ptr<Token>& tokens, std::string& systematic):
    BaseAnalyzer(),    
    era(era),
    tokens(tokens)
    {
        std::vector<std::string> validSystematics = {"energyScaleUp", "energyScaleDown", "energySigmaUp", "energySigmaDown"};
    
        if(std::find(validSystematics.begin(), validSystematics.end(), systematic) != validSystematics.end()){
            energyCorrection = systematic;
        }
        else energyCorrection = "ecalTrkEnergyPostCorr";
    }

ElectronAnalyzer::ElectronAnalyzer(const std::string& era, TTreeReader& reader):
    BaseAnalyzer(&reader),    
    era(era)
    {}

void ElectronAnalyzer::BeginJob(std::vector<TTree*>& trees, pt::ptree& skim, pt::ptree& sf){
    //Set data bool
    this->isData = skim.get<bool>("isData");
    this->isSyst = skim.get<bool>("isSyst");

    TH1::AddDirectory(false);

    //Kinematic cut
    ptCut = skim.get<float>("Analyzer.Electron.pt." + era);
    etaCut = skim.get<float>("Analyzer.Electron.eta." + era);

    //Hist with scale factors
    std::shared_ptr<TFile> recoSFfile = std::make_shared<TFile>((filePath + sf.get<std::string>("Electron.Reco." + era)).c_str());
    recoSFhist.reset((TH2F*)recoSFfile->Get("EGamma_SF2D")); 

    std::shared_ptr<TFile> looseSFfile = std::make_shared<TFile>((filePath + sf.get<std::string>("Electron.ID.Loose." + era)).c_str());
    looseSFhist.reset((TH2F*)looseSFfile->Get("EGamma_SF2D"));

    std::shared_ptr<TFile> mediumSFfile = std::make_shared<TFile>((filePath + sf.get<std::string>("Electron.ID.Medium." + era)).c_str());
    mediumSFhist.reset((TH2F*)mediumSFfile->Get("EGamma_SF2D"));

    std::shared_ptr<TFile> tightSFfile = std::make_shared<TFile>((filePath + sf.get<std::string>("Electron.ID.Tight." + era)).c_str());
    tightSFhist.reset((TH2F*)tightSFfile->Get("EGamma_SF2D"));

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
        std::map<std::string, float*> SFvariations = {
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
        tree->Branch("Electron_Size", &nElectrons, "Electron_Size/S");

        for(std::pair<const std::string, float*>& var: floatVar){
            std::fill_n(var.second, 20, 1.);
            tree->Branch(("Electron_" + var.first).c_str(), var.second, ("Electron_" + var.first + "[Electron_Size]/F").c_str());
        }

        for(std::pair<const std::string, short*>& var: intVar){
            std::fill_n(var.second, 20, 0);
            tree->Branch(("Electron_" + var.first).c_str(), var.second, ("Electron_" + var.first + "[Electron_Size]/S").c_str());
        }
    }
}

void ElectronAnalyzer::Analyze(std::vector<CutFlow>& cutflows, const edm::Event* event){
    //Get Event info is using MINIAOD
    std::vector<pat::Electron> electrons;
    std::vector<reco::GenParticle> genParts;

    if(!isNANO){
        electrons = Token::GetTokenValue<std::vector<pat::Electron>>(event, tokens->eleToken);
    }

    const int& eleSize = isNANO ? elePt->GetSize() : electrons.size();
    nElectrons = 0;
        
    //Loop over all electrons
    for(int i = 0; i < eleSize; ++i){
        if(nElectrons >= std::size(Pt)) break;

        const float& pt = isNANO ? elePt->At(i) : (electrons.at(i).p4()*electrons.at(i).userFloat(energyCorrection) / electrons.at(i).energy()).Pt();
        const float& eta = isNANO ? eleEta->At(i) : (electrons.at(i).p4()*electrons.at(i).userFloat(energyCorrection) / electrons.at(i).energy()).Eta();
        const float& phi = isNANO ? elePhi->At(i) : (electrons.at(i).p4()*electrons.at(i).userFloat(energyCorrection) / electrons.at(i).energy()).Phi();

        if(pt > ptCut && abs(eta) < etaCut){
            //Electron four momentum components
            Pt[nElectrons] = pt;
            Eta[nElectrons] =  eta;
            Phi[nElectrons] =  phi;

            Charge[nElectrons] = isNANO ? eleCharge->At(i) : electrons.at(i).charge();  //charge

            //Isolation
            const float& iso = isNANO ? eleIso->At(i) : (electrons.at(i).pfIsolationVariables().sumChargedHadronPt + std::max(electrons.at(i).pfIsolationVariables().sumNeutralHadronEt + electrons.at(i).pfIsolationVariables().sumPhotonEt - 0.5 * electrons.at(i).pfIsolationVariables().sumPUPt, 0.0)) / electrons.at(i).pt();
   
            Isolation[nElectrons] =  iso;

            //Electron ID
            if(isNANO ? eleID->At(i) == 4 : electrons.at(i).electronID("cutBasedElectronID-Fall17-94X-V2-tight")) ID[nElectrons] = 3;
            else if(isNANO ? eleID->At(i) == 3 : electrons.at(i).electronID("cutBasedElectronID-Fall17-94X-V2-medium")) ID[nElectrons] = 2;
            else if(isNANO ? eleID->At(i) == 2 : electrons.at(i).electronID("cutBasedElectronID-Fall17-94X-V2-loose")) ID[nElectrons] = 1;
            else ID[nElectrons] = 0;

            if(!isData){
               //Fill scale factors
                const Int_t& recoBin = recoSFhist->FindBin(eta, pt);
                const Int_t& looseBin = looseSFhist->FindBin(eta, pt);
                const Int_t& mediumBin = mediumSFhist->FindBin(eta, pt);
                const Int_t& tightBin = tightSFhist->FindBin(eta, pt);

                recoSF[nElectrons] = recoSFhist->GetBinContent(recoBin);
                looseSF[nElectrons] = looseSFhist->GetBinContent(looseBin);
                mediumSF[nElectrons] =  mediumSFhist->GetBinContent(mediumBin);
                tightSF[nElectrons] = tightSFhist->GetBinContent(tightBin);

                if(!isSyst){
                    recoSFUp[nElectrons] = recoSFhist->GetBinContent(recoBin) + recoSFhist->GetBinErrorUp(recoBin);
                    recoSFDown[nElectrons] = recoSFhist->GetBinContent(recoBin) - recoSFhist->GetBinErrorLow(recoBin);

                    looseSFUp[nElectrons] = looseSFhist->GetBinContent(looseBin) + looseSFhist->GetBinErrorUp(looseBin);
                    looseSFDown[nElectrons] = looseSFhist->GetBinContent(looseBin) - looseSFhist->GetBinErrorLow(looseBin);

                    mediumSFUp[nElectrons] = mediumSFhist->GetBinContent(mediumBin) + mediumSFhist->GetBinErrorUp(mediumBin);
                    mediumSFDown[nElectrons] = mediumSFhist->GetBinContent(mediumBin) - mediumSFhist->GetBinErrorLow(mediumBin);

                    tightSFUp[nElectrons] = tightSFhist->GetBinContent(tightBin) + tightSFhist->GetBinErrorUp(tightBin);
                    tightSFDown[nElectrons] = tightSFhist->GetBinContent(tightBin) - tightSFhist->GetBinErrorLow(tightBin);
                }

                //Save gen particle information
                std::tuple<int, int, int> IDs;

                if(isNANO){
                    IDs = SetGenParticles(pt, eta, phi, i, 13);
                    partID[nElectrons] = std::get<0>(IDs);
                    mothID[nElectrons] = std::get<1>(IDs);
                    grandID[nElectrons] = std::get<2>(IDs);
                }
    
                else{
                    genParts = tokens->GetTokenValue(event, tokens->genPartToken);

                    IDs = SetGenParticles(pt, eta, phi, i, 13, genParts);
                    partID[nElectrons] = std::get<0>(IDs);
                    mothID[nElectrons] = std::get<1>(IDs);
                    grandID[nElectrons] = std::get<2>(IDs);
                }
            }

            ++nElectrons; 
        }
    }

    for(CutFlow &cutflow: cutflows){
        if(cutflow.nMinEle <= nElectrons){
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
