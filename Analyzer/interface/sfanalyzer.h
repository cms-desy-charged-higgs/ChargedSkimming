#ifndef SFANALYZER_H
#define SFANALYZER_H

#include <memory>

#include <TH2F.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/btagcsvreader.h>

#include <correction.h>

template <typename T>
class SFAnalyzer : public BaseAnalyzer<T> {
    private:
        //Era information
        std::string era, eleEraAlias, muEraAlias;
        std::string muTriggName;
        std::string run;

        bool isData;

        //Hist of MC btag efficiency
        std::shared_ptr<TH2F> bTagEffBLooseDeepJet, bTagEffBMediumDeepJet, bTagEffBTightDeepJet, 
                              bTagEffCLooseDeepJet, bTagEffCMediumDeepJet, bTagEffCTightDeepJet, 
                              bTagEffLightLooseDeepJet, bTagEffLightMediumDeepJet, bTagEffLightTightDeepJet,
                              bTagEffBLooseDeepCSV, bTagEffBMediumDeepCSV, bTagEffBTightDeepCSV, 
                              bTagEffCLooseDeepCSV, bTagEffCMediumDeepCSV, bTagEffCTightDeepCSV, 
                              bTagEffLightLooseDeepCSV, bTagEffLightMediumDeepCSV, bTagEffLightTightDeepCSV,
                              bTotal, cTotal, lightTotal;
                              
        std::unique_ptr<correction::CorrectionSet> eleSF, muonSF, bTagSF;
        std::vector<std::string> bTagSyst, bTagSystLight;

    public:
        SFAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            TH1::AddDirectory(false);

            //Read information needed
            run = skim.get<std::string>("run");
            era = skim.get<std::string>("era");
            isData = run != "MC";

            eleSF = correction::CorrectionSet::from_file(this->CMSSWPath + sf.get<std::string>("Electron." + era + ".file"));
            muonSF = correction::CorrectionSet::from_file(this->CMSSWPath + sf.get<std::string>("Muon.SF." + era + ".file"));
            bTagSF = correction::CorrectionSet::from_file(this->CMSSWPath + sf.get<std::string>("Jet.BTag." + era));

            eleEraAlias = sf.get<std::string>("Electron." + era + ".eraAlias");
            muEraAlias = sf.get<std::string>("Muon.SF." + era + ".eraAlias");
            muTriggName = sf.get<std::string>("Muon.SF." + era + ".triggerName");

            bTagSyst = Util::GetVector<std::string>(skim, "Analyzer.Jet.BTagSyst");
            bTagSystLight = Util::GetVector<std::string>(skim, "Analyzer.Jet.BTagSystLight");

            //Set histograms
            float ptCut = skim.get<float>("Analyzer.Jet.pt." + era);
            float etaCut = skim.get<float>("Analyzer.Jet.eta." + era);

            std::vector<double> etaBins = {-etaCut, -1.4, 1.4, etaCut};
            std::vector<double> ptBins;

            ptBins = {ptCut, 50, 70, 90, 200};

            bTotal = std::make_shared<TH2F>("nTrueB", "TotalB", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            cTotal = std::make_shared<TH2F>("nTrueC", "TotalC", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            lightTotal = std::make_shared<TH2F>("nTrueLight", "TotalLight", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

            bTagEffBLooseDeepJet = std::make_shared<TH2F>("nLooseBbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffBMediumDeepJet = std::make_shared<TH2F>("nMediumBbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffBTightDeepJet = std::make_shared<TH2F>("nTightBbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

            bTagEffCLooseDeepJet = std::make_shared<TH2F>("nLooseCbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffCMediumDeepJet = std::make_shared<TH2F>("nMediumCbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffCTightDeepJet = std::make_shared<TH2F>("nTightCbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

            bTagEffLightLooseDeepJet = std::make_shared<TH2F>("nLooseLightbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffLightMediumDeepJet = std::make_shared<TH2F>("nMediumLightbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffLightTightDeepJet = std::make_shared<TH2F>("nTightLightbTagDeepJet", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

            bTagEffBLooseDeepCSV = std::make_shared<TH2F>("nLooseBbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffBMediumDeepCSV = std::make_shared<TH2F>("nMediumBbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffBTightDeepCSV  = std::make_shared<TH2F>("nTightBbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

            bTagEffCLooseDeepCSV = std::make_shared<TH2F>("nLooseCbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffCMediumDeepCSV = std::make_shared<TH2F>("nMediumCbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffCTightDeepCSV = std::make_shared<TH2F>("nTightCbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());

            bTagEffLightLooseDeepCSV = std::make_shared<TH2F>("nLooseLightbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffLightMediumDeepCSV = std::make_shared<TH2F>("nMediumLightbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
            bTagEffLightTightDeepCSV = std::make_shared<TH2F>("nTightLightbTagDeepCSV", "", ptBins.size() - 1, ptBins.data(), etaBins.size() - 1, etaBins.data());
        }

        void Analyze(T& input, Output& out){
            if(isData) return;
            
            float elePt, muEta, muPt;

            //Loop over selected electrons
            for(int i = 0; i < out.nElectrons; ++i){
                elePt = out.elePt[i] >= 30 ? out.elePt[i] : 30;
            
                out.eleRecoSF[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sf", "RecoAbove20", out.eleEta[i], elePt});
                out.eleLooseSF[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sf", "Loose", out.eleEta[i], elePt});
                out.eleMediumSF[i] =  eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sf", "Medium", out.eleEta[i], elePt});
                out.eleTightSF[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sf", "Tight", out.eleEta[i], elePt});
                out.eleMediumMVASF[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sf", "wp90iso", out.eleEta[i], elePt});
                out.eleTightMVASF[i] =  eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sf", "wp80iso", out.eleEta[i], elePt});

                out.eleRecoSFUp[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "RecoAbove20", out.eleEta[i], elePt});
                out.eleRecoSFDown[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfdown", "RecoAbove20", out.eleEta[i], elePt});

                out.eleLooseSFUp[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "Loose", out.eleEta[i], elePt});
                out.eleLooseSFDown[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfdown", "Loose", out.eleEta[i], elePt});

                out.eleMediumSFUp[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "Medium", out.eleEta[i], elePt});
                out.eleMediumSFDown[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfdown", "Medium", out.eleEta[i], elePt});

                out.eleTightSFUp[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "Tight", out.eleEta[i], elePt});
                out.eleTightSFDown[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfdown", "Tight", out.eleEta[i], elePt});

                out.eleMediumMVASFUp[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "wp90iso", out.eleEta[i], elePt});
                out.eleMediumMVASFDown[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfdown", "wp90iso", out.eleEta[i], elePt});

                out.eleTightMVASFUp[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "wp80iso", out.eleEta[i], elePt});
                out.eleTightMVASFDown[i] = eleSF->at("UL-Electron-ID-SF")->evaluate({eleEraAlias, "sfup", "wp80iso", out.eleEta[i], elePt});
            }

            //Loop over selected muons
            for(int i = 0; i < out.nMuons; ++i){
                muPt = out.muPt[i] >= 30 ? out.muPt[i] : 30;
                muEta = std::abs(out.muEta[i]);
            
                //Scale factors
                out.muLooseIsoSF[i] = muonSF->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muEraAlias, muEta, muPt, "sf"});
                out.muTightIsoSF[i] = muonSF->at("NUM_TightRelIso_DEN_TightIDandIPCut")->evaluate({muEraAlias, muEta, muPt, "sf"});
                out.muLooseSF[i] = muonSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "sf"});
                out.muMediumSF[i] = muonSF->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "sf"});
                out.muTightSF[i] = muonSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "sf"});
                out.muTriggerSF[i] = muonSF->at(muTriggName)->evaluate({muEraAlias, muEta, muPt, "sf"});

                out.muLooseIsoSFUp[i] = muonSF->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muEraAlias, muEta, muPt, "systup"});
                out.muLooseIsoSFDown[i] = muonSF->at("NUM_LooseRelIso_DEN_LooseID")->evaluate({muEraAlias, muEta, muPt, "systdown"});
                out.muTightIsoSFUp[i] = muonSF->at("NUM_TightRelIso_DEN_TightIDandIPCut")->evaluate({muEraAlias, muEta, muPt, "systup"});
                out.muTightIsoSFDown[i] = muonSF->at("NUM_TightRelIso_DEN_TightIDandIPCut")->evaluate({muEraAlias, muEta, muPt, "systdown"});

                out.muLooseSFUp[i] = muonSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "systup"});
                out.muLooseSFDown[i] = muonSF->at("NUM_LooseID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "systdown"});
                out.muMediumSFUp[i] = muonSF->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "systup"});
                out.muMediumSFDown[i] = muonSF->at("NUM_MediumID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "systdown"});
                out.muTightSFUp[i] = muonSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "systup"});
                out.muTightSFDown[i] = muonSF->at("NUM_TightID_DEN_TrackerMuons")->evaluate({muEraAlias, muEta, muPt, "systdown"});

                out.muTriggerSFUp[i] = muonSF->at(muTriggName)->evaluate({muEraAlias, muEta, muPt, "systup"});
                out.muTriggerSFDown[i] = muonSF->at(muTriggName)->evaluate({muEraAlias, muEta, muPt, "systdown"});
            }

            //Btag efficiency and SF
            for(int i = 0; i < out.nJets; ++i){
                if(std::abs(out.jetPartFlav[i]) == 5){
                    bTotal->Fill(out.jetPt[i], out.jetEta[i]);
                
                    if(out.jetDeepJetID[i] >= 3){
                        bTagEffBLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffBMediumDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffBTightDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepJetID[i] >= 2){
                        bTagEffBLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffBMediumDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepJetID[i] >= 1){
                        bTagEffBLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    } 

                    if(out.jetDeepCSVID[i] >= 3){
                        bTagEffBLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffBMediumDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffBTightDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepCSVID[i] >= 2){
                        bTagEffBLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffBMediumDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepCSVID[i] >= 1){
                        bTagEffBLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    } 
                }

                else if(std::abs(out.jetPartFlav[i]) == 4){
                    cTotal->Fill(out.jetPt[i], out.jetEta[i]);
                    
                    if(out.jetDeepJetID[i] >= 3){
                        bTagEffCLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffCMediumDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffCTightDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepJetID[i] >= 2){
                        bTagEffCLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffCMediumDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepJetID[i] >= 1){
                        bTagEffCLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    } 

                    if(out.jetDeepCSVID[i] >= 3){
                        bTagEffCLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffCMediumDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffCTightDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepCSVID[i] >= 2){
                        bTagEffCLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffCMediumDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepCSVID[i] >= 1){
                        bTagEffCLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    } 
                }

                else{
                    lightTotal->Fill(out.jetPt[i], out.jetEta[i]);
                    
                    if(out.jetDeepJetID[i] >= 3){
                        bTagEffLightLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffLightMediumDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffLightTightDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepJetID[i] >= 2){
                        bTagEffLightLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffLightMediumDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepJetID[i] >= 1){
                        bTagEffLightLooseDeepJet->Fill(out.jetPt[i], out.jetEta[i]);
                    } 

                    if(out.jetDeepCSVID[i] >= 3){
                        bTagEffLightLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffLightMediumDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffLightTightDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepCSVID[i] >= 2){
                        bTagEffLightLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                        bTagEffLightMediumDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    }

                    else if(out.jetDeepCSVID[i] >= 1){
                        bTagEffLightLooseDeepCSV->Fill(out.jetPt[i], out.jetEta[i]);
                    } 
                }
  
                //btag SF
                int flav = std::abs(out.jetPartFlav[i]) == 4 or std::abs(out.jetPartFlav[i]) == 5 ? std::abs(out.jetPartFlav[i]) : 0;
                std::string postFix = flav != 0 ? "_comb" : "_incl";        
                
                out.jetLooseDeepCSVSF[i] = bTagSF->at("deepCSV" + postFix)->evaluate({"central", "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                out.jetMediumDeepCSVSF[i] = bTagSF->at("deepCSV" + postFix)->evaluate({"central", "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                out.jetTightDeepCSVSF[i] = bTagSF->at("deepCSV" + postFix)->evaluate({"central", "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});

                out.jetLooseDeepJetSF[i] = bTagSF->at("deepJet" + postFix)->evaluate({"central", "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                out.jetMediumDeepJetSF[i] = bTagSF->at("deepJet" + postFix)->evaluate({"central", "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                out.jetTightDeepJetSF[i] = bTagSF->at("deepJet" + postFix)->evaluate({"central", "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});

                for(int bSyst = 0; bSyst < bTagSyst.size(); ++bSyst){
                    std::string shiftUp = flav != 0 ? "up_" + bTagSyst[bSyst] : "central",
                                shiftDown = flav != 0 ? "down_" + bTagSyst[bSyst] : "central";

                    out.jetLooseDeepCSVSFUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetLooseDeepCSVSFDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepCSVSFUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepCSVSFDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepCSVSFUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepCSVSFDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});

                    out.jetLooseDeepJetSFUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetLooseDeepJetSFDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepJetSFUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepJetSFDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepJetSFUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepJetSFDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                }

                for(int bSyst = 0; bSyst < bTagSystLight.size(); ++bSyst){
                    std::string shiftUp = flav == 0 ? "up_" + bTagSyst[bSyst] : "central",
                                shiftDown = flav == 0 ? "down_" + bTagSyst[bSyst] : "central";

                    out.jetLooseDeepCSVSFLightUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetLooseDeepCSVSFLightDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepCSVSFLightUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepCSVSFLightDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepCSVSFLightUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepCSVSFLightDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});

                    out.jetLooseDeepJetSFLightUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetLooseDeepJetSFLightDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepJetSFLightUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetMediumDeepJetSFLightDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepJetSFLightUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                    out.jetTightDeepJetSFLightDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.jetEta[i]), out.jetPt[i]});
                }
            }

            for(int i = 0; i < out.nSubJets; ++i){
                if(std::abs(out.subJetPartFlav[i]) == 5){
                    bTotal->Fill(out.subJetPt[i], out.subJetEta[i]);
                
                    if(out.subJetDeepJetID[i] >= 3){
                        bTagEffBLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffBMediumDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffBTightDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepJetID[i] >= 2){
                        bTagEffBLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffBMediumDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepJetID[i] >= 1){
                        bTagEffBLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    } 

                    if(out.subJetDeepCSVID[i] >= 3){
                        bTagEffBLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffBMediumDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffBTightDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepCSVID[i] >= 2){
                        bTagEffBLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffBMediumDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepCSVID[i] >= 1){
                        bTagEffBLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    } 
                }

                else if(std::abs(out.subJetPartFlav[i]) == 4){
                    cTotal->Fill(out.subJetPt[i], out.subJetEta[i]);
                    
                    if(out.subJetDeepJetID[i] >= 3){
                        bTagEffCLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffCMediumDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffCTightDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepJetID[i] >= 2){
                        bTagEffCLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffCMediumDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepJetID[i] >= 1){
                        bTagEffCLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    } 

                    if(out.subJetDeepCSVID[i] >= 3){
                        bTagEffCLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffCMediumDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffCTightDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepCSVID[i] >= 2){
                        bTagEffCLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffCMediumDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepCSVID[i] >= 1){
                        bTagEffCLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    } 
                }

                else{
                    lightTotal->Fill(out.subJetPt[i], out.subJetEta[i]);
                    
                    if(out.subJetDeepJetID[i] >= 3){
                        bTagEffLightLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffLightMediumDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffLightTightDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepJetID[i] >= 2){
                        bTagEffLightLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffLightMediumDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepJetID[i] >= 1){
                        bTagEffLightLooseDeepJet->Fill(out.subJetPt[i], out.subJetEta[i]);
                    } 

                    if(out.subJetDeepCSVID[i] >= 3){
                        bTagEffLightLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffLightMediumDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffLightTightDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepCSVID[i] >= 2){
                        bTagEffLightLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                        bTagEffLightMediumDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    }

                    else if(out.subJetDeepCSVID[i] >= 1){
                        bTagEffLightLooseDeepCSV->Fill(out.subJetPt[i], out.subJetEta[i]);
                    } 
                }
  
                //btag SF
                int flav = std::abs(out.subJetPartFlav[i]) == 4 or std::abs(out.subJetPartFlav[i]) == 5 ? std::abs(out.subJetPartFlav[i]) : 0;
                std::string postFix = flav != 0 ? "_comb" : "_incl";        
                
                out.subJetLooseDeepCSVSF[i] = bTagSF->at("deepCSV" + postFix)->evaluate({"central", "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                out.subJetMediumDeepCSVSF[i] = bTagSF->at("deepCSV" + postFix)->evaluate({"central", "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                out.subJetTightDeepCSVSF[i] = bTagSF->at("deepCSV" + postFix)->evaluate({"central", "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});

                out.subJetLooseDeepJetSF[i] = bTagSF->at("deepJet" + postFix)->evaluate({"central", "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                out.subJetMediumDeepJetSF[i] = bTagSF->at("deepJet" + postFix)->evaluate({"central", "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                out.subJetTightDeepJetSF[i] = bTagSF->at("deepJet" + postFix)->evaluate({"central", "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});

                for(int bSyst = 0; bSyst < bTagSyst.size(); ++bSyst){
                    std::string shiftUp = flav != 0 ? "up_" + bTagSyst[bSyst] : "central",
                                shiftDown = flav != 0 ? "down_" + bTagSyst[bSyst] : "central";

                    out.subJetLooseDeepCSVSFUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetLooseDeepCSVSFDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepCSVSFUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepCSVSFDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepCSVSFUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepCSVSFDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});

                    out.subJetLooseDeepJetSFUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetLooseDeepJetSFDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepJetSFUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepJetSFDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepJetSFUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepJetSFDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                }

                for(int bSyst = 0; bSyst < bTagSystLight.size(); ++bSyst){
                    std::string shiftUp = flav == 0 ? "up_" + bTagSyst[bSyst] : "central",
                                shiftDown = flav == 0 ? "down_" + bTagSyst[bSyst] : "central";

                    out.subJetLooseDeepCSVSFLightUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetLooseDeepCSVSFLightDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepCSVSFLightUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepCSVSFLightDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepCSVSFLightUp[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepCSVSFLightDown[bSyst][i] = bTagSF->at("deepCSV" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});

                    out.subJetLooseDeepJetSFLightUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetLooseDeepJetSFLightDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "L", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepJetSFLightUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetMediumDeepJetSFLightDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "M", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepJetSFLightUp[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftUp, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                    out.subJetTightDeepJetSFLightDown[bSyst][i] = bTagSF->at("deepJet" + postFix)->evaluate({shiftDown, "T", flav, std::abs(out.subJetEta[i]), out.subJetPt[i]});
                }
            }
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){
            if(isData) return;
            outFile->cd();

            bTotal->Write();
            cTotal->Write();
            lightTotal->Write();

            bTagEffBLooseDeepJet->Write();
            bTagEffBMediumDeepJet->Write();
            bTagEffBTightDeepJet->Write();

            bTagEffBLooseDeepCSV->Write();
            bTagEffBMediumDeepCSV->Write();
            bTagEffBTightDeepCSV->Write();

            bTagEffCLooseDeepJet->Write();
            bTagEffCMediumDeepJet->Write();
            bTagEffCTightDeepJet->Write();

            bTagEffCLooseDeepCSV->Write();
            bTagEffCMediumDeepCSV->Write();
            bTagEffCTightDeepCSV->Write();

            bTagEffLightLooseDeepJet->Write();
            bTagEffLightMediumDeepJet->Write();
            bTagEffLightTightDeepJet->Write();

            bTagEffLightLooseDeepCSV->Write();
            bTagEffLightMediumDeepCSV->Write();
            bTagEffLightTightDeepCSV->Write();
        }
};

#endif
