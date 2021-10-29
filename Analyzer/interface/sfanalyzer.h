#ifndef SFANALYZER_H
#define SFANALYZER_H

#include <memory>

#include <TH2F.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/btagcsvreader.h>

template <typename T>
class SFAnalyzer : public BaseAnalyzer<T> {
    private:
        //Era information
        std::string era;
        std::string run;

        bool isSyst, isData;

        //Hist with electron scale factors
        std::shared_ptr<TH2F> eleLooseSFhist, eleMediumSFhist, eleTightSFhist, eleMediumMVASFhist, eleTightMVASFhist, eleRecoSFhist;

        //Hist with electron scale factorsMuon
        std::shared_ptr<TH2F> muLooseSFhist, muMediumSFhist, muTightSFhist, muTriggerSFhist, muLooseIsoSFhist, muTightIsoSFhist;

        //Classes for reading btag SF
        BTagCSVReader DeepCSVReader, DeepJetReader;

        //Hist of MC btag efficiency
        std::shared_ptr<TH2F> bTagEffBLooseDeepJet, bTagEffBMediumDeepJet, bTagEffBTightDeepJet, 
                              bTagEffCLooseDeepJet, bTagEffCMediumDeepJet, bTagEffCTightDeepJet, 
                              bTagEffLightLooseDeepJet, bTagEffLightMediumDeepJet, bTagEffLightTightDeepJet,
                              bTagEffBLooseDeepCSV, bTagEffBMediumDeepCSV, bTagEffBTightDeepCSV, 
                              bTagEffCLooseDeepCSV, bTagEffCMediumDeepCSV, bTagEffCTightDeepCSV, 
                              bTagEffLightLooseDeepCSV, bTagEffLightMediumDeepCSV, bTagEffLightTightDeepCSV,
                              bTotal, cTotal, lightTotal;

    public:
        SFAnalyzer(){}

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            TH1::AddDirectory(false);

            //Read information needed
            run = skim.get<std::string>("run");
            era = skim.get<std::string>("era");
            isSyst = skim.get<bool>("isSyst");
            isData = run != "MC";

            //Electron SF hists
            std::shared_ptr<TFile> eleRecoSFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Electron.Reco." + era)).c_str());
            eleRecoSFhist.reset((TH2F*)eleRecoSFfile->Get("EGamma_SF2D")); 

            std::shared_ptr<TFile> eleLooseSFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Electron.ID.Loose." + era)).c_str());
            eleLooseSFhist.reset((TH2F*)eleLooseSFfile->Get("EGamma_SF2D"));

            std::shared_ptr<TFile> eleMediumSFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Electron.ID.Medium." + era)).c_str());
            eleMediumSFhist.reset((TH2F*)eleMediumSFfile->Get("EGamma_SF2D"));

            std::shared_ptr<TFile> eleTightSFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Electron.ID.Tight." + era)).c_str());
            eleTightSFhist.reset((TH2F*)eleTightSFfile->Get("EGamma_SF2D"));

            std::shared_ptr<TFile> eleMediumMVASFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Electron.ID.MVAMedium." + era)).c_str());
            eleMediumMVASFhist.reset((TH2F*)eleMediumMVASFfile->Get("EGamma_SF2D"));

            std::shared_ptr<TFile> eleTightMVASFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Electron.ID.MVATight." + era)).c_str());
            eleTightMVASFhist.reset((TH2F*)eleTightMVASFfile->Get("EGamma_SF2D"));

            //Muon SF hists
            std::shared_ptr<TFile> muTriggerSFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Muon.Trigger.File." + era)).c_str());
            muTriggerSFhist.reset((TH2F*)muTriggerSFfile->Get(sf.get<std::string>("Muon.Trigger.Histogram." + era).c_str()));

            std::shared_ptr<TFile> muIsoSFfile = std::make_shared<TFile>((this->filePath + sf.get<std::string>("Muon.Isolation." + era)).c_str());

            muLooseIsoSFhist.reset((TH2F*)muIsoSFfile->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_abseta_pt"));
            muTightIsoSFhist.reset((TH2F*)muIsoSFfile->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"));

            std::shared_ptr<TFile> muIDSFfile = std::make_shared<TFile>((this->filePath  + sf.get<std::string>("Muon.ID.File." + era)).c_str());
            muLooseSFhist.reset((TH2F*)muIDSFfile->Get(sf.get<std::string>("Muon.ID.Histogram.Loose." + era).c_str()));
            muMediumSFhist.reset((TH2F*)muIDSFfile->Get(sf.get<std::string>("Muon.ID.Histogram.Medium." + era).c_str()));
            muTightSFhist.reset((TH2F*)muIDSFfile->Get(sf.get<std::string>("Muon.ID.Histogram.Tight." + era).c_str()));

            //Set configuration for bTagSF reader  
            DeepCSVReader = BTagCSVReader(this->filePath + sf.get<std::string>("Jet.BTag.CSV." + era));
            DeepJetReader = BTagCSVReader(this->filePath + sf.get<std::string>("Jet.BTag.DeepJet." + era));

            //Set histograms
            float ptCut = skim.get<float>("Analyzer.Jet.pt." + era);
            float etaCut = skim.get<float>("Analyzer.Jet.eta." + era);

            std::vector<double> etaBins = {-etaCut, -1.4, 1.4, etaCut};
            std::vector<double> ptBins;

            for(std::size_t i = 0; i <= 20; ++i){
                ptBins.push_back(ptCut + i*10.);
            } 

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

            //Loop over selected electrons
            for(int i = 0; i < out.nElectrons; ++i){
                const int recoBin = eleRecoSFhist->FindBin(out.eleEta[i], out.elePt[i]);
                const int looseBin = eleLooseSFhist->FindBin(out.eleEta[i], out.elePt[i]);
                const int mediumBin = eleMediumSFhist->FindBin(out.eleEta[i], out.elePt[i]);
                const int tightBin = eleTightSFhist->FindBin(out.eleEta[i], out.elePt[i]);
                const int mediumMVABin = eleMediumMVASFhist->FindBin(out.eleEta[i], out.elePt[i]);
                const int tightMVABin = eleTightMVASFhist->FindBin(out.eleEta[i], out.elePt[i]);

                out.eleRecoSF[i] = eleRecoSFhist->GetBinContent(recoBin);
                out.eleLooseSF[i] = eleLooseSFhist->GetBinContent(looseBin);
                out.eleMediumSF[i] =  eleMediumSFhist->GetBinContent(mediumBin);
                out.eleTightSF[i] = eleTightSFhist->GetBinContent(tightBin);
                out.eleMediumMVASF[i] =  eleMediumMVASFhist->GetBinContent(mediumMVABin);
                out.eleTightMVASF[i] =  eleTightMVASFhist->GetBinContent(tightMVABin);


                if(!isSyst){
                    out.eleRecoSFUp[i] = eleRecoSFhist->GetBinContent(recoBin) + eleRecoSFhist->GetBinErrorUp(recoBin);
                    out.eleRecoSFDown[i] = eleRecoSFhist->GetBinContent(recoBin) - eleRecoSFhist->GetBinErrorLow(recoBin);

                    out.eleLooseSFUp[i] = eleLooseSFhist->GetBinContent(looseBin) + eleLooseSFhist->GetBinErrorUp(looseBin);
                    out.eleLooseSFDown[i] = eleLooseSFhist->GetBinContent(looseBin) - eleLooseSFhist->GetBinErrorLow(looseBin);

                    out.eleMediumSFUp[i] = eleMediumSFhist->GetBinContent(mediumBin) + eleMediumSFhist->GetBinErrorUp(mediumBin);
                    out.eleMediumSFDown[i] = eleMediumSFhist->GetBinContent(mediumBin) - eleMediumSFhist->GetBinErrorLow(mediumBin);

                    out.eleTightSFUp[i] = eleTightSFhist->GetBinContent(tightBin) + eleTightSFhist->GetBinErrorUp(tightBin);
                    out.eleTightSFDown[i] = eleTightSFhist->GetBinContent(tightBin) - eleTightSFhist->GetBinErrorLow(tightBin);

                    out.eleMediumMVASFUp[i] = eleMediumMVASFhist->GetBinContent(mediumMVABin) + eleMediumMVASFhist->GetBinErrorUp(mediumMVABin);
                    out.eleMediumMVASFDown[i] = eleMediumMVASFhist->GetBinContent(mediumMVABin) - eleMediumMVASFhist->GetBinErrorLow(mediumMVABin);

                    out.eleTightMVASFUp[i] = eleTightMVASFhist->GetBinContent(tightMVABin) + eleTightMVASFhist->GetBinErrorUp(tightMVABin);
                    out.eleTightMVASFDown[i] = eleTightMVASFhist->GetBinContent(tightMVABin) - eleTightMVASFhist->GetBinErrorLow(tightMVABin);
                }
            }

            //Loop over selected muons
            for(int i = 0; i < out.nMuons; ++i){
                //Scale factors
                const int looseIsoBin = muLooseIsoSFhist->FindBin(abs(out.muEta[i]), out.muPt[i]);
                const int tightIsoBin = muTightIsoSFhist->FindBin(abs(out.muEta[i]), out.muPt[i]);
                const int looseIDBin = muLooseSFhist->FindBin(abs(out.muEta[i]), out.muPt[i]);
                const int mediumIDBin = muMediumSFhist->FindBin(abs(out.muEta[i]), out.muPt[i]);
                const int tightIDBin = muTightSFhist->FindBin(abs(out.muEta[i]), out.muPt[i]);
                const int triggerBin = muTriggerSFhist->FindBin(abs(out.muEta[i]), out.muPt[i]);

                out.muLooseIsoSF[i] = muLooseIsoSFhist->GetBinContent(looseIsoBin);
                out.muTightIsoSF[i] = muTightIsoSFhist->GetBinContent(tightIsoBin);
                out.muLooseSF[i] = muLooseSFhist->GetBinContent(looseIDBin);
                out.muMediumSF[i] = muMediumSFhist->GetBinContent(mediumIDBin);
                out.muTightSF[i] = muTightSFhist->GetBinContent(triggerBin);
                out.muTriggerSF[i] = muTriggerSFhist->GetBinContent(triggerBin);

                if(!isSyst){
                    out.muLooseIsoSFUp[i] = muLooseIsoSFhist->GetBinContent(looseIsoBin) + muLooseIsoSFhist->GetBinErrorUp(looseIsoBin);
                    out.muLooseIsoSFDown[i] = muLooseIsoSFhist->GetBinContent(looseIsoBin) - muLooseIsoSFhist->GetBinContent(looseIsoBin);
                    out.muTightIsoSFUp[i] = muTightIsoSFhist->GetBinContent(tightIsoBin) + muTightIsoSFhist->GetBinErrorUp(tightIsoBin);
                    out.muTightIsoSFDown[i] = muTightIsoSFhist->GetBinContent(tightIsoBin) - muTightIsoSFhist->GetBinErrorLow(tightIsoBin);

                    out.muLooseSFUp[i] = muLooseSFhist->GetBinContent(looseIDBin) + muLooseSFhist->GetBinErrorUp(looseIDBin);
                    out.muLooseSFDown[i] = muLooseSFhist->GetBinContent(looseIDBin) - muLooseSFhist->GetBinErrorLow(looseIDBin);
                    out.muMediumSFUp[i] = muMediumSFhist->GetBinContent(mediumIDBin) + muMediumSFhist->GetBinErrorUp(mediumIDBin);
                    out.muMediumSFDown[i] = muMediumSFhist->GetBinContent(mediumIDBin) - muMediumSFhist->GetBinErrorLow(mediumIDBin);
                    out.muTightSFUp[i] = muTightSFhist->GetBinContent(tightIDBin) + muTightSFhist->GetBinErrorUp(tightIDBin);
                    out.muTightSFDown[i] = muTightSFhist->GetBinContent(tightIDBin) - muTightSFhist->GetBinErrorLow(tightIDBin);

                    out.muTriggerSFUp[i] = muTriggerSFhist->GetBinContent(triggerBin) + muTriggerSFhist->GetBinErrorUp(triggerBin);
                    out.muTriggerSFDown[i] = muTriggerSFhist->GetBinContent(triggerBin) - muTriggerSFhist->GetBinErrorLow(triggerBin);
                }
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
                out.jetLooseDeepCSVSF[i] = DeepCSVReader.Get(out.jetPt[i], 0, std::abs(out.jetPartFlav[i]));
                out.jetMediumDeepCSVSF[i] = DeepCSVReader.Get(out.jetPt[i], 1, std::abs(out.jetPartFlav[i]));
                out.jetTightDeepCSVSF[i] = DeepCSVReader.Get(out.jetPt[i], 2, std::abs(out.jetPartFlav[i]));

                out.jetLooseDeepJetSF[i] = DeepJetReader.Get(out.jetPt[i], 0, std::abs(out.jetPartFlav[i]));
                out.jetMediumDeepJetSF[i] = DeepJetReader.Get(out.jetPt[i], 1, std::abs(out.jetPartFlav[i]));
                out.jetTightDeepJetSF[i] = DeepJetReader.Get(out.jetPt[i], 2, std::abs(out.jetPartFlav[i]));

                if(!isSyst){
                    out.jetLooseDeepCSVSFUp[i] = DeepCSVReader.GetUp(out.jetPt[i], 0, std::abs(out.jetPartFlav[i]));
                    out.jetLooseDeepCSVSFDown[i] = DeepCSVReader.GetDown(out.jetPt[i], 0, std::abs(out.jetPartFlav[i]));
                    out.jetMediumDeepCSVSFUp[i] = DeepCSVReader.GetUp(out.jetPt[i], 1, std::abs(out.jetPartFlav[i]));
                    out.jetMediumDeepCSVSFDown[i] = DeepCSVReader.GetDown(out.jetPt[i], 1, std::abs(out.jetPartFlav[i]));
                    out.jetTightDeepCSVSFUp[i] = DeepCSVReader.GetUp(out.jetPt[i], 2, std::abs(out.jetPartFlav[i]));
                    out.jetTightDeepCSVSFDown[i] = DeepCSVReader.GetDown(out.jetPt[i], 2, std::abs(out.jetPartFlav[i]));

                    out.jetLooseDeepJetSFUp[i] = DeepJetReader.GetUp(out.jetPt[i], 0, std::abs(out.jetPartFlav[i]));
                    out.jetLooseDeepJetSFDown[i] = DeepJetReader.GetDown(out.jetPt[i], 0, std::abs(out.jetPartFlav[i]));
                    out.jetMediumDeepJetSFUp[i] = DeepJetReader.GetUp(out.jetPt[i], 1, std::abs(out.jetPartFlav[i]));
                    out.jetMediumDeepJetSFDown[i] = DeepJetReader.GetDown(out.jetPt[i], 1, std::abs(out.jetPartFlav[i]));
                    out.jetTightDeepJetSFUp[i] = DeepJetReader.GetUp(out.jetPt[i], 2, std::abs(out.jetPartFlav[i]));
                    out.jetTightDeepJetSFDown[i] = DeepJetReader.GetDown(out.jetPt[i], 2, std::abs(out.jetPartFlav[i]));
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
                out.subJetLooseDeepCSVSF[i] = DeepCSVReader.Get(out.subJetPt[i], 0, std::abs(out.subJetPartFlav[i]));
                out.subJetMediumDeepCSVSF[i] = DeepCSVReader.Get(out.subJetPt[i], 1, std::abs(out.subJetPartFlav[i]));
                out.subJetTightDeepCSVSF[i] = DeepCSVReader.Get(out.subJetPt[i], 2, std::abs(out.subJetPartFlav[i]));

                out.subJetLooseDeepJetSF[i] = DeepJetReader.Get(out.subJetPt[i], 0, std::abs(out.subJetPartFlav[i]));
                out.subJetMediumDeepJetSF[i] = DeepJetReader.Get(out.subJetPt[i], 1, std::abs(out.subJetPartFlav[i]));
                out.subJetTightDeepJetSF[i] = DeepJetReader.Get(out.subJetPt[i], 2, std::abs(out.subJetPartFlav[i]));

                if(!isSyst){
                    out.subJetLooseDeepCSVSFUp[i] = DeepCSVReader.GetUp(out.subJetPt[i], 0, std::abs(out.subJetPartFlav[i]));
                    out.subJetLooseDeepCSVSFDown[i] = DeepCSVReader.GetDown(out.subJetPt[i], 0, std::abs(out.subJetPartFlav[i]));
                    out.subJetMediumDeepCSVSFUp[i] = DeepCSVReader.GetUp(out.subJetPt[i], 1, std::abs(out.subJetPartFlav[i]));
                    out.subJetMediumDeepCSVSFDown[i] = DeepCSVReader.GetDown(out.subJetPt[i], 1, std::abs(out.subJetPartFlav[i]));
                    out.subJetTightDeepCSVSFUp[i] = DeepCSVReader.GetUp(out.subJetPt[i], 2, std::abs(out.subJetPartFlav[i]));
                    out.subJetTightDeepCSVSFDown[i] = DeepCSVReader.GetDown(out.subJetPt[i], 2, std::abs(out.subJetPartFlav[i]));

                    out.subJetLooseDeepJetSFUp[i] = DeepJetReader.GetUp(out.subJetPt[i], 0, std::abs(out.subJetPartFlav[i]));
                    out.subJetLooseDeepJetSFDown[i] = DeepJetReader.GetDown(out.subJetPt[i], 0, std::abs(out.subJetPartFlav[i]));
                    out.subJetMediumDeepJetSFUp[i] = DeepJetReader.GetUp(out.subJetPt[i], 1, std::abs(out.subJetPartFlav[i]));
                    out.subJetMediumDeepJetSFDown[i] = DeepJetReader.GetDown(out.subJetPt[i], 1, std::abs(out.subJetPartFlav[i]));
                    out.subJetTightDeepJetSFUp[i] = DeepJetReader.GetUp(out.subJetPt[i], 2, std::abs(out.subJetPartFlav[i]));
                    out.subJetTightDeepJetSFDown[i] = DeepJetReader.GetDown(out.subJetPt[i], 2, std::abs(out.subJetPartFlav[i]));
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
