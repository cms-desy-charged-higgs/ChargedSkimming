#ifndef JETANALYZER_H
#define JETANALYZER_H

#include <cmath>
#include <random>
#include <limits>
#include <functional>
#include <experimental/filesystem>

#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <JetMETCorrections/Modules/interface/JetCorrectionProducer.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectorParameters.h>
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>

#include <ChargedSkimming/Analyzer/interface/baseanalyzer.h>
#include <ChargedSkimming/Skimming/interface/util2.h>
#include <ChargedSkimming/Skimming/interface/btagcsvreader.h>

template <typename T>
class JetAnalyzer : public BaseAnalyzer<T> {
    private:
        //Systematic
        std::string shift;
        bool isJERSyst, isJECSyst;

        //Input information
        std::string era, run;
        bool isData, isSyst;

        //Kinematic cut criteria
        float ptCut, etaCut;

        //JEC corrector
        std::shared_ptr<FactorizedJetCorrector> jetCorrectorAK4, jetCorrectorAK8;

        //JEC uncertainty
        std::shared_ptr<JetCorrectionUncertainty> jecUncAK4, jecUncAK8;

        //Classes for reading jet energy SF 
        JME::JetParameters jetParameter;
        JME::JetResolution resolutionAK4, resolutionAK8;
        JME::JetResolutionScaleFactor resolutionSFAK4, resolutionSFAK8;

        std::default_random_engine generator;

        //BTag cuts
        std::function<bool(const float&)> isDeepCSVLoose, isDeepCSVMedium, isDeepCSVTight,
                                          isDeepJetLoose, isDeepJetMedium, isDeepJetTight;

    public:
        JetAnalyzer(const std::string& systematic, const std::string& shift) : shift(shift) {
            isJERSyst = systematic.find("JER") != std::string::npos;
            isJECSyst = systematic.find("JEC") != std::string::npos;
        }

        //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
        float CorrectEnergy(const float& pt, const float& eta, const float& phi, const float& area, const float& rho, const bool& isAK4){
            std::shared_ptr<FactorizedJetCorrector>& corr = isAK4 ? jetCorrectorAK4 : jetCorrectorAK8;

            corr->setJetPt(pt);
            corr->setJetEta(eta);
            corr->setJetPhi(phi);
            corr->setRho(rho);
            corr->setJetA(area);
            
            return corr->getCorrection();
        }

        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
        //https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
        float SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, T& input, const bool& isAK4){
            //Configure smearer class
            jetParameter.setJetPt(pt).setJetEta(eta).setRho(rho);
            JME::JetResolution& resolution = isAK4 ? resolutionAK4 : resolutionAK8;
            JME::JetResolutionScaleFactor& resolutionSF = isAK4 ? resolutionSFAK4 : resolutionSFAK8;

            float reso = resolution.getResolution(jetParameter);
            float resoSF; 
            if(!isJERSyst) resoSF = resolutionSF.getScaleFactor(jetParameter);
            else resoSF = resolutionSF.getScaleFactor(jetParameter, shift == "Up" ? Variation::UP : Variation::DOWN);

            float JME = 1., 
                  genPt = -1., 
                  dRmin = std::numeric_limits<float>::max(), 
                  coneSize = isAK4 ? 0.2 : 0.4;

            int size = isAK4 ? input.genJetSize : input.genFatJetSize;

            //Gen jet matching
            for(int i = 0; i < size; ++i){
                float dR, dPT, genJetPt;

                if(isAK4){
                    input.GetGenJet(i);
                    dR = Util2::DeltaR(eta, phi, input.genJetEta, input.genJetPhi);
                    genJetPt = input.genJetPt;
                    dPT = std::abs(pt - genJetPt);
                }

                else{
                    input.GetGenFatJet(i);
                    dR = Util2::DeltaR(eta, phi, input.genFatJetEta, input.genFatJetPhi);
                    genJetPt = input.genFatJetPt;
                    dPT = std::abs(pt - genJetPt);
                }

                if(dR > dRmin) continue;
           
                if(dR < coneSize and dPT < 3.*reso*pt){
                    dRmin = dR;
                    genPt = genJetPt;
                }
            }

            //If gen matched
            if(genPt > 0.){
                JME = 1. + (resoSF-1.)*(pt - genPt)/pt; 
            }

            else if (resoSF > 1.){
                std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));

                JME = 1. + gaus(generator);
            }

            //Check if direction of jet not changed
            if(pt*JME < 1e-2){
                JME = 1e-2/pt;
            }

            return JME;
        }

        void BeginJob(const pt::ptree& skim, const pt::ptree& sf){
            //Read information needed
            era = skim.get<std::string>("era");
            run = skim.get<std::string>("run");

            if(run != "MC" and era.find("2016") != std::string::npos){
                if(era == "2016Post") run = "FGH";
                else{
                    if(std::string("BCD").find(run) != std::string::npos) run = "BCD";
                    else run = "EF";
                }
            }

            isData = skim.get<bool>("isData");
            isSyst = skim.get<bool>("isSyst");

            ptCut = skim.get<float>("Analyzer.Jet.pt." + era);
            etaCut = skim.get<float>("Analyzer.Jet.eta." + era);

            //Set JEC/JME classes
            for(const std::string type : {"AK4", "AK8"}){
                //Class for JEC
                std::vector<JetCorrectorParameters> corrVec;

                for(std::string jecFile: Util2::GetVector<std::string>(sf, (!isData ? "Jet.JEC.MC." : "Jet.JEC.DATA.") + era)){
                    if(run != "MC") jecFile.replace(jecFile.find("@"), 1, run);
                    jecFile.replace(jecFile.find("&"), 1, type);
                    jecFile = this->filePath + jecFile;

                    if(!std::experimental::filesystem::exists(jecFile)){
                        throw std::runtime_error("File not exists" + jecFile);
                    }

                    corrVec.push_back(JetCorrectorParameters(jecFile));
                }

                if(type == "AK4") jetCorrectorAK4 = std::make_shared<FactorizedJetCorrector>(corrVec);
                else jetCorrectorAK8 = std::make_shared<FactorizedJetCorrector>(corrVec);

                //Set object to get JEC uncertainty
                std::string jecUncFile = this->filePath + sf.get<std::string>("Jet.JECUNC." + era);
                jecUncFile.replace(jecUncFile.find("&"), 1, type);  

                if(type == "AK4") jecUncAK4 = std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecUncFile, "Total"));
                else jecUncAK8 = std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecUncFile, "Total"));

                //Class for JME
                std::string JMEResoFile = this->filePath + sf.get<std::string>("Jet.JMEPtReso." + era);
                JMEResoFile.replace(JMEResoFile.find("&"), 1, type);
                if(type == "AK4") resolutionAK4 = JME::JetResolution(JMEResoFile);
                else resolutionAK8 = JME::JetResolution(JMEResoFile);

                std::string JMEFile = this->filePath + sf.get<std::string>("Jet.JME." + era);
                JMEFile.replace(JMEFile.find("&"), 1, type);
                if(type == "AK4") resolutionSFAK4 = JME::JetResolutionScaleFactor(JMEFile);
                else resolutionSFAK8 = JME::JetResolutionScaleFactor(JMEFile);
            }

            //BTag cuts
            float looseDeepCSVThr = skim.get<float>("Analyzer.Jet.btag.DeepCSV." + era + ".loose");
            float mediumDeepCSVThr = skim.get<float>("Analyzer.Jet.btag.DeepCSV." + era + ".medium");
            float tightDeepCSVThr = skim.get<float>("Analyzer.Jet.btag.DeepCSV." + era + ".tight");

            float looseDeepJetThr = skim.get<float>("Analyzer.Jet.btag.DeepJet." + era + ".loose");
            float mediumDeepJetThr = skim.get<float>("Analyzer.Jet.btag.DeepJet." + era + ".medium");
            float tightDeepJetThr = skim.get<float>("Analyzer.Jet.btag.DeepJet." + era + ".tight");

            isDeepCSVLoose = [=](const float& s){return s > looseDeepCSVThr;};
            isDeepCSVMedium = [=](const float& s){return s > mediumDeepCSVThr;};
            isDeepCSVTight = [=](const float& s){return s > tightDeepCSVThr;};

            isDeepJetLoose = [=](const float& s){return s > looseDeepJetThr;};
            isDeepJetMedium = [=](const float& s){return s > mediumDeepJetThr;};
            isDeepJetTight = [=](const float& s){return s > tightDeepJetThr;};
        }

        void Analyze(T& input, Output& out){
            out.nJets = 0, out.nSubJets = 0, out.nFatJets = 0;
            input.ReadJetEntry(isData);
            
            //MET
            float metPx = input.metPt * std::cos(input.metPhi);
            float metPy = input.metPt * std::sin(input.metPhi);

            //Same needed definition
            float jetJEC = 1., fatJetJEC = 1.,
                  jetJME = 1., fatJetJME = 1.;
        
            //Loop over all fat jets
            for(int i = 0; i < input.fatJetSize; ++i){
                if(out.nFatJets >= fatJetMax) break;
                input.GetFatJet(i);

                fatJetJEC = CorrectEnergy(input.fatJetPtRaw, input.fatJetEta, input.fatJetPhi, input.fatJetArea, input.rho, false);

                if(isJECSyst){
                    jecUncAK8->setJetPt(fatJetJEC*input.fatJetPtRaw);
                    jecUncAK8->setJetEta(input.fatJetEta);
                    jecUncAK8->setJetPhi(input.fatJetPhi);

                    fatJetJEC *= shift == "Up" ? 1 + jecUncAK8->getUncertainty(1) : 1 - jecUncAK8->getUncertainty(0);
                }

            
                if(!isData) fatJetJME = SmearEnergy(input.fatJetPtRaw*fatJetJEC, input.fatJetEta, input.fatJetPhi, input.rho, input, false);

                if(input.fatJetPtRaw*fatJetJEC*fatJetJME > 170. and input.fatJetMassRaw*fatJetJEC*fatJetJME > 40. and std::abs(input.fatJetEta) < etaCut){
                    out.fatJetPt[out.nFatJets] = input.fatJetPtRaw*fatJetJEC*fatJetJME;
                    out.fatJetMass[out.nFatJets] =  input.fatJetMassRaw*fatJetJEC*fatJetJME;
                    out.fatJetPhi[out.nFatJets] =  input.fatJetPhi;
                    out.fatJetEta[out.nFatJets] =  input.fatJetEta;
                    out.fatJetTau1[out.nFatJets] =  input.fatJetTau1;
                    out.fatJetTau2[out.nFatJets] =  input.fatJetTau2;
                    out.fatJetTau3[out.nFatJets] =  input.fatJetTau3;
                    out.fatJetDAK8ID[out.nFatJets] =  input.fatJetDAK8ID;
                    out.fatJetJEC[out.nFatJets] = fatJetJEC;
                    out.fatJetJME[out.nFatJets] = fatJetJME;

                    ++out.nFatJets;
                };
            }

            //Loop over all jets
            for(int i = 0; i < input.jetSize; ++i){
                if(out.nJets >= jetMax) break;
                input.GetJet(i);

                jetJEC = CorrectEnergy(input.jetPtRaw, input.jetEta, input.jetPhi, input.jetArea, input.rho, true);

                if(isJECSyst){
                    jecUncAK4->setJetPt(jetJEC*input.jetPtRaw);
                    jecUncAK4->setJetEta(input.jetEta);
                    jecUncAK4->setJetPhi(input.jetPhi);

                    jetJEC *= shift == "Up" ? 1 + jecUncAK4->getUncertainty(1) : 1 - jecUncAK4->getUncertainty(0);
                }

                if(!isData) jetJME = SmearEnergy(input.jetPtRaw*jetJEC, input.jetEta, input.jetPhi, input.rho, input, true);

                if(input.jetPtRaw*jetJEC*jetJME > ptCut and std::abs(input.jetEta) < etaCut){
                    //Propagate JEC/JME to met
                    metPx += input.jetPt * std::cos(input.jetPhi) - input.jetPtRaw*jetJEC*jetJME * std::cos(input.jetPhi);
                    metPy += input.jetPt * std::sin(input.jetPhi) - input.jetPtRaw*jetJEC*jetJME * std::sin(input.jetPhi);

                    out.fatJetIdx[out.nSubJets] = -1;

                    for(int j = 0; j < out.nFatJets; ++j){
                        if(Util2::DeltaR(out.fatJetEta[j], out.fatJetPhi[j], input.jetEta, input.jetPhi) < 1.2){
                            out.fatJetIdx[out.nSubJets] = j;
                            break;
                        }
                    }

                    if(out.fatJetIdx[out.nSubJets] == -1){
                        out.jetPt[out.nJets] = input.jetPtRaw*jetJEC*jetJME;
                        out.jetMass[out.nJets] =  input.jetMassRaw*jetJEC*jetJME;
                        out.jetPhi[out.nJets] =  input.jetPhi;
                        out.jetEta[out.nJets] =  input.jetEta;
                        out.jetDeepJet[out.nJets] =  input.jetDeepJet;
                        out.jetDeepJetID[out.nJets] =  isDeepJetLoose(input.jetDeepJet) + isDeepJetMedium(input.jetDeepJet) +  isDeepJetTight(input.jetDeepJet);
                        out.jetDeepCSV[out.nJets] =  input.jetDeepCSV;
                        out.jetDeepCSVID[out.nJets] =  isDeepCSVLoose(input.jetDeepCSV) + isDeepCSVMedium(input.jetDeepCSV) +  isDeepCSVTight(input.jetDeepCSV);
                        out.jetPartFlav[out.nJets] = input.jetPartFlav;
                        out.jetJEC[out.nJets] = jetJEC;
                        out.jetJME[out.nJets] = jetJME;

                        ++out.nJets;
                    }

                    else{
                        out.subJetPt[out.nSubJets] = input.jetPtRaw*jetJEC*jetJME;
                        out.subJetMass[out.nSubJets] =  input.jetMassRaw*jetJEC*jetJME;
                        out.subJetPhi[out.nSubJets] =  input.jetPhi;
                        out.subJetEta[out.nSubJets] =  input.jetEta;
                        out.subJetDeepJet[out.nSubJets] =  input.jetDeepJet;
                        out.subJetDeepJetID[out.nSubJets] =  isDeepJetLoose(input.jetDeepJet) + isDeepJetMedium(input.jetDeepJet) +  isDeepJetTight(input.jetDeepJet);
                        out.subJetDeepCSV[out.nSubJets] =  input.jetDeepCSV;
                        out.subJetDeepCSVID[out.nSubJets] =  isDeepCSVLoose(input.jetDeepCSV) + isDeepCSVMedium(input.jetDeepCSV) +  isDeepCSVTight(input.jetDeepCSV);
                        out.subJetPartFlav[out.nSubJets] = input.jetPartFlav;
                        out.subJetJEC[out.nSubJets] = jetJEC;
                        out.subJetJME[out.nSubJets] = jetJME;

                        ++out.nSubJets;
                    }
                }
            }

            out.metPt = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
            out.metPtDown = std::sqrt(std::pow(metPx - input.metDeltaUnClustX, 2) + std::pow(metPy - input.metDeltaUnClustY, 2));
            out.metPtUp = std::sqrt(std::pow(metPx + input.metDeltaUnClustX, 2) + std::pow(metPy + input.metDeltaUnClustY, 2));
            out.metPhi = std::atan2(metPy, metPx);
            out.metPhiDown = std::atan2(metPy - input.metDeltaUnClustY, metPx - input.metDeltaUnClustX);
            out.metPhiUp = std::atan2(metPy + input.metDeltaUnClustY, metPx + - input.metDeltaUnClustX);
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
