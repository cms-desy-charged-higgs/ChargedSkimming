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
#include <ChargedSkimming/Skimming/interface/btagcsvreader.h>
#include <ChargedSkimming/Skimming/interface/util.h>

template <typename T>
class JetAnalyzer : public BaseAnalyzer<T> {
    private:
        //Input information
        std::string era, run;
        bool isData;

        //Kinematic cut criteria
        float ptCut, etaCut;

        //JEC corrector
        std::shared_ptr<FactorizedJetCorrector> jetCorrectorAK4, jetCorrectorAK8;

        //JEC uncertainty
        std::vector<std::string> JECSysts;
        std::vector<std::shared_ptr<JetCorrectionUncertainty>> jecUncAK4, jecUncAK8;

        //Classes for reading jet energy SF 
        JME::JetParameters jetParameter;
        JME::JetResolution resolutionAK4, resolutionAK8;
        JME::JetResolutionScaleFactor resolutionSFAK4, resolutionSFAK8;

        std::default_random_engine generator;

        //BTag cuts
        std::function<bool(const float&)> isDeepCSVLoose, isDeepCSVMedium, isDeepCSVTight,
                                          isDeepJetLoose, isDeepJetMedium, isDeepJetTight;

        template <typename V>
        void resortByIndex(V& arr, const std::vector<int>& index){
            V tmp = arr;
                    
            for(int i = 0; i < index.size(); ++i){
                arr.at(i) = tmp[index.at(i)];
            } 
        }

    public:
        JetAnalyzer(){}

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
        std::tuple<float, float, float> SmearEnergy(const float& pt, const float& eta, const float& phi, const float& rho, T& input, const bool& isAK4){
            //Configure smearer class
            jetParameter.setJetPt(pt).setJetEta(eta).setRho(rho);
            JME::JetResolution& resolution = isAK4 ? resolutionAK4 : resolutionAK8;
            JME::JetResolutionScaleFactor& resolutionSF = isAK4 ? resolutionSFAK4 : resolutionSFAK8;

            float reso = resolution.getResolution(jetParameter);
            float resoSF = resolutionSF.getScaleFactor(jetParameter);
            float resoSFUp = resolutionSF.getScaleFactor(jetParameter, Variation::UP);
            float resoSFDown = resolutionSF.getScaleFactor(jetParameter, Variation::DOWN);
            
            float JME = 1., JMEUp = 1., JMEDown = 1., 
                  genPt = -1., 
                  dRmin = std::numeric_limits<float>::max(), 
                  coneSize = isAK4 ? 0.2 : 0.4;

            int size = isAK4 ? input.genJetSize : input.genFatJetSize;

            //Gen jet matching
            for(int i = 0; i < size; ++i){
                float dR, dPT, genJetPt;

                if(isAK4){
                    input.GetGenJet(i);
                    dR = Util::DeltaR(eta, phi, input.genJetEta, input.genJetPhi);
                    genJetPt = input.genJetPt;
                    dPT = std::abs(pt - genJetPt);
                }

                else{
                    input.GetGenFatJet(i);
                    dR = Util::DeltaR(eta, phi, input.genFatJetEta, input.genFatJetPhi);
                    genJetPt = input.genFatJetPt;
                    dPT = std::abs(pt - genJetPt);
                }

                if(dR > dRmin) continue;
           
                if(dR < coneSize and dPT < 3.*reso*pt){
                    dRmin = dR;
                    genPt = genJetPt;
                }
            }

            //If gen matcheds
            if(genPt > 0.){
                JME = 1. + (resoSF-1.)*(pt - genPt)/pt; 
                JMEUp = 1. + (resoSFUp-1.)*(pt - genPt)/pt; 
                JMEDown = 1. + (resoSFDown-1.)*(pt - genPt)/pt; 
            }

            else {
                if(resoSF > 1.){
                    std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));
                    JME = 1. + gaus(generator);
                } 
                    
                if(resoSFUp > 1.){ 
                    std::normal_distribution<> gausUp(0, reso * std::sqrt(resoSFUp * resoSFUp - 1));
                    JMEUp = 1. + gausUp(generator);
                }    
                    
                if(resoSFDown > 1.){ 
                    std::normal_distribution<> gausDown(0, reso * std::sqrt(resoSFDown * resoSFDown - 1));
                    JMEDown = 1. + gausDown(generator);
                }
            }

            //Check if direction of jet not changed
            if(pt*JME < 1e-2){
                JME = 1e-2/pt;
            }

            if(pt*JMEUp < 1e-2){
                JMEUp = 1e-2/pt;
            }

            if(pt*JMEDown < 1e-2){
                JMEDown = 1e-2/pt;
            }

            return {JME, JMEUp, JMEDown};
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

            ptCut = skim.get<float>("Analyzer.Jet.pt." + era);
            etaCut = skim.get<float>("Analyzer.Jet.eta." + era);

            //Set JEC/JME classes
            for(const std::string type : {"AK4", "AK8"}){
                //Class for JEC
                std::vector<JetCorrectorParameters> corrVec;

                for(std::string jecFile: Util::GetVector<std::string>(sf, (!isData ? "Jet.JEC.MC." : "Jet.JEC.DATA.") + era)){
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
                JECSysts = !isData ? Util::GetVector<std::string>(skim, "Analyzer.Jet.JECSyst") : std::vector<std::string>{};
                
                std::string jecUncFile = this->filePath + sf.get<std::string>("Jet.JECUNC." + era);
                jecUncFile.replace(jecUncFile.find("&"), 1, type);  

                for(std::string& JECSyst: JECSysts){
                    if(type == "AK4") jecUncAK4.push_back(std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecUncFile, JECSyst)));
                    else jecUncAK8.push_back(std::make_shared<JetCorrectionUncertainty>(JetCorrectorParameters(jecUncFile, JECSyst)));
                }

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
            if(!isData) input.ReadGenEntry();
            
            //MET
            float metPx = input.metPt * std::cos(input.metPhi);
            std::vector<float> metPxJECUp(JECSysts.size(), input.metPt * std::cos(input.metPhi));
            std::vector<float> metPxJECDown(JECSysts.size(), input.metPt * std::cos(input.metPhi));
            float metPxJMEUp = input.metPt * std::cos(input.metPhi);
            float metPxJMEDown = input.metPt * std::cos(input.metPhi);
            float metPy = input.metPt * std::sin(input.metPhi);
            std::vector<float> metPyJECUp(JECSysts.size(), input.metPt * std::sin(input.metPhi));
            std::vector<float> metPyJECDown(JECSysts.size(), input.metPt * std::sin(input.metPhi));
            float metPyJMEUp = input.metPt * std::sin(input.metPhi);
            float metPyJMEDown = input.metPt * std::sin(input.metPhi);

            //Same needed definition
            float jetJEC = 1.,
                  jetJME = 1., jetJMEUp = 1., jetJMEDown = 1.,
                  fatJetJEC = 1.,
                  fatJetJME = 1., fatJetJMEUp = 1., fatJetJMEDown = 1.;

            int fIdx = -1;
                  
            std::vector<float> jetJECUp(JECSysts.size(), 1.), jetJECDown(JECSysts.size(), 1.),
                               fatJetJECUp(JECSysts.size(), 1.), fatJetJECDown(JECSysts.size(), 1.);
        
            //Loop over all fat jets
            for(int i = 0; i < input.fatJetSize; ++i){
                if(out.nFatJets >= fatJetMax) break;
                input.GetFatJet(i);

                fatJetJEC = CorrectEnergy(input.fatJetPtRaw, input.fatJetEta, input.fatJetPhi, input.fatJetArea, input.rho, false);

                for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                    jecUncAK8.at(JEC)->setJetPt(fatJetJEC*input.fatJetPtRaw);
                    jecUncAK8.at(JEC)->setJetEta(input.fatJetEta);
                    jecUncAK8.at(JEC)->setJetPhi(input.fatJetPhi);
                
                    fatJetJECUp.at(JEC) = fatJetJEC*(1 + jecUncAK8.at(JEC)->getUncertainty(1));
                    
                    jecUncAK8.at(JEC)->setJetPt(fatJetJEC*input.fatJetPtRaw);
                    jecUncAK8.at(JEC)->setJetEta(input.fatJetEta);
                    jecUncAK8.at(JEC)->setJetPhi(input.fatJetPhi);
                    
                    fatJetJECDown.at(JEC) = fatJetJEC*(1 - jecUncAK8.at(JEC)->getUncertainty(0));
                }

                if(!isData) std::tie(fatJetJME, fatJetJMEUp, fatJetJMEDown) = SmearEnergy(input.fatJetPtRaw*fatJetJEC, input.fatJetEta, input.fatJetPhi, input.rho, input, false);
                
                float maxJEC = isData ? fatJetJEC : std::max(fatJetJEC, std::max(*std::max_element(fatJetJECDown.begin(), fatJetJECDown.end()), *std::max_element(fatJetJECUp.begin(), fatJetJECUp.end())));
                float maxJME = isData ? 1. : std::max(fatJetJME, std::max(fatJetJMEDown, fatJetJMEUp));

                if(input.fatJetPtRaw*maxJEC*maxJME > 170. and input.fatJetMassRaw*maxJEC*maxJME > 40. and std::abs(input.fatJetEta) < etaCut){
                    out.fatJetPt.at(out.nFatJets) = input.fatJetPtRaw*fatJetJEC*fatJetJME;
                    out.fatJetPtJMEUp.at(out.nFatJets) = input.fatJetPtRaw*fatJetJEC*fatJetJMEUp;
                    out.fatJetPtJMEDown.at(out.nFatJets) = input.fatJetPtRaw*fatJetJEC*fatJetJMEDown;
                    out.fatJetMass.at(out.nFatJets) =  input.fatJetMassRaw*fatJetJEC*fatJetJME;
                    out.fatJetMassJMEUp.at(out.nFatJets) = input.fatJetMassRaw*fatJetJEC*fatJetJMEUp;
                    out.fatJetMassJMEDown.at(out.nFatJets) = input.fatJetMassRaw*fatJetJEC*fatJetJMEDown;
                    out.fatJetPhi.at(out.nFatJets) =  input.fatJetPhi;
                    out.fatJetEta.at(out.nFatJets) =  input.fatJetEta;
                    out.fatJetTau1.at(out.nFatJets) =  input.fatJetTau1;
                    out.fatJetTau2.at(out.nFatJets) =  input.fatJetTau2;
                    out.fatJetTau3.at(out.nFatJets) =  input.fatJetTau3;
                    out.fatJetDAK8ID.at(out.nFatJets) =  input.fatJetDAK8ID;
                    out.fatJetJEC.at(out.nFatJets) = fatJetJEC;
                    out.fatJetJME.at(out.nFatJets) = fatJetJME;

                    for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                        out.fatJetPtJECUp.at(JEC).at(out.nFatJets) = input.fatJetPtRaw*fatJetJECUp.at(JEC)*fatJetJME;
                        out.fatJetPtJECDown.at(JEC).at(out.nFatJets) = input.fatJetPtRaw*fatJetJECDown.at(JEC)*fatJetJME;
                        
                        out.fatJetMassJECUp.at(JEC).at(out.nFatJets) = input.fatJetMassRaw*fatJetJECUp.at(JEC)*fatJetJME;
                        out.fatJetMassJECDown.at(JEC).at(out.nFatJets) = input.fatJetMassRaw*fatJetJECDown.at(JEC)*fatJetJME;   
                    }

                    ++out.nFatJets;
                };
            }

            //Loop over all jets
            for(int i = 0; i < input.jetSize; ++i){
                if(out.nJets >= jetMax) break;
                input.GetJet(i);

                jetJEC = CorrectEnergy(input.jetPtRaw, input.jetEta, input.jetPhi, input.jetArea, input.rho, true);

                for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                    jecUncAK4.at(JEC)->setJetPt(jetJEC*input.jetPtRaw);
                    jecUncAK4.at(JEC)->setJetEta(input.jetEta);
                    jecUncAK4.at(JEC)->setJetPhi(input.jetPhi);

                    jetJECUp.at(JEC) = jetJEC*(1 + jecUncAK4.at(JEC)->getUncertainty(1));

                    jecUncAK4.at(JEC)->setJetPt(jetJEC*input.jetPtRaw);
                    jecUncAK4.at(JEC)->setJetEta(input.jetEta);
                    jecUncAK4.at(JEC)->setJetPhi(input.jetPhi);

                    jetJECDown.at(JEC) = jetJEC*(1 - jecUncAK4.at(JEC)->getUncertainty(0));
                }

                if(!isData) std::tie(jetJME, jetJMEUp, jetJMEDown) = SmearEnergy(input.jetPtRaw*jetJEC, input.jetEta, input.jetPhi, input.rho, input, true);

                float maxJEC = isData ? jetJEC : std::max(jetJEC, std::max(*std::max_element(jetJECDown.begin(), jetJECDown.end()), *std::max_element(jetJECUp.begin(), jetJECUp.end())));
                float maxJME = isData ? 1. : std::max(jetJME, std::max(jetJMEUp, jetJMEDown));

                if(input.jetPtRaw*maxJEC*maxJME > ptCut and std::abs(input.jetEta) < etaCut){
                    //Propagate JEC/JME to met
                    metPx += input.jetPt * std::cos(input.jetPhi) - input.jetPtRaw*jetJEC*jetJME * std::cos(input.jetPhi);
                    metPxJMEUp += input.jetPt * std::cos(input.jetPhi) - input.jetPtRaw*jetJEC*jetJMEUp * std::cos(input.jetPhi);
                    metPxJMEDown += input.jetPt * std::cos(input.jetPhi) - input.jetPtRaw*jetJEC*jetJMEDown * std::cos(input.jetPhi);
                    metPy += input.jetPt * std::sin(input.jetPhi) - input.jetPtRaw*jetJEC*jetJME * std::sin(input.jetPhi);
                    metPyJMEUp += input.jetPt * std::sin(input.jetPhi) - input.jetPtRaw*jetJEC*jetJMEUp * std::sin(input.jetPhi);
                    metPyJMEDown += input.jetPt * std::sin(input.jetPhi) - input.jetPtRaw*jetJEC*jetJMEDown * std::sin(input.jetPhi);
                    
                    for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                        metPxJECUp.at(JEC) += input.jetPt * std::cos(input.jetPhi) - input.jetPtRaw*jetJECUp.at(JEC)*jetJME * std::cos(input.jetPhi);
                        metPxJECDown.at(JEC) += input.jetPt * std::cos(input.jetPhi) - input.jetPtRaw*jetJECDown.at(JEC)*jetJME * std::cos(input.jetPhi);
                        
                        metPyJECUp.at(JEC) += input.jetPt * std::sin(input.jetPhi) - input.jetPtRaw*jetJECUp.at(JEC)*jetJME * std::sin(input.jetPhi);
                        metPyJECDown.at(JEC) += input.jetPt * std::sin(input.jetPhi) - input.jetPtRaw*jetJECDown.at(JEC)*jetJME * std::sin(input.jetPhi);
                    }

                    fIdx = -1;

                    for(int j = 0; j < out.nFatJets; ++j){
                        if(Util::DeltaR(out.fatJetEta[j], out.fatJetPhi[j], input.jetEta, input.jetPhi) < 1.2){
                            fIdx = j;
                            break;
                        }
                    }

                    if(fIdx == -1){
                        out.jetPt.at(out.nJets) = input.jetPtRaw*jetJEC*jetJME;
                        out.jetPtJMEUp.at(out.nJets) = input.jetPtRaw*jetJEC*jetJMEUp;
                        out.jetPtJMEDown.at(out.nJets) = input.jetPtRaw*jetJEC*jetJMEDown;
                        out.jetMass.at(out.nJets) =  input.jetMassRaw*jetJEC*jetJME;
                        out.jetMassJMEUp.at(out.nJets) = input.jetMassRaw*jetJEC*jetJMEUp;
                        out.jetMassJMEDown.at(out.nJets) = input.jetMassRaw*jetJEC*jetJMEDown;
                        out.jetPhi.at(out.nJets) =  input.jetPhi;
                        out.jetEta.at(out.nJets) =  input.jetEta;
                        out.jetDeepJet.at(out.nJets) =  input.jetDeepJet;
                        out.jetDeepJetID.at(out.nJets) =  isDeepJetLoose(input.jetDeepJet) + isDeepJetMedium(input.jetDeepJet) +  isDeepJetTight(input.jetDeepJet);
                        out.jetDeepCSV.at(out.nJets) =  input.jetDeepCSV;
                        out.jetDeepCSVID.at(out.nJets) =  isDeepCSVLoose(input.jetDeepCSV) + isDeepCSVMedium(input.jetDeepCSV) +  isDeepCSVTight(input.jetDeepCSV);
                        out.jetPartFlav.at(out.nJets) = input.jetPartFlav;
                        out.jetJEC.at(out.nJets) = jetJEC;
                        out.jetJME.at(out.nJets) = jetJME;
                        out.jetID.at(out.nJets) = input.jetID;
                        out.jetPUID.at(out.nJets) = input.jetPUID;
                        
                        for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                            out.jetPtJECUp.at(JEC).at(out.nJets) = input.jetPtRaw*jetJECUp.at(JEC)*jetJME;
                            out.jetPtJECDown.at(JEC).at(out.nJets) = input.jetPtRaw*jetJECDown.at(JEC)*jetJME;
                            
                            out.jetMassJECUp.at(JEC).at(out.nJets) = input.jetMassRaw*jetJECUp.at(JEC)*jetJME;
                            out.jetMassJECDown.at(JEC).at(out.nJets) = input.jetMassRaw*jetJECDown.at(JEC)*jetJME;
                        }

                        if(!isData){
                            int genIdx = input.GenMatch(input, out.jetPt.at(out.nJets), out.jetPhi.at(out.nJets), out.jetEta.at(out.nJets), 5, 0.4, 3.);

                            if(genIdx != -1){
                                input.alreadyMatchedIdx.push_back(genIdx);

                                input.GetGenPart(genIdx);
                                out.jetGenPt.at(out.nJets) = input.genPt;
                                out.jetGenPhi.at(out.nJets) = input.genPhi;
                                out.jetGenEta.at(out.nJets) = input.genEta;
                                out.jetGenID.at(out.nJets) = input.genPDG;

                                input.GetGenPart(input.genMotherIdx);
                                out.jetGenMotherID.at(out.nJets) = input.genPDG;

                                input.GetGenPart(input.genMotherIdx);
                                out.jetGenGrandMotherID.at(out.nJets) = input.genPDG;
                            }

                            else{
                                out.jetGenPt.at(out.nJets) = -999.;
                                out.jetGenPhi.at(out.nJets) = -999.;
                                out.jetGenEta.at(out.nJets) = -999.;
                                out.jetGenID.at(out.nJets) = -999;
                                out.jetGenMotherID.at(out.nJets) = -999;
                                out.jetGenGrandMotherID.at(out.nJets) = -999;
                            }
                        }

                        ++out.nJets;
                    }

                    else{
                        if(out.nSubJets >= jetMax) continue;

                        out.fatJetIdx.at(out.nSubJets) = fIdx;
                        out.subJetPt.at(out.nSubJets) = input.jetPtRaw*jetJEC*jetJME;
                        out.subJetPtJMEUp.at(out.nSubJets) = input.jetPtRaw*jetJEC*jetJMEUp;
                        out.subJetPtJMEDown.at(out.nSubJets) = input.jetPtRaw*jetJEC*jetJMEDown;
                        out.subJetMass.at(out.nSubJets) =  input.jetMassRaw*jetJEC*jetJME;
                        out.subJetMassJMEUp.at(out.nSubJets) = input.jetMassRaw*jetJEC*jetJMEUp;
                        out.subJetMassJMEDown.at(out.nSubJets) = input.jetMassRaw*jetJEC*jetJMEDown;
                        out.subJetPhi.at(out.nSubJets) =  input.jetPhi;
                        out.subJetEta.at(out.nSubJets) =  input.jetEta;
                        out.subJetDeepJet.at(out.nSubJets) =  input.jetDeepJet;
                        out.subJetDeepJetID.at(out.nSubJets) =  isDeepJetLoose(input.jetDeepJet) + isDeepJetMedium(input.jetDeepJet) +  isDeepJetTight(input.jetDeepJet);
                        out.subJetDeepCSV.at(out.nSubJets) =  input.jetDeepCSV;
                        out.subJetDeepCSVID.at(out.nSubJets) =  isDeepCSVLoose(input.jetDeepCSV) + isDeepCSVMedium(input.jetDeepCSV) +  isDeepCSVTight(input.jetDeepCSV);
                        out.subJetPartFlav.at(out.nSubJets) = input.jetPartFlav;
                        out.subJetJEC.at(out.nSubJets) = jetJEC;
                        out.subJetJME.at(out.nSubJets) = jetJME;

                        for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                            out.subJetPtJECUp.at(JEC).at(out.nSubJets) = input.jetPtRaw*jetJECUp.at(JEC)*jetJME;
                            out.subJetPtJECDown.at(JEC).at(out.nSubJets) = input.jetPtRaw*jetJECDown.at(JEC)*jetJME;
                            
                            out.subJetMassJECUp.at(JEC).at(out.nSubJets) = input.jetMassRaw*jetJECUp.at(JEC)*jetJME;
                            out.subJetMassJECDown.at(JEC).at(out.nSubJets) = input.jetMassRaw*jetJECDown.at(JEC)*jetJME;   
                        }

                        if(!isData){
                            int genIdx = input.GenMatch(input, out.subJetPt.at(out.nSubJets), out.subJetPhi.at(out.nSubJets), out.subJetEta.at(out.nSubJets), 5, 0.4, 3.);

                            if(genIdx != -1){
                                input.alreadyMatchedIdx.push_back(genIdx);

                                input.GetGenPart(genIdx);
                                out.subJetGenPt.at(out.nSubJets) = input.genPt;
                                out.subJetGenPhi.at(out.nSubJets) = input.genPhi;
                                out.subJetGenEta.at(out.nSubJets) = input.genEta;
                                out.subJetGenID.at(out.nSubJets) = input.genPDG;

                                input.GetGenPart(input.genMotherIdx);
                                out.subJetGenMotherID.at(out.nSubJets) = input.genPDG;

                                input.GetGenPart(input.genMotherIdx);
                                out.subJetGenGrandMotherID.at(out.nSubJets) = input.genPDG;
                            }

                            else{
                                out.subJetGenPt.at(out.nSubJets) = -999.;
                                out.subJetGenPhi.at(out.nSubJets) = -999.;
                                out.subJetGenEta.at(out.nSubJets) = -999.;
                                out.subJetGenID.at(out.nSubJets)= -999;
                                out.subJetGenMotherID.at(out.nSubJets) = -999;
                                out.subJetGenGrandMotherID.at(out.nSubJets) = -999;
                            }
                        }

                        ++out.nSubJets;
                    }
                }
            }

            out.metPt = std::sqrt(std::pow(metPx, 2) + std::pow(metPy, 2));
            out.metPtJMEUp = std::sqrt(std::pow(metPxJMEUp, 2) + std::pow(metPyJMEUp, 2));
            out.metPtJMEDown = std::sqrt(std::pow(metPxJMEDown, 2) + std::pow(metPyJMEDown, 2));
            out.metPtUnclusteredUp = std::sqrt(std::pow(metPx + input.metDeltaUnClustX, 2) + std::pow(metPy + input.metDeltaUnClustY, 2));
            out.metPtUnclusteredDown = std::sqrt(std::pow(metPx - input.metDeltaUnClustX, 2) + std::pow(metPy - input.metDeltaUnClustY, 2));

            out.metPhi = std::atan2(metPy, metPx);
            out.metPhiJMEUp = std::atan2(metPyJMEUp, metPxJMEUp);
            out.metPhiJMEDown = std::atan2(metPyJMEDown, metPxJMEDown);
            out.metPhiUnclusteredDown = std::atan2(metPy - input.metDeltaUnClustY, metPx - input.metDeltaUnClustX);
            out.metPhiUnclusteredUp = std::atan2(metPy + input.metDeltaUnClustY, metPx + input.metDeltaUnClustX);
            
            for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                out.metPtJECUp.at(JEC) = std::sqrt(std::pow(metPxJECUp.at(JEC), 2) + std::pow(metPyJECUp.at(JEC), 2));
                out.metPtJECDown.at(JEC) = std::sqrt(std::pow(metPxJECDown.at(JEC), 2) + std::pow(metPyJECDown.at(JEC), 2));
                            
                out.metPhiJECUp.at(JEC) = std::atan2(metPyJECUp.at(JEC), metPxJECUp.at(JEC));
                out.metPhiJECDown.at(JEC) = std::atan2(metPyJECDown.at(JEC), metPxJECDown.at(JEC));
            }
            
            std::vector<int> jetIdx(out.nJets, 0), subJetIdx(out.nSubJets, 0), fatJetIdx(out.nFatJets, 0);
            std::iota(jetIdx.begin(), jetIdx.end(), 0);
            std::iota(subJetIdx.begin(), subJetIdx.end(), 0);
            std::iota(fatJetIdx.begin(), fatJetIdx.end(), 0);
            
            std::sort(jetIdx.begin(), jetIdx.end(), [&](const int& i1, const int& i2){return out.jetPt[i1] > out.jetPt[i2];});
            std::sort(subJetIdx.begin(), subJetIdx.end(), [&](const int& i1, const int& i2){return out.subJetPt[i1] > out.subJetPt[i2];});
            std::sort(fatJetIdx.begin(), fatJetIdx.end(), [&](const int& i1, const int& i2){return out.fatJetPt[i1] > out.fatJetPt[i2];});

            for(int JEC = 0; JEC < JECSysts.size(); ++JEC){
                resortByIndex(out.jetPtJECUp.at(JEC), jetIdx);
                resortByIndex(out.jetPtJECDown.at(JEC), jetIdx);
                resortByIndex(out.jetMassJECUp.at(JEC), jetIdx);
                resortByIndex(out.jetMassJECDown.at(JEC), jetIdx);

                resortByIndex(out.subJetPtJECUp.at(JEC), subJetIdx);
                resortByIndex(out.subJetPtJECDown.at(JEC), subJetIdx);
                resortByIndex(out.subJetMassJECUp.at(JEC), subJetIdx);
                resortByIndex(out.subJetMassJECDown.at(JEC), subJetIdx);

                resortByIndex(out.fatJetPtJECUp.at(JEC), fatJetIdx);
                resortByIndex(out.fatJetPtJECDown.at(JEC), fatJetIdx);
                resortByIndex(out.fatJetMassJECUp.at(JEC), fatJetIdx);
                resortByIndex(out.fatJetMassJECDown.at(JEC), fatJetIdx);
            }
            
            resortByIndex(out.jetPt, jetIdx);
            resortByIndex(out.jetPtJMEUp, jetIdx);
            resortByIndex(out.jetPtJMEDown, jetIdx);
            resortByIndex(out.jetMass, jetIdx);
            resortByIndex(out.jetMassJMEUp, jetIdx);
            resortByIndex(out.jetMassJMEDown, jetIdx);
            resortByIndex(out.jetEta, jetIdx);
            resortByIndex(out.jetPhi, jetIdx);
            resortByIndex(out.jetDeepJet, jetIdx);
            resortByIndex(out.jetDeepCSV, jetIdx);
            resortByIndex(out.jetJEC, jetIdx);
            resortByIndex(out.jetJME, jetIdx);
            resortByIndex(out.jetGenPt, jetIdx);
            resortByIndex(out.jetGenEta, jetIdx);
            resortByIndex(out.jetGenPhi, jetIdx);
            resortByIndex(out.jetDeepJetID, jetIdx);
            resortByIndex(out.jetDeepCSVID, jetIdx);
            resortByIndex(out.jetPartFlav, jetIdx);
            resortByIndex(out.jetGenID, jetIdx);
            resortByIndex(out.jetGenMotherID, jetIdx);
            resortByIndex(out.jetGenGrandMotherID, jetIdx);
            resortByIndex(out.jetID, jetIdx);
            resortByIndex(out.jetPUID, jetIdx);

            resortByIndex(out.subJetPt, subJetIdx);
            resortByIndex(out.subJetPtJMEUp, subJetIdx);
            resortByIndex(out.subJetPtJMEDown, subJetIdx);
            resortByIndex(out.subJetMass, subJetIdx);
            resortByIndex(out.subJetMassJMEUp, subJetIdx);
            resortByIndex(out.subJetMassJMEDown, subJetIdx);
            resortByIndex(out.subJetEta, subJetIdx);
            resortByIndex(out.subJetPhi, subJetIdx);
            resortByIndex(out.subJetDeepJet, subJetIdx);
            resortByIndex(out.subJetDeepCSV, subJetIdx);
            resortByIndex(out.subJetJEC, subJetIdx);
            resortByIndex(out.subJetJME, subJetIdx);
            resortByIndex(out.subJetGenPt, subJetIdx);
            resortByIndex(out.subJetGenEta, subJetIdx);
            resortByIndex(out.subJetGenPhi, subJetIdx);
            resortByIndex(out.fatJetIdx, subJetIdx);
            resortByIndex(out.subJetDeepJetID, subJetIdx);
            resortByIndex(out.subJetDeepCSVID, subJetIdx);
            resortByIndex(out.subJetPartFlav, subJetIdx);
            resortByIndex(out.subJetGenID, subJetIdx);
            resortByIndex(out.subJetGenMotherID, subJetIdx);
            resortByIndex(out.subJetGenGrandMotherID, subJetIdx);

            resortByIndex(out.fatJetPt, fatJetIdx);
            resortByIndex(out.fatJetPtJMEUp, fatJetIdx);
            resortByIndex(out.fatJetPtJMEDown, fatJetIdx);
            resortByIndex(out.fatJetMass, fatJetIdx);
            resortByIndex(out.fatJetMassJMEUp, fatJetIdx);
            resortByIndex(out.fatJetMassJMEDown, fatJetIdx);
            resortByIndex(out.fatJetEta, fatJetIdx);
            resortByIndex(out.fatJetPhi, fatJetIdx);
            resortByIndex(out.fatJetTau1, fatJetIdx);
            resortByIndex(out.fatJetTau2, fatJetIdx);
            resortByIndex(out.fatJetTau3, fatJetIdx);
            resortByIndex(out.fatJetJEC, fatJetIdx);
            resortByIndex(out.fatJetJME, fatJetIdx);
            resortByIndex(out.fatJetDAK8ID, fatJetIdx);
        }

        void EndJob(const std::shared_ptr<TFile>& outFile){};
};

#endif
