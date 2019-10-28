#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

class MuonAnalyzer: public BaseAnalyzer{
    private:
        //Check if data
        bool isData;

        //Map for SF files
        std::map<int, std::string> isoSFfiles;
        std::map<int, std::string> triggerSFfiles;
        std::map<int, std::string> IDSFfiles;

        //Hist with scale factors
        TH2F* triggerSFhist;

        std::vector<TH2F*> IsoHist;
        std::vector<TH2F*> IDHist;

        //Kinematic cut criteria
        int era;
        float ptCut;
        float etaCut;

        //Vector with output varirables of the output tree
        std::vector<std::string> floatNames;
        std::vector<std::string> boolNames;

        std::vector<std::vector<float>> floatVariables;
        std::vector<std::vector<bool>> boolVariables;

        //EDM Token for MINIAOD analysis
        muToken muonToken;
        trigObjToken triggerObjToken; 
        genPartToken genParticleToken;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> muonPt;
        std::unique_ptr<TTreeReaderArray<float>> muonEta;
        std::unique_ptr<TTreeReaderArray<float>> muonPhi;
        std::unique_ptr<TTreeReaderArray<float>> muonIso;
        std::unique_ptr<TTreeReaderArray<int>> muonCharge;
        std::unique_ptr<TTreeReaderArray<bool>> muonLooseID;
        std::unique_ptr<TTreeReaderArray<bool>> muonTightID;
        std::unique_ptr<TTreeReaderArray<int>> muonGenIdx;

    public:
        MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader);
        MuonAnalyzer(const int &era, const float &ptCut, const float &etaCut, muToken &muonToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken);

        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
