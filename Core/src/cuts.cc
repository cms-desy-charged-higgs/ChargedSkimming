#include <ChargedSkimming/Core/interface/cuts.h>
#include <iostream>

Cuts::Cuts(const std::shared_ptr<TFile>& outFile, const std::string& channel){
    cutFlow = std::make_shared<TH1F>();
    cutFlow->SetName(("Cutflow_" + channel).c_str());
    cutFlow->SetTitle(("Cutflow_" + channel).c_str());
    cutFlow->SetDirectory(outFile.get());
}

std::function<bool()> Cuts::ConstructCut(short& value, const std::string& op, const short& threshold){
    if(op == "=="){
        return [&value, threshold](){return value == threshold;};
    }

    else if(op == ">="){
        return [&value, threshold](){return value >= threshold;};
    } 

    else if(op == "<="){
        return [&value, threshold](){return value <= threshold;};
    } 

    else throw std::runtime_error(("Unknow cut operator: '" + op + "'").c_str());
}

void Cuts::AddCut(const std::string& part, Output& out, const std::string& op, const short& threshold){
    if(part == "Electron"){
        cuts.push_back(ConstructCut(out.nElectrons, op, threshold));
        cutNames.push_back("N_{e} " + op + std::to_string(threshold) + " (No ID.)");
    }

    else if(part == "Muon"){
        cuts.push_back(ConstructCut(out.nMuons, op, threshold));
        cutNames.push_back("N_{#mu} " + op + std::to_string(threshold) + " (No ID.)");
    }

    else if(part == "Jet"){
        cuts.push_back(ConstructCut(out.nJets, op, threshold));
        cutNames.push_back("N_{j} " + op + std::to_string(threshold) + " (Not clean)");
    }

    else if(part == "FatJet"){
        cuts.push_back(ConstructCut(out.nFatJets, op, threshold));
        cutNames.push_back("N_{fj} " + op + std::to_string(threshold));
    }

    else throw std::runtime_error(("Unknow particle: '" + op + "'").c_str());
}

bool Cuts::Passed(){
    for(std::size_t i = 0; i < cuts.size(); ++i){
        if(!cuts[i]()) return false;
    }

    return true;
}

void Cuts::FillCutflow(){
    for(std::size_t i = 0; i < cuts.size(); ++i){
        if(!cuts[i]()) return;
        cutFlow->Fill(cutNames[i].c_str(), 1);
    }
}
