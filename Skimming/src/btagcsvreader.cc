#include <ChargedSkimming/Skimming/interface/btagcsvreader.h>

BTagCSVReader::BTagCSVReader(const std::string& fileName){
    //Open file and check if it is open
    std::ifstream CSV(fileName);
    if(!CSV.is_open()) throw std::runtime_error("File not exists: " + fileName);

    //Header, ignofire it
    int lineNr = 0;
    std::string line;
    std::getline(CSV, line);

    //Loop over file
    while(std::getline(CSV, line)){
        //Split line and get respektive parameter
        std::vector<std::string> splittedLine = SplitString<std::string>(line, ",");

        std::string type = splittedLine.at(1), sysType = splittedLine.at(2); 
        int wp = std::stoi(splittedLine.at(0)), jetFlavour = std::stoi(splittedLine.at(3));
        std::string formula = splittedLine.at(10).substr(1, splittedLine.at(10).size()-2);
      //  float etaMin = std::stof(splittedLine.at(4)), etaMax = std::stof(splittedLine.at(5));
        float ptMin = std::stof(splittedLine.at(6)), ptMax = std::stof(splittedLine.at(7));

        if(type != "comb" or jetFlavour != 0) continue;
        if(type == "iterativefit") break;

        std::shared_ptr<TF1> func = std::make_shared<TF1>(("Func" + std::to_string(lineNr)).c_str(), formula.c_str(), ptMin, ptMax);
        ptRange[{wp, sysType == "central" ? 0 : sysType == "up" ? 1 : 2}].push_back({ptMin, ptMax});
        SF[{wp, sysType == "central" ? 0 : sysType == "up" ? 1 : 2}].push_back(std::move(func));

        lineNr++;
    }
}

double BTagCSVReader::Get(const float& pt, const int& wp){
    int i=0;
    
    for(const std::pair<float, float> range : ptRange[{wp, 0}]){
        if(range.first < pt and pt < range.second){
            return SF[{wp, 0}][i]->Eval(pt);
        }   

        i++;
    }
  
    return 1.;
}

double BTagCSVReader::GetUp(const float& pt, const int& wp){
    int i=0;
    
    for(const std::pair<float, float> range : ptRange[{wp, 1}]){
        if(range.first < pt and pt < range.second){
            return SF[{wp, 1}][i]->Eval(pt);
        }   

        i++;
    }

    return 1.;
}

double BTagCSVReader::GetDown(const float& pt, const int& wp){
    int i=0;
    
    for(const std::pair<float, float> range : ptRange[{wp, 2}]){
        if(range.first < pt and pt < range.second){
            return SF[{wp, 2}][i]->Eval(pt);
        }   

        i++;
    }

    return 1.;
}
