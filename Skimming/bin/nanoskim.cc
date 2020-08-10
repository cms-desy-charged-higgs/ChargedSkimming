#include <ChargedSkimming/Skimming/interface/nanoskimmer.h>

#include <vector>
#include <string>

std::vector<std::string> SplitString(const std::string& splitString, const std::string& delimeter){
    //Function which handles splitting of string input
    std::vector<std::string> splittedString;
    std::string string;
    std::istringstream splittedStream(splitString);
    while (std::getline(splittedStream, string, delimeter.c_str()[0])){
        splittedString.push_back(string);
    }

    return splittedString;
}


int main(int argc, char* argv[]){
  //Extract informations of command line
    std::string fileName = std::string(argv[1]);
    bool isData = std::string(argv[2]) == "True" ? true : false;
    std::vector<std::string> channels = SplitString(std::string(argv[3]), " ");
    float xSec = std::stof(std::string(argv[4]));
    std::string outName = std::string(argv[5]);
    int era = std::stoi(argv[6]);   

    NanoSkimmer skimmer(fileName, outName, channels, xSec, era, isData);
    skimmer.EventLoop();
    skimmer.WriteOutput();
}


