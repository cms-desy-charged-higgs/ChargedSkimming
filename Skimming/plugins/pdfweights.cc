#include <ChargedSkimming/Skimming/interface/pdfweights.h>

PDFWeights::PDFWeights(const edm::ParameterSet& iConfig) :
    LHETag(iConfig.getParameter<edm::InputTag>("LHE")),
    LHAID(iConfig.getParameter<int>("LHAID")) {
        consumes<LHERunInfoProduct, edm::InRun>(LHETag);
        LHEToken = consumes<LHEEventProduct>(LHETag);
        PDFToken = produces<std::vector<float>>("pdfVariations").setBranchAlias("pdfVariations");
        ScaleToken = produces<std::vector<float>>("scaleVariations").setBranchAlias("scaleVariations");
    }

PDFWeights::~PDFWeights(){}

std::vector<std::string> PDFWeights::SplitString(const std::string& splitString, const std::string& delimeter){
    //Function which handles splitting of string input
    std::vector<std::string> splittedString;
    std::string string;
    std::istringstream splittedStream(splitString);
    while (std::getline(splittedStream, string, delimeter.c_str()[0])){
        splittedString.push_back(string);
    }

    return splittedString;
}

void PDFWeights::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){
    //Get information of LHE files about pdf weights
    edm::Handle<LHERunInfoProduct> run; 
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
 
    iRun.getByLabel(LHETag, run);
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
 
    for(headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
        std::vector<std::string> lines = iter->lines();
    
        //Loop over each line in the LHE file information
        for(unsigned int iLine = 0; iLine<lines.size(); iLine++) {
            //102 pdf variations
            for(int i = 1; i <= 102; i++){
                //Check if line contains wished information with for PDF with LHAID
                if(lines.at(iLine).find(std::to_string(LHAID +  i)) !=  std::string::npos){
                    //Split line in strings by empty space
                    std::vector<std::string> splittedLine = this->SplitString(lines.at(iLine), " ");

                    //Search for 'id="some number"' and extract number
                    for(const std::string& lineInfo: splittedLine){
                        std::size_t idStart = lineInfo.find("id=\"");

                        if(idStart != std::string::npos){
                            std::string partialNumber = lineInfo.substr(idStart+4);
                            std::string number = partialNumber.substr(0, partialNumber.find("\""));
                            variationIDs.push_back(number);
                        }
                    }
                }
            }

            //Renormalization and factorizlation scale
            std::vector<std::vector<std::string>> possibleScales = {
                {"1", "0.5"},
                {"1", "2"},
                {"2", "0.5"},
            };

            if(scaleIDs.size() == 8) continue;

            for(std::vector<std::string> scales: possibleScales){
               if(lines.at(iLine).find(scales[0]) !=  std::string::npos and lines.at(iLine).find(scales[1]) !=  std::string::npos){
                    //Split line in strings by empty space
                    std::vector<std::string> splittedLine = this->SplitString(lines.at(iLine), " ");

                    //Search for 'id="some number"' and extract number
                    for(const std::string& lineInfo: splittedLine){
                        std::size_t idStart = lineInfo.find("id=\"");

                        if(idStart != std::string::npos){
                            std::string partialNumber = lineInfo.substr(idStart+4);
                            std::string number = partialNumber.substr(0, partialNumber.find("\""));
                            scaleIDs.push_back(number);
                        }
                    }
                } 
            }
        }
    } 
}

void PDFWeights::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    std::unique_ptr<std::vector<float>> pdfVariations(new std::vector<float>);
    std::unique_ptr<std::vector<float>> scaleVariations(new std::vector<float>);

    edm::Handle<LHEEventProduct> LHEInfo;
    iEvent.getByToken(LHEToken, LHEInfo);

    for(const gen::WeightsInfo& weight: LHEInfo->weights()){
        if(std::find(variationIDs.begin(), variationIDs.end(), weight.id) != variationIDs.end()){
            pdfVariations->push_back(weight.wgt/LHEInfo->originalXWGTUP());
        }

        if(std::find(scaleIDs.begin(), scaleIDs.end(), weight.id) != variationIDs.end()){
            scaleVariations->push_back(weight.wgt/LHEInfo->originalXWGTUP());
        }
    }
    
    iEvent.put(PDFToken, std::move(pdfVariations));
    iEvent.put(ScaleToken, std::move(scaleVariations));
}

void PDFWeights::endRun(edm::Run const&, edm::EventSetup const&) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PDFWeights);
