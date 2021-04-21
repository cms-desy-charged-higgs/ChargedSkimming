import FWCore.ParameterSet.Config as cms
from FWCore.PythonUtilities.LumiList import LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.AlCa.GlobalTag import GlobalTag
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.MXNet.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsProbs, _pfDeepBoostedJetTagsMetaDiscrs
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

import yaml
import os

##Argument parsing
options = VarParsing()
options.register("channel", "MuonIncl", VarParsing.multiplicity.list, VarParsing.varType.string,
"Channel names")
options.register("era", 2017, VarParsing.multiplicity.singleton, VarParsing.varType.int,
"Channel names")
options.register("filename", "", VarParsing.multiplicity.list, VarParsing.varType.string,
"Name of file for skimming")
options.register("outname", "outputSkim.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "Name of file for output")

options.parseArguments()

##Get xSec
xSecFile = yaml.load(file("{}/src/ChargedSkimming/Skimming/data/xsec.yaml".format(os.environ["CMSSW_BASE"]), "r"))

xSec = 1.

for key in xSecFile.keys():
    if key in options.outname:
        xSec = xSecFile[key]

##Check if file is true data file
isData = True in [name in options.outname for name in ["Electron", "Muon", "MET", "JetHT", "EGamma"]]

process = cms.Process("MiniSkimming")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

##Load necessary modules
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

##Global tag https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
tags = {2016: "102X_mcRun2_asymptotic_v8", 2017: "102X_mc2017_realistic_v8", 2018: "102X_upgrade2018_realistic_v21"}
dataTag = "102X_dataRun2_v13" if not "Run2018D" in options.outname else "102X_dataRun2_Prompt_v16"

tag = dataTag if isData else tags[options.era]
process.GlobalTag = GlobalTag(process.GlobalTag, tag, '')

##Input file
goodLumiFile = "{}/src/ChargedSkimming/Skimming/data/goldenJSON/{}/goodLumi.txt".format(os.environ["CMSSW_BASE"], options.era)

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(["file:{}".format(f) for f in options.filename]),
                            lumisToProcess = LumiList(filename = goodLumiFile).getVLuminosityBlockRange() if isData else cms.untracked.VLuminosityBlockRange()
)

"""process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100),
)"""

##Calculate deep flavour discriminator
updateJetCollection(
    process,
    postfix = "WithDeepB",
    jetSource = cms.InputTag('slimmedJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
    printWarning = False
)

process.goodPatJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("slimmedJets"),
                                        filter = cms.bool(True)
)

##Deep AK8
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsAK8'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    rParam = 0.8,
    jetCorrections = ('AK8PFchs', cms.vstring([]), 'None'),
    btagDiscriminators = _pfDeepBoostedJetTagsProbs + _pfDeepBoostedJetTagsMetaDiscrs,
    postfix='AK8WithDeepTags',
    printWarning = False
)

##Prefiring weight https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
tags = {2016: "2016BtoH", 2017: "2017BtoF", 2018: ""}
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
                            DataEra = cms.string(tags[options.era]),
                            UseJetEMPt = cms.bool(False),
                            PrefiringRateSystematicUncty = cms.double(0.2),
                            SkipWarnings = False
)

##https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
tags = {2016: "2016-Legacy", 2017: "2017-Nov17ReReco", 2018: "2018-Prompt"}
setupEgammaPostRecoSeq(process, era= tags[options.era])

##Producer to get pdf weights
process.pdfweights = cms.EDProducer("PDFWeights",
                                    LHE = cms.InputTag("externalLHEProducer"),
                                    LHAID=cms.int32(306000),
                                    isData = cms.bool(isData),
)

##Mini Skimmer class which does the skimming
process.skimmer = cms.EDAnalyzer("MiniSkimmer", 
                                jets = cms.InputTag("selectedUpdatedPatJetsWithDeepB"),
                                fatjets = cms.InputTag("selectedUpdatedPatJetsAK8WithDeepTags"),
                                genjets = cms.InputTag("slimmedGenJets"),
                                genfatjets = cms.InputTag("slimmedGenJetsAK8"),
                                mets = cms.InputTag("slimmedMETs"),
                                electrons = cms.InputTag("slimmedElectrons"), 
                                muons = cms.InputTag("slimmedMuons"),
                                isoTrack = cms.InputTag("isolatedTracks"), 
                                trigger = cms.InputTag("TriggerResults", "", "HLT"),
                                pileUp = cms.InputTag("slimmedAddPileupInfo"),
                                genInfo = cms.InputTag("generator"),
                                genPart = cms.InputTag("prunedGenParticles"),
                                rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                svtx = cms.InputTag("slimmedSecondaryVertices"),
                                pdf = cms.InputTag("pdfweights","pdfVariations"),
                                scale = cms.InputTag("pdfweights", "scaleVariations"),
                                lhe = cms.InputTag("externalLHEProducer"),
                                channels = cms.vstring(options.channel),
                                xSec = cms.double(xSec),
                                era = cms.int32(options.era),
                                outFile = cms.string(options.outname),
                                isData = cms.bool(isData)
)

##Let it run baby
process.p = cms.Path(
                     process.goodPatJetsPFlow*

                     process.patJetCorrFactorsWithDeepB * 
                     process.updatedPatJetsWithDeepB *
                     process.selectedUpdatedPatJetsWithDeepB *

                     process.patJetCorrFactorsAK8WithDeepTags *
                     process.updatedPatJetsAK8WithDeepTags *
                     process.patJetCorrFactorsTransientCorrectedAK8WithDeepTags *
                     process.pfDeepBoostedJetTagInfosAK8WithDeepTags *
                     process.pfDeepBoostedJetTagsAK8WithDeepTags *
                     process.pfDeepBoostedDiscriminatorsJetTagsAK8WithDeepTags *
                     process.updatedPatJetsTransientCorrectedAK8WithDeepTags *
                     process.selectedUpdatedPatJetsAK8WithDeepTags *

                     process.egammaPostRecoSeq*
                     (process.prefiringweight if options.era != 2018 else cms.Sequence())*
                     process.pdfweights*
                     process.skimmer
)
