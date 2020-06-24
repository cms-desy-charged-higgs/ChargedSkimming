import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.AlCa.GlobalTag import GlobalTag
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer

import yaml
import os

##Argument parsing
options = VarParsing()
options.register("channel", "MuonIncl", VarParsing.multiplicity.list, VarParsing.varType.string,
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
isData = True in [name in options.outname for name in ["Electron", "Muon", "MET", "JetHT"]]
isSignal = "HPlus" in options.outname

process = cms.Process("MiniSkimming")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

##Load necessary modules
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

##Global tag
tag = "102X_dataRun2_v12" if isData else "102X_mc2017_realistic_v7"
process.GlobalTag = GlobalTag(process.GlobalTag, tag, '')

##Input file
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(["file:{}".format(f) for f in options.filename]))
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

"""
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(1),
)
"""

##Calculate deep flavour discriminator
updateJetCollection(
    process,
    postfix = "WithDeepB",
    jetSource = cms.InputTag('slimmedJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
    btagDiscriminators = [
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
    ],
    printWarning = False
)

##Deep AK8
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsAK8'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    rParam = 0.8,
    jetCorrections = ('AK8PFchs', cms.vstring([]), 'None'),
    btagDiscriminators = [
        "pfDeepBoostedJetTags:probHbb",
        "pfDeepBoostedJetTags:probTbcq",
        "pfDeepBoostedJetTags:probTbqq",
        "pfDeepBoostedJetTags:probTbc",
        "pfDeepBoostedJetTags:probTbq",
        "pfDeepBoostedJetTags:probWcq",
        "pfDeepBoostedJetTags:probWqq", 
        "pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD", 
    ],
    postfix='AK8WithDeepTags',
    printWarning = False
)

##Prefiring weight https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
                            DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
                            UseJetEMPt = cms.bool(False),
                            PrefiringRateSystematicUncty = cms.double(0.2),
                            SkipWarnings = False
)

##https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
setupEgammaPostRecoSeq(process, era='2017-Nov17ReReco')

##Producer to get pdf weights
process.pdfweights = cms.EDProducer("PDFWeights",
                                    LHE = cms.InputTag("externalLHEProducer" if not isSignal else "source"),
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
                                trigger = cms.InputTag("TriggerResults","","HLT"),
                                pileUp = cms.InputTag("slimmedAddPileupInfo"),
                                genInfo = cms.InputTag("generator"),
                                genPart = cms.InputTag("prunedGenParticles"),
                                rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                svtx = cms.InputTag("slimmedSecondaryVertices"),
                                pdf = cms.InputTag("pdfweights","pdfVariations"),
                                scale = cms.InputTag("pdfweights", "scaleVariations"),
                                channels = cms.vstring(options.channel),
                                xSec = cms.double(xSec),
                                outFile = cms.string(options.outname),
                                isData = cms.bool(isData)
)

##Let it run baby
process.p = cms.Path(
                     process.patJetCorrFactorsWithDeepB * 
                     process.updatedPatJetsWithDeepB *
                     process.pfImpactParameterTagInfosWithDeepB *
                     process.pfInclusiveSecondaryVertexFinderTagInfosWithDeepB *
                     process.pfDeepCSVTagInfosWithDeepB * 
                     process.pfDeepFlavourTagInfosWithDeepB * 
                     process.pfDeepFlavourJetTagsWithDeepB *
                     process.patJetCorrFactorsTransientCorrectedWithDeepB * 
                     process.updatedPatJetsTransientCorrectedWithDeepB *
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
                     process.prefiringweight*
                     process.pdfweights*
                     process.skimmer
)
