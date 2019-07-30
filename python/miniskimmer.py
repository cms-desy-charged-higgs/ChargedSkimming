import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.AlCa.GlobalTag import GlobalTag
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

import yaml
import os

##Argument parsing
options = VarParsing()
options.register("channel", ["mu4j", "e4j", "mu2j1f", "e2j1f", "mu2f", "e2f"], VarParsing.multiplicity.list, VarParsing.varType.string,
"Name of file for skimming")
options.register("filename", "", VarParsing.multiplicity.list, VarParsing.varType.string,
"Name of file for skimming")
options.register("outname", "outputSkim.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "Name of file for output")
options.register("outdir", "{}/src".format(os.environ["CMSSW_BASE"]), VarParsing.multiplicity.singleton, VarParsing.varType.string, "Dir of file for output")

options.parseArguments()

##Get xSec
xSecFile = yaml.load(file("{}/src/ChargedHiggs/Skimming/data/xsec.yaml".format(os.environ["CMSSW_BASE"]), "r"))

xSec = 1.

for key in xSecFile.keys():
    if key in options.outname:
        xSec = xSecFile[key]["xsec"]

##Check if file is true data file
isData = True in [name in options.outname for name in ["Electron", "Muon", "MET"]]

process = cms.Process("MiniSkimming")
##Load necessary modules
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

##Global tag
tag = "94X_dataRun2_v11" if isData else "94X_mc2017_realistic_v17"
process.GlobalTag = GlobalTag(process.GlobalTag, tag, '')

##Input file
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(["file:{}".format(f) for f in options.filename]))

##Calculate deep flavour discriminator
updateJetCollection(
    process,
    postfix = "RAW",
    jetSource = cms.InputTag('slimmedJets'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring([]), 'None'),
    btagDiscriminators = [
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
    ],
)

updateJetCollection(
    process,
    postfix = 'AK8RAW',
    jetSource = cms.InputTag('slimmedJetsAK8'),
    jetCorrections = ('AK8PFchs', cms.vstring([]), 'None')
)

##https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
setupEgammaPostRecoSeq(process, era='2017-Nov17ReReco')

##Mini Skimmer class which does the skimming
process.skimmer = cms.EDAnalyzer("MiniSkimmer", 
                                jets = cms.InputTag("selectedUpdatedPatJetsRAW"),
                                fatjets = cms.InputTag("updatedPatJetsAK8RAW"),
                                genjets = cms.InputTag("slimmedGenJets"),
                                genfatjets = cms.InputTag("slimmedGenJetsAK8"),
                                mets = cms.InputTag("slimmedMETs"),
                                electrons = cms.InputTag("slimmedElectrons"), 
                                muons = cms.InputTag("slimmedMuons"),
                                trigger = cms.InputTag("TriggerResults","","HLT"),
                                triggerObjects = cms.InputTag("slimmedPatTrigger"),
                                pileUp = cms.InputTag("slimmedAddPileupInfo"),
                                genInfo = cms.InputTag("generator"),
                                genPart = cms.InputTag("prunedGenParticles"),
                                rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                vtx = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                channels = cms.vstring(options.channel[0]),
                                xSec = cms.double(xSec),
                                outFile = cms.string(options.outname),
                                isData = cms.bool(isData),
                )

##Let it run baby
process.p = cms.Path(
                     process.patJetCorrFactorsRAW*process.updatedPatJetsRAW*
                     process.patJetCorrFactorsTransientCorrectedRAW*
                     process.pfImpactParameterTagInfosRAW*
                     process.pfInclusiveSecondaryVertexFinderTagInfosRAW*
                     process.pfDeepCSVTagInfosRAW*process.pfDeepFlavourTagInfosRAW*
                     process.pfDeepFlavourJetTagsRAW*
                     process.updatedPatJetsTransientCorrectedRAW*process.selectedUpdatedPatJetsRAW*
                     process.patJetCorrFactorsAK8RAW*process.updatedPatJetsAK8RAW*
                     process.egammaPostRecoSeq*
                     process.skimmer
            )
