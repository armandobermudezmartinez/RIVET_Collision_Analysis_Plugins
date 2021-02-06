import sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("runMinBiasRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

process.source = cms.Source("EmptySource")

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
        comEnergy = cms.double(13000.0),
        crossSection = cms.untracked.double(78400000000),
        filterEfficiency = cms.untracked.double(1),
        maxEventsToPrint = cms.untracked.int32(0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        PythiaParameters = cms.PSet(
            processParameters = cms.vstring(
                'SoftQCD:singleDiffractive = on',
        'SoftQCD:doubleDiffractive = on',
        'SoftQCD:centralDiffractive = on',
        'Tune:pp 14',
        'Tune:ee 7',
        'SigmaDiffractive:PomFlux=4',
        'SigmaDiffractive:PomFluxEpsilon=0.1195',
        'SigmaDiffractive:PomFluxAlphaPrime=0.1417',
        'MultipartonInteractions:ecmRef=13000.0',
        'MultipartonInteractions:ecmPow=0.25208',
        'SpaceShower:alphaSvalue=0.1108',
        'SigmaTotal:zeroAXB=off',
        'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
        'MultipartonInteractions:bProfile=2',
        'MultipartonInteractions:pT0Ref=2.481',
        'MultipartonInteractions:coreFraction=0.5457',
        'MultipartonInteractions:coreRadius=0.7195',
        'ColourReconnection:range=2.921',
            ),
    parameterSets = cms.vstring('processParameters')
	)
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring(
    'CMS_2018_I1653948'
)
process.rivetAnalyzer.OutputFile      = "CMS_2018_I1653948_MinBias_CUETP8M4_Pythia8_13TeV.yoda"
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = cms.double(78400000000)

process.p = cms.Path(process.generator*process.rivetAnalyzer)
