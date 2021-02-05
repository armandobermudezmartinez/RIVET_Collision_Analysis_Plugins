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
                'Main:timesAllowErrors    = 10000',
                'ParticleDecays:limitTau0 = on',
                'ParticleDecays:tauMax = 10',
                'SoftQCD:nonDiffractive = on',
                'SoftQCD:singleDiffractive = on',
                'SoftQCD:centralDiffractive = on',
                'SoftQCD:doubleDiffractive = on',
                'Tune:pp 14',
                'Tune:ee 7',
                'MultipartonInteractions:pT0Ref=2.4024',
                'MultipartonInteractions:ecmPow=0.25208',
                'MultipartonInteractions:expPow=1.6',
            ),
    parameterSets = cms.vstring('processParameters')
	)
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring(
    'CMS_2018_I1653948'
)
process.rivetAnalyzer.OutputFile      = "CMS_2018_I1653948_MinBias_CUETP8M1_Pythia8_13TeV.yoda"
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = cms.double(78400000000)

process.p = cms.Path(process.generator*process.rivetAnalyzer)
