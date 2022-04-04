import sys
import FWCore.ParameterSet.Config as cms

process = cms.Process("runQCDRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

process.source = cms.Source("EmptySource")

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.generator = cms.EDFilter("Pythia8GeneratorFilter",
        comEnergy = cms.double(2760.0),
        crossSection = cms.untracked.double(253500000),
        filterEfficiency = cms.untracked.double(1),
        maxEventsToPrint = cms.untracked.int32(0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        PythiaParameters = cms.PSet(
            processParameters = cms.vstring('Main:timesAllowErrors = 10000',
            'ParticleDecays:limitTau0 = on',
            'ParticleDecays:tauMax = 10',
            'HardQCD:all = on',
            'PhaseSpace:pTHatMin = 15',
            'PhaseSpace:pTHatMax = 1000',
            'Tune:pp 14',
            'Tune:ee 7'),
    parameterSets = cms.vstring('processParameters')
	)
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring(
    'CMS_2021_I1963239' # Dijet cross sections and ratios

)
process.rivetAnalyzer.OutputFile      = "QCD-Pt-15to1000-2p76TeV-Pythia8CUETP8M1.yoda"
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = cms.double(253500000)

process.p = cms.Path(process.generator*process.rivetAnalyzer)
