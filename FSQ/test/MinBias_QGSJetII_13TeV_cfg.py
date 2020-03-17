import sys
import FWCore.ParameterSet.Config as cms




process = cms.Process("runMinBiasRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

process.source = cms.Source("EmptySource")

#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
from GeneratorInterface.ReggeGribovPartonMCInterface.ReggeGribovPartonMC_AdvancedParameters_cfi import *
process.generator = cms.EDFilter("ReggeGribovPartonMCGeneratorFilter",
                    ReggeGribovPartonMCAdvancedParameters,
                    beammomentum = cms.double(6500),
                    targetmomentum = cms.double(-6500),
                    beamid = cms.int32(1),
                    targetid = cms.int32(1),
                    model = cms.int32(7)
                    )


process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring(
    'CMS_2017_I1511284', # CASTOR energy spectra
    'CMS_2019_I1747892', # CASTOR mulitiplicity-dependent energy spectra
)
process.rivetAnalyzer.OutputFile      = "MinBiasCRMCQGSJetII.yoda"
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = cms.double(78300000000)

process.p = cms.Path(process.generator*process.rivetAnalyzer)
