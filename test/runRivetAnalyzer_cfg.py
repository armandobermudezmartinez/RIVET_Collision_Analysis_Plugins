import FWCore.ParameterSet.Config as cms

process = cms.Process("runRivetAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator"),
)
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_TOP_12_028', 'CMS_TOP_12_028_Parton',)

#process.p = cms.Path(process.generator*process.rivetAnalyzer)
process.p = cms.Path(process.rivetAnalyzer)

process.source.fileNames = [
    '/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM/v2/00000/6478E706-FF29-E211-A8EE-001A92811714.root',
    '/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM/v2/00000/6ACD603E-F329-E211-A450-003048679164.root',
    '/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM/v2/00000/6E4F77D7-F929-E211-A049-001A92971B90.root',
    '/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM/v2/00000/BAD0CC54-FF29-E211-A1FF-00261894396D.root',
    '/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM/v2/00000/DE03BB7E-F429-E211-A0B4-001A928116CC.root',
]
