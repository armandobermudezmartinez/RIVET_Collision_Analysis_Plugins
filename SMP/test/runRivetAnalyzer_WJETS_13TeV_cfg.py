import FWCore.ParameterSet.Config as cms

process = cms.Process("runRivetAnalysis")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

									  
									  #'file:gen.root'
									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0022F532-C2CE-E111-B33F-003048D476A6.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00929F05-20CF-E111-AC61-003048D45FEC.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00935A20-A9CE-E111-AA3C-001E67398223.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0099FCD0-BECE-E111-AE0E-003048673F2C.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00A1C606-70CF-E111-A038-003048D47770.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00DEB323-6FCF-E111-B7C4-001E67396E64.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00E56CEE-BCCE-E111-920B-00304867407E.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/020463F0-0CCF-E111-A032-002481E14F1E.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0218D322-6FCF-E111-95F3-003048D4601E.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/022A564C-24CF-E111-B6C8-0025B3E05DDA.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/024AF908-18CF-E111-825F-003048D4610E.root'


)
)

#process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")

process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2017_I1491953')

process.rivetAnalyzer.OutputFile= cms.string('wjets.yoda')

process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")

process.p = cms.Path(process.generator*process.rivetAnalyzer)

#process.p = cms.Path(process.rivetAnalyzer)


