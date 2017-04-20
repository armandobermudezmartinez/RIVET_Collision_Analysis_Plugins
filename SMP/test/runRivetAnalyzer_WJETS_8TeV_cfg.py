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
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/024AF908-18CF-E111-825F-003048D4610E.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/025DF613-94CE-E111-9A0B-001E67397747.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0266BBD6-B2CE-E111-B574-001E67397C33.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/027CD229-56CE-E111-92EF-001E67396874.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02B3C436-DECE-E111-942F-002590200A88.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0407242D-92CE-E111-A960-001E67397701.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/04202A87-FFCE-E111-8401-001E67398E49.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/043B6F9E-C9CE-E111-BCD4-003048D45FA0.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0466D3F6-AACE-E111-A90E-003048D460DE.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/04A151B9-D9CE-E111-A189-001E67396685.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/04BBF1FC-EACE-E111-AEB6-0030486709B4.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/04C882F1-0CCF-E111-806E-003048D46022.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/04D31037-15CF-E111-90F8-001E67396897.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/04E40FB9-60CE-E111-97E4-001E67397D55.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/061C162E-83CE-E111-8A3C-001E67396BAD.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0654BF4A-A3CE-E111-9FF3-001E67396EAA.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06564B29-0ECF-E111-8AB5-001E67398BE7.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/0656A3BB-65CE-E111-BF18-001E67398CE1.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06608230-9ECE-E111-A968-001E67396653.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/066F610C-A2CE-E111-ADED-001E67397AD0.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/067FC761-08CF-E111-8B95-001E67398412.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06C5A18C-A4CE-E111-B9A0-001E67397238.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06D1E378-D8CE-E111-AA99-001E67396E7D.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06D8F243-75CE-E111-9792-001E67397431.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06DB573A-67CE-E111-8F9C-003048D45FEC.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/06E16028-19CF-E111-87DE-002481E150A2.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/080D7867-08CF-E111-ADAC-003048D476C4.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/081A3F96-1ECF-E111-A66B-001E673966B2.root',
#									  'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/083442F0-7FCE-E111-8057-003048D437CA.root'
									  
)
)

#process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")

process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("genParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2016_I1491953')

process.rivetAnalyzer.OutputFile= cms.string('wjets.yoda')

process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")

process.p = cms.Path(process.generator*process.rivetAnalyzer)

#process.p = cms.Path(process.rivetAnalyzer)


