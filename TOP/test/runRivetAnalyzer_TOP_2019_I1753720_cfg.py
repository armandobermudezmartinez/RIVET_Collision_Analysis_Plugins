import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('standard')
options.register('runOnly', '', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Run only specified analysis")
options.register('yodafile', 'test.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Name of yoda output file")
options.setDefault('maxEvents', 1000000)
if(hasattr(sys, "argv")):
    options.parseArguments()
print(options)

process = cms.Process("runRivetAnalysis")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.generator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator", "", "SIM"),
    signalParticlePdgIds = cms.vint32()
)

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

if options.runOnly:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(options.runOnly)
else:
    process.rivetAnalyzer.AnalysisNames = cms.vstring(
        #'CMS_2016_I1434354', # diff xs lepton+jets
        #'MC_TTBAR', # MC analysis for lepton+jets
        #'MC_TOPMASS_LJETS', # MC analysis for lepton+jets top mass
        #'CMS_LesHouches2015', # MC analysis for dilepton
        #'MC_GENERIC', # MC generic analysis
        #'MC_XS', # MC xs analysis
        'CMS_2019_I1753720',  # all-had ttbb
    )
process.rivetAnalyzer.OutputFile      = options.yodafile
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = cms.double(731) # powheg NLO
# process.rivetAnalyzer.CrossSection    = 13.93 # ttbb 4FS MG5
# process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)

process.p = cms.Path(process.mergedGenParticles*process.generator*process.rivetAnalyzer)

process.source.fileNames = [

# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/FC4581E1-28FE-E611-9B1B-0025901D4898.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/F2E6D27E-FDFD-E611-AF19-00266CF9B878.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/EAF4C7DF-28FE-E611-83CD-3417EBE644C2.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/D83DDB86-19FE-E611-8BB1-68B59972BFD8.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/D639CCFE-0BFE-E611-B50C-C45444922BD1.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/D6268336-03FE-E611-A52E-70106F48BB92.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/CA7A11E5-28FE-E611-A5E7-0CC47A7EEE32.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/B86D89FE-0BFE-E611-8564-008CFAF73190.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/9AC62565-05FE-E611-9BC9-001D09FDD7C5.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/4016DE82-11FE-E611-B133-848F69FD297F.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/3619D1E5-03FE-E611-B7AB-782BCB161FC2.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/2A8951C8-28FE-E611-81F9-70106F4D2364.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/1E11E7AA-07FE-E611-AC1D-00266CF9B878.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/1AB0C8D9-06FE-E611-BF68-C45444922949.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/145D51EC-28FE-E611-87DF-90E2BAD814F4.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/08EBAFDD-06FE-E611-8557-001EC9B215E2.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/90000/04E85F2B-0BFE-E611-9CB8-0CC47A7EEE0E.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/DC04A6DF-1FFF-E611-B2B9-0CC47A706CF0.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/D003B1AF-D5FD-E611-A688-848F69FD4C94.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/B0D43593-D5FD-E611-848F-441EA1733CCC.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/A25B1D66-CFFD-E611-BA0D-00266CF32F18.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/92938F96-D7FD-E611-8703-0CC47A706CD6.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/8A13B7D5-1FFF-E611-9694-1866DA89082B.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/78729F87-9DFE-E611-B8DA-047D7B881D3A.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/6C3FD7DE-1FFF-E611-9479-70106F4A8E9C.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/5C17855F-D0FD-E611-86C8-C45444922BFE.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/545885CB-DAFD-E611-8BCF-001D09FDD7C5.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/50579F44-97FE-E611-9C2F-B499BAABD8C6.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/4E6242B1-A1FE-E611-8486-0CC47A706D18.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/4AA8FE75-DFFD-E611-9C4A-001D09FDD7C5.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/4A7A50C6-CFFD-E611-86F8-7845C4FC3623.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/480FB8B0-D1FD-E611-80A2-441EA1733CCC.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/3ACAD00B-CFFD-E611-8BAD-1CB72C1B2D80.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/130000/2415A0D3-1FFF-E611-9266-008CFA001B18.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/F63C3B5D-02FE-E611-B904-70106F4A9948.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/F0FBBA17-E6FE-E611-871B-848F69FD2961.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/C4742C5A-E1FD-E611-9806-002590800C6C.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/C46BD5BD-CDFE-E611-B419-0CC47AC08BF8.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/C0200CA0-09FE-E611-B88B-842B2B760921.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/A67EA779-FBFD-E611-B2AB-008CFAFBF2BE.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/A4C8BD7A-FBFD-E611-9458-1866DA87A664.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/9A4A3D21-E4FE-E611-AB16-1CC1DE18CB10.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/8A776A16-E6FE-E611-9992-68B59972BF62.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/76753DA9-CDFE-E611-BFAA-008CFAFBDC0E.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/7409ED9E-D2FE-E611-81FD-B499BAABCFA4.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/722FF91E-E6FE-E611-AA8A-002590DE6E86.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/5256CB3C-EBFD-E611-A9B0-008CFA001F78.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/5225A232-EBFD-E611-960F-848F69FD4C94.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/4C8CB60A-F3FD-E611-800E-0CC47A706DC0.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/4C27BF54-E6FE-E611-8629-0242AC110005.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/3A6DA46D-CDFE-E611-BB46-70106F4A9384.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/160AEC11-E6FE-E611-8D7E-1CC1DE18CB10.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/0899C78D-E1FD-E611-B980-E41D2D08DF80.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/FEF31B0F-FBFE-E611-887D-002590DE38C8.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/F271A41A-DDFD-E611-BA02-70106F4A9414.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/E84E10B1-CBFD-E611-BC54-F04DA275BFF2.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/E0ECF5DB-D4FD-E611-A7AC-0CC47A706F42.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/DA9905C6-D9FD-E611-82FD-0CC47A7E6A74.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/D05812AE-D5FD-E611-AED9-848F69FD4FC1.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/D005024A-D4FD-E611-8A9E-002590FD5838.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/CC794253-D4FD-E611-9DD7-E41D2D08DDE0.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/C29DD631-CFFD-E611-B102-1CB72C1B64E6.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/C26891FC-FAFE-E611-BBF8-70106F4D2588.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/B45B8485-D7FD-E611-A644-002590FD5838.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/AA275CC3-D7FD-E611-A08C-70106F4A9414.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/A8FEF1B7-D7FD-E611-9F59-002590DE6E30.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/9E0BBDB8-DBFD-E611-9A67-002590DE6E34.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/9C02DB20-FBFE-E611-991B-0CC47A7EEC70.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/90ECA376-DCFD-E611-A04C-C45444922C46.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/7C6C59E8-D0FD-E611-9AD4-0CC47A7E6A4C.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/760AE738-CEFD-E611-9252-C45444922D6C.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/6EB3545F-D0FD-E611-86DF-008CFAFBEA7E.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/629FCF15-18FF-E611-99EE-008CFAFBEFE8.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/623F48E6-D5FD-E611-B3CA-002590DE6E1E.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/58134ECB-DFFD-E611-A48C-002590DE6E34.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/424F731E-18FF-E611-A4F5-0CC47A7EED28.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/3C70031D-FBFE-E611-8591-001E68864A92.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/3C276135-D4FD-E611-9BC2-3417EBE7062D.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/3A551970-DDFD-E611-B847-002590DE6E52.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/386179DB-DAFD-E611-9EDF-1CB72C1B2EF4.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/382D0FF8-E8FD-E611-97C5-7845C4FC3A0D.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/36389C18-18FF-E611-9954-0CC47A7034D2.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/30E70401-E3FD-E611-880E-002590DE3AC0.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/308C7E15-18FF-E611-9607-70106F4D2588.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/110000/0ECB1A20-18FF-E611-ACF2-D8D385FF17CA.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/EA3BB207-C1FD-E611-B52D-F45214938700.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/E070494F-FDFD-E611-9A49-0242AC110004.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/C0FD029D-CAFD-E611-A489-848F69FD4FC1.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/C06476CA-ADFE-E611-A8C1-008CFA000F5C.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/C002BE58-5CFF-E611-8A25-0CC47AC08816.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/BA5F5149-8CFE-E611-9157-1CB72C1B64E2.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/B2DB9C09-8FFE-E611-B64A-002590DE6E78.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/AAAB343B-E4FE-E611-A5F7-1866DA87C2CD.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/A47033C1-C4FD-E611-9B76-0CC47A706FFE.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/A2B7B100-8EFF-E611-8D4D-848F69FD29AF.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/A0B5CDC9-98FF-E611-BD8F-68B59972BF62.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/8C59EDAB-C2FD-E611-BDD3-848F69FD2A30.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/8A6785FC-C3FD-E611-A034-0CC47A7E0106.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/88E16AA0-D1FD-E611-B638-7845C4FC3A0D.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/884C6B23-CFFD-E611-9C4B-70106F4D254C.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/386ED805-03FE-E611-8E44-1866DA8907CB.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/323F5141-06FE-E611-A04B-1866DA8907CB.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/089A5E1C-8CFF-E611-8EEB-0242AC110003.root",
# "/store/mc/RunIISummer16MiniAODv2/ttbb_4FS_ckm_amcatnlo_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/00000/087D1E20-C9FD-E611-AECE-0242AC110005.root",


'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0806AB92-99BE-E611-9ECD-0025905A6138.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/18E31463-B3BE-E611-B6A3-0CC47A4D7678.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22B7C5F1-A4BE-E611-AC0A-0025905B85AE.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/26ABF488-A0BE-E611-BEEB-0CC47A4D7640.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/38AB0251-9FBE-E611-A8C2-0025905A60DE.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/3CCF34DF-9DBE-E611-9512-0025905B858E.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/3E18521B-A4BE-E611-8843-0025905A607E.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4072DE86-A0BE-E611-8C32-0025905A497A.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/445E1780-ABBE-E611-B6D6-0CC47A78A360.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/485B62A7-9CBE-E611-B628-0025905A607A.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4A5DD249-9FBE-E611-8E5F-0CC47A7C345C.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4C7B302E-A1BE-E611-ACD0-0CC47A4D7668.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4E043DF2-A4BE-E611-84AF-0025905B85C0.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4E1AB394-99BE-E611-A6C3-0CC47A4D760A.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/5C4AED1B-A4BE-E611-972C-0CC47A7C3472.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/6225D2FC-9DBE-E611-82E5-0025905A6064.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/62EF3B3B-A9BE-E611-9899-0CC47A7C34E6.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/669298F7-9DBE-E611-9EDD-0CC47A4D768E.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/66A40847-9FBE-E611-82D9-003048FFD71C.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/6C0FE789-99BE-E611-876D-0CC47A7C3430.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/6E6CCEA4-A1BE-E611-BBED-0025905A609E.root',
'/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/6EE6C71F-A4BE-E611-9466-0CC47A4D7640.root',
]
