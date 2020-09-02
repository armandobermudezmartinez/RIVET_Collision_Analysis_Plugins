import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register(
    'runGridpack', 
    False, 
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, 
    "Run gridpack instead of miniaod file"
)
options.register(
    'runTop', 
    True, 
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, 
    "Run on top quark events (top antiquark events otherwise)"
)
options.setDefault('maxEvents', 1000)
options.parseArguments()

process = cms.Process("runRivetAnalysis")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

if options.runGridpack:
    process.source = cms.Source("EmptySource")
    process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
        nEvents = cms.untracked.uint32(options.maxEvents),
        numberOfParameters = cms.uint32(1),
        outputFile = cms.string('cmsgrid_final.lhe'),
        scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
    )
    
    if options.runTop:
        process.externalLHEProducer.args = cms.vstring(
            '/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/powheg/V2/ST_tch_4f_13TeV_top_powheg_madspin/st_tch_4f_ckm_NLO_top_powheg_madspin_tarball.tar.xz'
        )
    else:
        process.externalLHEProducer.args = cms.vstring(
            '/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/powheg/V2/ST_tch_4f_13TeV_antitop_powheg_madspin/st_tch_4f_ckm_NLO_antitop_powheg_madspin_tarball.tar.xz'
        )

    from Configuration.Generator.Pythia8CommonSettings_cfi import *
    from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
    from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
    
    process.generator = cms.EDFilter("Pythia8HadronizerFilter",
        maxEventsToPrint = cms.untracked.int32(1),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        filterEfficiency = cms.untracked.double(1.0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        comEnergy = cms.double(13000.),
        PythiaParameters = cms.PSet(
            pythia8CommonSettingsBlock,
            pythia8CUEP8M1SettingsBlock,
            pythia8PowhegEmissionVetoSettingsBlock,
            processParameters = cms.vstring(
                'POWHEG:nFinal = 3', ## Number of final state particles
                ## (BEFORE THE DECAYS) in the LHE
                ## other than emitted extra parton
                'TimeShower:mMaxGamma = 1.0',#cutting off lepton-pair production
                ##in the electromagnetic shower
                ##to not overlap with ttZ/gamma* samples
                'Tune:pp 14',
                'Tune:ee 7',
                'MultipartonInteractions:ecmPow=0.25208',
                'SpaceShower:alphaSvalue=0.1108',
                'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
                'MultipartonInteractions:pT0Ref=2.197139e+00',
                'MultipartonInteractions:expPow=1.608572e+00',
                'ColourReconnection:range=6.593269e+00',
            ),
            parameterSets = cms.vstring(
                'pythia8CommonSettings',
                'pythia8PowhegEmissionVetoSettings',
                'processParameters'
            )
        )
    )
else:
    if options.runTop:
        # example file from /ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
        process.source.fileNames = cms.untracked.vstring(
            'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/28B73FDD-B1B9-E611-AFED-002590494C82.root',
        )
    else:
        # example file from /ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM
        process.source.fileNames = cms.untracked.vstring(
            'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/801355FB-ACBD-E611-A24D-00266CF3323C.root',
        )
    
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.genMerger = cms.EDProducer("MergedGenParticleProducer",
            inputPruned = cms.InputTag("prunedGenParticles"),
            inputPacked = cms.InputTag("packedGenParticles")
    )
    process.hepmcProducer = cms.EDProducer("GenParticles2HepMCConverter",
        genParticles = cms.InputTag("genMerger"),
        genEventInfo = cms.InputTag("generator"),
        signalParticlePdgIds = cms.vint32(6,-6)
    )
    
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
if options.runGridpack:
    process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
else:
    process.rivetAnalyzer.HepMCCollection = cms.InputTag("hepmcProducer:unsmeared")
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2019_I1744604')
if options.runTop:
    process.rivetAnalyzer.OutputFile = cms.string('out_top.yoda')
    process.rivetAnalyzer.CrossSection = cms.double(136.02)
else:
    process.rivetAnalyzer.OutputFile = cms.string('out_antitop.yoda')
    process.rivetAnalyzer.CrossSection = cms.double(80.95)
    
if options.runGridpack:
    process.lhe_step = cms.Path(process.externalLHEProducer)
    process.generation_step = cms.Path(process.generator)
    process.analysis_step = cms.Path(process.rivetAnalyzer)
    process.schedule = cms.Schedule(process.lhe_step,process.generation_step,process.analysis_step)
else:
    process.p = cms.Path(process.genMerger*process.hepmcProducer*process.rivetAnalyzer)

