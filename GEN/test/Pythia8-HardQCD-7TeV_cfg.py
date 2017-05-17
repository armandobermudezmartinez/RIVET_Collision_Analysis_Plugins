# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/MinBias_TuneZ2_7TeV_pythia6_cff.py -s GEN --datatier=GEN-SIM-RAW --conditions auto:mc --eventcontent RAWSIM --no_exec -n 10000 --python_filename=rivet_cfg.py --customise=Configuration/GenProduction/rivet_customize.py
import FWCore.ParameterSet.Config as cms
import os

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200000)
)

process.RandomNumberGeneratorService.generator.initialSeed = int(os.getenv('SEED'))

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('PYTHIA6-MinBias TuneZ2 at 7TeV'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/CMSSW/Configuration/GenProduction/python/Attic/MinBias_TuneZ2_7TeV_pythia6_cff.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MinBias_TuneZ2_7TeV_pythia6_cff_py_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Other statements
process.GlobalTag.globaltag = 'PHYS14_25_V2'

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
				 comEnergy = cms.double(7000.0),
				 crossSection = cms.untracked.double(9.573213e+08),
				 filterEfficiency = cms.untracked.double(1),
				 maxEventsToPrint = cms.untracked.int32(0),
				 pythiaHepMCVerbosity = cms.untracked.bool(False),
				 pythiaPylistVerbosity = cms.untracked.int32(0),
				 PythiaParameters = cms.PSet(
	processParameters = cms.vstring(
	'Main:timesAllowErrors    = 10000',
	'ParticleDecays:limitTau0 = on',
	'ParticleDecays:tauMax = 10',
	'HardQCD:all = on',
	'PhaseSpace:pTHatMin = 45', 
	'Tune:pp 14',
	'Tune:ee 7',
	'MultipartonInteractions:ecmPow=0.25208',
	'SpaceShower:alphaSvalue=0.1108',
        'PDF:pSet=LHAPDF6:NNPDF30_lo_as_0130',
        'ColourReconnection:range=6.59e+00',
	#'MultipartonInteractions:pT0Ref=2.20e+00',
        #'MultipartonInteractions:expPow=1.60e+00',
	#M1
	#'MultipartonInteractions:pT0Ref=2.4024',
        #'MultipartonInteractions:ecmPow=0.25208',
        #'MultipartonInteractions:expPow=1.6',
	#DPS
        'MultipartonInteractions:pT0Ref=2.182',  
        'MultipartonInteractions:expPow=1.208',  	
	),
	parameterSets = cms.vstring('processParameters')
	)
				 )

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.
process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')

# Automatic addition of the customisation function from Configuration.GenProduction.rivet_customize
#from Configuration.GenProduction.rivet_customize import customise 

def customise(process):
	process.load('GeneratorInterface.RivetInterface.rivetAnalyzer_cfi')
	process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_FSQ_12_013')
	process.rivetAnalyzer.CrossSection = cms.double(975700)
	process.rivetAnalyzer.OutputFile = cms.string(os.getenv('OUTPUT'))
	process.generation_step+=process.rivetAnalyzer
	process.schedule.remove(process.RAWSIMoutput_step)
        return(process)      
#call to customisation function customise imported from Configuration.GenProduction.rivet_customize
process = customise(process)

# End of customisation functions
