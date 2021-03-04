import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys

options = VarParsing.VarParsing ('standard')
options.register('pdf', 'off', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "PDF:lepton")
options.register('isr', 'off', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "SpaceShower:QEDshowerByL")
options.register('photos', 'off', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "Use Photos")
options.setDefault('maxEvents', 100)
if(hasattr(sys, "argv")):
    options.parseArguments()
print(options)

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = options.maxEvents/100


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# Input source
process.source = cms.Source("EmptySource")

# Input source

process.options = cms.untracked.PSet(

)

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    comEnergy = cms.double(91.2),
    ElectronPositronInitialState = cms.untracked.bool(True),
    crossSection = cms.untracked.double(1),
    filterEfficiency = cms.untracked.double(1),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring(
            'Beams:idA =  11',
            'Beams:idB = -11',
            'Beams:eCM =  91.2',
            # Pythia 8 settings for LEP
            # Hadronic decays including b quarks, with ISR photons switched off
            'WeakSingleBoson:ffbar2gmZ = on',
            '23:onMode = off',
            '23:onIfAny = 13',
            'PDF:lepton = ' + options.pdf,
            'SpaceShower:QEDshowerByL = ' + options.isr,
            # LEP : Make particles with c*tau > 100 mm stable:
            'ParticleDecays:limitTau0 = On',
            'ParticleDecays:tau0Max = 100.0',
        ),
        parameterSets = cms.vstring(
        'processParameters',
        )
    )
)

if options.photos == 'on':
    process.generator.PythiaParameters.processParameters.append(
        'TimeShower:QEDshowerByL = off',
    )
    process.generator.ExternalDecays = cms.PSet(
        Photospp = cms.untracked.PSet(
            parameterSets = cms.vstring("setExponentiation", "setInfraredCutOff", "setMeCorrectionWtForW", "setMeCorrectionWtForZ", "setMomentumConservationThreshold", "setPairEmission", "setPhotonEmission", "setStopAtCriticalError", "suppressAll", "forceBremForDecay"),
            setExponentiation = cms.bool(True),
            setMeCorrectionWtForW = cms.bool(True),
            setMeCorrectionWtForZ = cms.bool(True),
            setInfraredCutOff = cms.double(0.00011),
            setMomentumConservationThreshold = cms.double(0.1),
            setPairEmission = cms.bool(True),
            setPhotonEmission = cms.bool(True),
            setStopAtCriticalError = cms.bool(False),
            # Use Photos only for W/Z decays
            suppressAll = cms.bool(True),
            forceBremForDecay = cms.PSet(
                parameterSets = cms.vstring("Z", "Wp", "Wm"),
                Z = cms.vint32(0, 23),
                Wp = cms.vint32(0, 24),
                Wm = cms.vint32(0, -24),
            ),
        ),
        parameterSets = cms.vstring("Photospp")
    )

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.AnalysisNames = cms.vstring('DELPHI_1994_I375963')
process.rivetAnalyzer.OutputFile = cms.string('test.yoda')

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
#process.lhe_step = cms.Path(process.externalLHEProducer)
process.generation_step = cms.Path(process.generator*process.rivetAnalyzer)

# Schedule definition
#process.schedule = cms.Schedule(process.lhe_step,process.generation_step,process.endjob_step)
process.schedule = cms.Schedule(process.generation_step)

