#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {

  // This analysis template can be used to obtain the same physics objects
  // as the ParticleLevelProducer added to GeneratorInterface/RivetInterface.
  // See CMSSW pull requests #18402 (80X) and #18404 (master).
  class ParticleLevelProducerTemplate : public Analysis {
  private:
    bool _usePromptFinalStates;
    bool _excludePromptLeptonsFromJetClustering;
    bool _excludeNeutrinosFromJetClustering;
    
    double _particleMinPt, _particleMaxEta;
    double _lepConeSize, _lepMinPt, _lepMaxEta;
    double _jetConeSize, _jetMinPt, _jetMaxEta;
    double _fatJetConeSize, _fatJetMinPt, _fatJetMaxEta;
    
    std::vector<DressedLepton> _leptons;
    ParticleVector _photons, _neutrinos;
    Jets _jets, _fatjets;
    Vector3 _met;

  public:
    ParticleLevelProducerTemplate() : Analysis("ParticleLevelProducerTemplate"),
    _usePromptFinalStates(true),
    _excludePromptLeptonsFromJetClustering(true),
    _excludeNeutrinosFromJetClustering(true),

    _particleMinPt  (0.),
    _particleMaxEta (5.),
    
    _lepConeSize (0.1),
    _lepMinPt    (15.),
    _lepMaxEta   (2.5),
    
    _jetConeSize (0.4),
    _jetMinPt    (30.),
    _jetMaxEta   (2.4),
    
    _fatJetConeSize (0.8),
    _fatJetMinPt    (200.),
    _fatJetMaxEta   (2.4)
    {
    }

    // Initialize Rivet projections
    void init() {
      // Cuts
      Cut particle_cut = (Cuts::abseta < _particleMaxEta) and (Cuts::pT > _particleMinPt*GeV);
      Cut lepton_cut   = (Cuts::abseta < _lepMaxEta)      and (Cuts::pT > _lepMinPt*GeV);
      
      // Generic final state
      FinalState fs(particle_cut);
      
      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);
      
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);
      
      // useDecayPhotons=true allows for photons with tau ancestor,
      // photons from hadrons are vetoed by the PromptFinalState;
      // will be default DressedLeptons behaviour for Rivet >= 2.5.4
      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, _lepConeSize, 
                     lepton_cut, /*cluster*/ true, /*useDecayPhotons*/ true);
      if (not _usePromptFinalStates)
        dressed_leptons = DressedLeptons(photons, charged_leptons, _lepConeSize, 
                          lepton_cut, /*cluster*/ true, /*useDecayPhotons*/ true);
      addProjection(dressed_leptons, "DressedLeptons");
      
      // Photons
      if (_usePromptFinalStates) {
        // We remove the photons used up for lepton dressing in this case
        VetoedFinalState vetoed_prompt_photons(prompt_photons);
        vetoed_prompt_photons.addVetoOnThisFinalState(dressed_leptons);
        addProjection(vetoed_prompt_photons, "Photons");
      }
      else
        addProjection(photons, "Photons");
      
      // Jets
      VetoedFinalState fsForJets(fs);
      if (_usePromptFinalStates and _excludePromptLeptonsFromJetClustering)
        fsForJets.addVetoOnThisFinalState(dressed_leptons);
      JetAlg::InvisiblesStrategy invisiblesStrategy = JetAlg::DECAY_INVISIBLES;
      if (_excludeNeutrinosFromJetClustering)
        invisiblesStrategy = JetAlg::NO_INVISIBLES;
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, _jetConeSize,
                             JetAlg::ALL_MUONS, invisiblesStrategy), "Jets");
      
      // FatJets
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, _fatJetConeSize), "FatJets");
      
      // Neutrinos
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      if (_usePromptFinalStates) {
        PromptFinalState prompt_neutrinos(neutrinos);
        prompt_neutrinos.acceptMuonDecays(true);
        prompt_neutrinos.acceptTauDecays(true);
        addProjection(prompt_neutrinos, "Neutrinos");
      }
      else
        addProjection(neutrinos, "Neutrinos");
      
      // MET
      addProjection(MissingMomentum(fs), "MET");
    };

    // Apply Rivet projections
    void analyze(const Event& event) {
      _jets.clear();
      _fatjets.clear();
      _leptons.clear();
      _photons.clear();
      _neutrinos.clear();
      
      // Get analysis objects from projections
      Cut jet_cut    = (Cuts::abseta < _jetMaxEta)    and (Cuts::pT > _jetMinPt*GeV);
      Cut fatjet_cut = (Cuts::abseta < _fatJetMaxEta) and (Cuts::pT > _fatJetMinPt*GeV);
      
      _leptons   = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      for (const Particle& lepton : _leptons) std::cout << "lepton pt = " << lepton.pt() << std::endl;
      _jets      = applyProjection<FastJets>(event, "Jets").jetsByPt(jet_cut);
      for (const Jet& jet : _jets) std::cout << "jet pt = " << jet.pt() << std::endl;
      _fatjets   = applyProjection<FastJets>(event, "FatJets").jetsByPt(fatjet_cut);
      _photons   = applyProjection<FinalState>(event, "Photons").particlesByPt();
      _neutrinos = applyProjection<FinalState>(event, "Neutrinos").particlesByPt();
      _met       = -applyProjection<MissingMomentum>(event, "MET").vectorEt();
    };

    // Do nothing here
    void finalize() {};

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ParticleLevelProducerTemplate);
}
