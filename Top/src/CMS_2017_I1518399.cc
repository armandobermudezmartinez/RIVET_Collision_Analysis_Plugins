// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Tools/Cuts.fhh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {

  class CMS_2017_I1518399 : public Analysis {
  public:

    /// Constructor
    CMS_2017_I1518399()
      : Analysis("CMS_2017_I1518399")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);
      addProjection(fs, "FS");

      //dressed leptons
      IdentifiedFinalState photons(fs);       
      photons.acceptIdPair(PID::PHOTON);             
       
      ChargedLeptons charged_leptons(fs);
      PromptFinalState prompt_leptons(charged_leptons);

      Cut leptonCuts = Cuts::pt > 45*GeV && Cuts::abseta < 2.1;          
   
      DressedLeptons dressed_leptons(photons, prompt_leptons, 0.1, leptonCuts, true, false);       
      declare(dressed_leptons, "DressedLeptons");             

      //jets
      VetoedFinalState fs_jets(FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV));     
      fs_jets.vetoNeutrinos();
      declare(FastJets(fs_jets, FastJets::CAM, 1.2), "JetsCA12");

      //partonic top for decay channel defintion
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicTops");       
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicTops");

      //main histograms 
      _hist_mass        = bookHisto1D("d01-x01-y01");
      _hist_mass_norm   = bookHisto1D("d02-x01-y01");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

      //decay mode
      const Particles leptonicTops = apply<PartonicTops>(event, "LeptonicTops").particlesByPt();
      const Particles hadronicTops = apply<PartonicTops>(event, "HadronicTops").particlesByPt();  
      if (leptonicTops.size() != 1 || hadronicTops.size() != 1) vetoEvent;  

      //get the leptons
      const DressedLeptons& dressed_leptons = apply<DressedLeptons>(event, "DressedLeptons");       

      //leading dressed lepton
      const std::vector<DressedLepton> leptons = dressed_leptons.dressedLeptons();
      if ( leptons.size() == 0 ) vetoEvent;

      Particle lepton;
      double max_lepton_pt = 0.;
      for (unsigned int i = 0; i < leptons.size(); ++i) {
	if (leptons.at(i).pt() > max_lepton_pt) {
	  max_lepton_pt = leptons.at(i).pt();
	  lepton = leptons.at(i);
	}
      }

      //get the jets
      Cut jetCuts = Cuts::pt > 50*GeV;
      const Jets& psjetsCA12 = applyProjection<FastJets>(event, "JetsCA12").jetsByPt( jetCuts );

      // subtract the lepton four vector from a jet in case of overlap and clean jets
      Jets cleanedJets;
  
      for (unsigned int i = 0; i < psjetsCA12.size(); ++i) {
	Jet jet = psjetsCA12.at(i);
	if (deltaR(jet.momentum(), lepton.momentum()) < 1.2 ) {
	  const FourMomentum cleanedMom = jet.momentum() - lepton.momentum();
	  Jet cleanedjet(cleanedMom);
	  jet = cleanedjet;
	}
	if (fabs(jet.eta()) < 2.5) cleanedJets.push_back(jet);
      }

      //sort the cleaned jets
      std::sort(cleanedJets.begin(), cleanedJets.end(), higherPt);
    
      // jet pt cuts
      if (cleanedJets.size() < 2) vetoEvent;
      if (cleanedJets.at(0).pt() < 400) vetoEvent;
      if (cleanedJets.at(1).pt() < 150) vetoEvent;

      // jet veto 
      if (cleanedJets.size() > 2 && cleanedJets.at(2).pt() > 150) vetoEvent;

      // small distance between 2nd jet and lepton
      if (deltaR(cleanedJets.at(1).momentum(), lepton.momentum()) > 1.2 ) vetoEvent;
        
      // m(jet1) > m(jet2 +lepton)
      FourMomentum secondJetLepton = cleanedJets.at(1).momentum() + lepton.momentum();
      if (cleanedJets.at(0).mass() < secondJetLepton.mass()) vetoEvent;

      // fill histograms
      _hist_mass->fill(cleanedJets.at(0).mass(), weight);
      _hist_mass_norm->fill(cleanedJets.at(0).mass(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double sf = 252.89 * 1000 / sumOfWeights();

      scale(_hist_mass, sf);
      normalize(_hist_mass_norm, 1.0, false);
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here

    struct higherPt {
      bool operator() (Jet j1, Jet j2) {return (j1.pt() > j2.pt());}
    } higherPt;

    Histo1DPtr _hist_mass, _hist_mass_norm;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1518399);


}


