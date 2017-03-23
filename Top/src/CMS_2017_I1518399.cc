// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/Cuts.fhh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"

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

      // Jets collection 
      addProjection(FastJets(fs, FastJets::CAM, 1.2), "JetsCA12");

      //partonic top for decay channel defintion
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicTops");       
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicTops");

      //main histograms 
      _hist_mass        = bookHisto1D("d01-x01-y01");
      _hist_mass_norm   = bookHisto1D("d02-x01-y01");

      //histogram for not merged fraction (just in simulation)
      std::vector<double> binning {140, 170, 200, 240, 290, 350};
      _hist_mass_not_merged = bookHisto1D("mass_not_merged", binning, "", "Leading jet mass [GeV]", "");
     
      //control histograms
      _hist_pt1         = bookHisto1D("pt1", 100, 0., 1000., "", "Leading jet p_{T} [GeV]", "");
      _hist_pt2         = bookHisto1D("pt2", 100, 0., 1000., "", "2nd jet p_{T} [GeV]", "");
      _hist_pt3         = bookHisto1D("pt3", 100, 0., 1000., "", "3rd jet p_{T} [GeV]", "");
      _hist_ptlep       = bookHisto1D("ptlep", 100, 0., 1000., "", "Lepton p_{T} [GeV]", "");

      _hist_eta1        = bookHisto1D("eta1", 60, -3., 3., "", "Leading jet #eta [GeV]", "");
      _hist_eta2        = bookHisto1D("eta2", 60, -3., 3., "", "2nd jet #eta [GeV]", "");
      _hist_eta3        = bookHisto1D("eta3", 60, -3., 3., "", "3rd jet #eta [GeV]", "");
      _hist_etalep      = bookHisto1D("etalep", 60, -3., 3., "", "Lepton #eta [GeV]", "");

      _hist_dR_jet2_lep = bookHisto1D("dR_jet2_lep", 50, 0., 5., "", "Lepton #eta [GeV]", "");

      _hist_m1_m2lep    = bookHisto1D("m1_m2lep", 20, 0., 2., "", "m_{jet1}/m_{jet2+lepton} [GeV]", "");

      _hist_Nlep        = bookHisto1D("NLep", 10, 0., 10., "", "Number of promt leptons (e,mu)", "");
      _hist_NQuarks     = bookHisto1D("NQuarks", 10, 0., 10., "", "Number of light Quark candidates", "");
      _hist_pt_top      = bookHisto1D("pt_top", 20, 0., 1000., "", "hadronic top quark pt", "");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

      //decay mode
      const Particles leptonicTops = apply<PartonicTops>(event, "LeptonicTops").particlesByPt();
      const Particles hadronicTops = apply<PartonicTops>(event, "HadronicTops").particlesByPt();  
      if (leptonicTops.size() != 1 || hadronicTops.size() != 1) vetoEvent;  

      //get the lepton
      const Particle lepTop = leptonicTops.at(0);
      const auto isLeptonfromW = [](const Particle& p){return (isChargedLepton(p) && !fromDecay(p));};
      Particles lepton_cadidates = lepTop.allDescendants(firstParticleWith(isLeptonfromW), false);

      // In some cases there is no lepton from the W decay but only leptons from a decay of a radiated gamma. 
      // These events contain fully hadronic ttbar decays and need to be rejected
      Particles leptons;
      for (unsigned int i = 0; i < lepton_cadidates.size(); ++i) {
	Particle lepton_candidate = lepton_cadidates.at(i);
	if ( hasParentWith(lepton_candidate, Cuts::abspid == 22)) continue;
	leptons.push_back(lepton_candidate);
      }
      _hist_Nlep->fill(leptons.size(), weight);
      if(!leptons.size()) vetoEvent;

      Particle lepton = leptons.back();

      //plot the hadronic top pT
      _hist_pt_top->fill(hadronicTops.at(0).pt(), weight);

      //lepton cuts
      if (lepton.pt() < 45 || fabs(lepton.eta()) > 2.1) vetoEvent; 

      //get the jets
      const PseudoJets& psjetsCA12 = applyProjection<FastJets>(event, "JetsCA12").pseudojetsByPt( 50.0*GeV );

      // subtract the lepton four vector from a jet in case of overlap and clean jets
      PseudoJets cleanedJets;
      const fastjet::PseudoJet& pseudoLepton = lepton.pseudojet();
      for (unsigned int i = 0; i < psjetsCA12.size(); ++i) {
	fastjet::PseudoJet jet = psjetsCA12.at(i);
	if (jet.delta_R(pseudoLepton) < 1.2 ) {
	  jet -= pseudoLepton;
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
      if (cleanedJets.at(1).delta_R(pseudoLepton) > 1.2 ) vetoEvent;
        
      // m(jet1) > m(jet2 +lepton)
      const fastjet::PseudoJet& secondJetLepton = cleanedJets.at(1) + pseudoLepton;
      if (cleanedJets.at(0).m() < secondJetLepton.m()) vetoEvent;
   
      //========merged events==================
      //just used as a check (not used for data comparison)
      //======================================
      bool fully_merged = true;
    
      const Particle hadTop = hadronicTops.at(0);
      const auto isQuarkfromTopDecay = [](const Particle& p){return (p.abspid() < 6 && !fromDecay(p) && (hasParentWith(p, Cuts::abspid == 24) || hasParentWith(p, Cuts::abspid == 6)));};
      Particles quarks = hadTop.allDescendants(firstParticleWith(isQuarkfromTopDecay), false);
      _hist_NQuarks->fill(quarks.size(),weight);

      for (unsigned int i = 0; (i < 3 && i < quarks.size()); ++i) {
	const fastjet::PseudoJet& pseudoQuark = quarks.at(i).pseudojet();
	if (cleanedJets.at(0).delta_R(pseudoQuark) > 1.2) fully_merged = false;
      }
      //=======================================

      // fill histograms
      _hist_mass->fill(cleanedJets.at(0).m(), weight);
      if (!fully_merged) _hist_mass_not_merged->fill(cleanedJets.at(0).m(), weight);
      _hist_mass_norm->fill(cleanedJets.at(0).m(), weight);

      _hist_pt1->fill(cleanedJets.at(0).pt(), weight);
      _hist_pt2->fill(cleanedJets.at(1).pt(), weight);
      if (cleanedJets.size() > 2) _hist_pt3->fill(cleanedJets.at(2).pt(), weight);
      _hist_ptlep->fill(lepton.pt(), weight);

      _hist_eta1->fill(cleanedJets.at(0).eta(), weight);
      _hist_eta2->fill(cleanedJets.at(1).eta(), weight);
      if (cleanedJets.size() > 2) _hist_eta3->fill(cleanedJets.at(2).eta(), weight);
      _hist_etalep->fill(lepton.eta(), weight);

      _hist_dR_jet2_lep->fill(cleanedJets.at(1).delta_R(pseudoLepton), weight);

      _hist_m1_m2lep->fill(cleanedJets.at(0).m() / secondJetLepton.m(), weight); 
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
      bool operator() (fastjet::PseudoJet j1, fastjet::PseudoJet j2) {return (j1.pt() > j2.pt());}
    } higherPt;

    Histo1DPtr _hist_mass, _hist_mass_norm, _hist_mass_not_merged;
    Histo1DPtr _hist_pt1, _hist_pt2, _hist_pt3, _hist_ptlep, _hist_pt_top;
    Histo1DPtr _hist_eta1, _hist_eta2, _hist_eta3, _hist_etalep;
    Histo1DPtr _hist_dR_jet2_lep, _hist_m1_m2lep;
    Histo1DPtr _hist_Nlep, _hist_NQuarks;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1518399);


}


