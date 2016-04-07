// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/PseudoBoostedTop.hh"


namespace Rivet {


  class CMS_2015_I1388555 : public Analysis {
  public:

    /// Constructor
    CMS_2015_I1388555()
      : Analysis("CMS_2015_I1388555")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      addProjection(PseudoBoostedTop(45., 2.1, 30., 2.4, 0.5, 0.8), "ttbar");      
      
      _hEl_topPt_parton        = bookHisto1D("d01-x01-y01"); // dsigma/dpt(top)
      _hEl_topPt_particle      = bookHisto1D("d02-x01-y01"); // dsigma/dpt(top)
      _hEl_topY_parton         = bookHisto1D("d03-x01-y01"); // dsigma/dy(top)
      _hEl_topY_particle       = bookHisto1D("d04-x01-y01"); // dsigma/dy(top)
      
      _hMu_topPt_parton        = bookHisto1D("d05-x01-y01"); // dsigma/dpt(top)
      _hMu_topPt_particle      = bookHisto1D("d06-x01-y01"); // dsigma/dpt(top)
      _hMu_topY_parton         = bookHisto1D("d07-x01-y01"); // dsigma/dy(top)
      _hMu_topY_particle       = bookHisto1D("d08-x01-y01"); // dsigma/dy(top)
      
      _hComb_topPt_parton        = bookHisto1D("d09-x01-y01"); // dsigma/dpt(top)
      _hComb_topPt_particle      = bookHisto1D("d10-x01-y01"); // dsigma/dpt(top)
      _hComb_topY_parton         = bookHisto1D("d11-x01-y01"); // dsigma/dy(top)
      _hComb_topY_particle       = bookHisto1D("d12-x01-y01"); // dsigma/dy(top)

      _hEl_topPt_parton_norm        = bookHisto1D("d13-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hEl_topPt_particle_norm      = bookHisto1D("d14-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hEl_topY_parton_norm         = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dy(top)
      _hEl_topY_particle_norm       = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dy(top)
      
      _hMu_topPt_parton_norm        = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hMu_topPt_particle_norm      = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hMu_topY_parton_norm         = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dy(top)
      _hMu_topY_particle_norm       = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dy(top)

      _hComb_topPt_parton_norm        = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hComb_topPt_particle_norm      = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hComb_topY_parton_norm         = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dy(top)
      _hComb_topY_particle_norm       = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dy(top)

      //Diagnostic histograms (not compared against data)
      _hMu_cutflow = bookHisto1D("mu_cutflow",6,-0.5,5.5);
      _hEl_cutflow = bookHisto1D("el_cutflow",6,-0.5,5.5);
      _hChannel_truth = bookHisto1D("truth_channel",5,-0.5,4.5);
      _hEventsAnalyzed = bookHisto1D("EventsAnalyzed",1,-0.5,0.5);

      nEvents = 0;
      nSemilepTau = 0;
      nAllHadronic = 0;
      nDilep = 0;
      nMu = 0;
      nEl = 0;
      nPassParton_mu = 0;
      nParticleLep_mu = 0;
      nIsBJet_mu = 0;
      nIsTJet_mu = 0;
      nPassParticle_mu = 0;
      nPassParton_el = 0;
      nParticleLep_el = 0;
      nIsBJet_el = 0;
      nIsTJet_el = 0;
      nPassParticle_el = 0;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      nEvents += 1;
      //cout << endl << "Event " << nEvents << endl;
      
      const double weight = event.weight();

      // Get the ttbar candidate
      const PseudoBoostedTop& ttbar = applyProjection<PseudoBoostedTop>(event, "ttbar");

      _hEventsAnalyzed->fill(0.);

      if (ttbar.mode() == PseudoBoostedTop::CH_DILEP) {
	nDilep += 1;
	_hChannel_truth->fill(0.);
      }
      if (ttbar.mode() == PseudoBoostedTop::CH_ALLHADRON) {
	nAllHadronic += 1;
	_hChannel_truth->fill(1.);
      }
      if (ttbar.mode() == PseudoBoostedTop::CH_SEMILEP_TAU) {
	nSemilepTau += 1;
	_hChannel_truth->fill(2.);
      }
      
      if ( ttbar.mode() == PseudoBoostedTop::CH_DILEP || ttbar.mode() == PseudoBoostedTop::CH_ALLHADRON || ttbar.mode() == PseudoBoostedTop::CH_SEMILEP_TAU) vetoEvent;

      // Do selection
      const FourMomentum& partonTopP4 = ttbar.partonTop().momentum();
      
      if ( ttbar.mode() == PseudoBoostedTop::CH_SEMILEP_MU ) { //pass loose parton, muon channel

	nMu += 1;
	_hChannel_truth->fill(3.);
	_hMu_cutflow->fill(0.);
	_hMu_topPt_parton->fill(partonTopP4.pT(), weight);
      	_hMu_topPt_parton_norm->fill(partonTopP4.pT(), weight);
	_hComb_topPt_parton->fill(partonTopP4.pT(), weight);
      	_hComb_topPt_parton_norm->fill(partonTopP4.pT(), weight);

	if (partonTopP4.pT() >= 400.){
	  nPassParton_mu += 1;
	  _hMu_cutflow->fill(1.);
	  _hMu_topY_parton->fill(partonTopP4.rapidity(), weight);
	  _hMu_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
	  _hComb_topY_parton->fill(partonTopP4.rapidity(), weight);
	  _hComb_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
	}

	if (ttbar.isParticleLep()) {
	  nParticleLep_mu += 1;
	  _hMu_cutflow->fill(2.);
	}
	if (ttbar.isGenBJet()) {
	  nIsBJet_mu += 1;
	  _hMu_cutflow->fill(3.);
	}
	
	if (ttbar.isGenTopJet()){
	  nIsTJet_mu += 1;
	  _hMu_cutflow->fill(4.);
	  const FourMomentum& particleTopP4 = ttbar.particleTop().momentum();

	  _hMu_topPt_particle->fill(particleTopP4.pT(), weight);
	  _hMu_topPt_particle_norm->fill(particleTopP4.pT(), weight);
	  _hComb_topPt_particle->fill(particleTopP4.pT(), weight);
	  _hComb_topPt_particle_norm->fill(particleTopP4.pT(), weight);

	  if (ttbar.passParticle()){
	    nPassParticle_mu += 1;
	    _hMu_cutflow->fill(5.);
	    _hMu_topY_particle->fill(particleTopP4.rapidity(), weight);
	    _hMu_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
	    _hComb_topY_particle->fill(particleTopP4.rapidity(), weight);
	    _hComb_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
	  }
	}
      }
      
      else if ( ttbar.mode() == PseudoBoostedTop::CH_SEMILEP_EL ) {
	
	nEl += 1;
	_hChannel_truth->fill(4.);
	_hEl_cutflow->fill(0.);
	_hEl_topPt_parton->fill(partonTopP4.pT(), weight);
      	_hEl_topPt_parton_norm->fill(partonTopP4.pT(), weight);
	_hComb_topPt_parton->fill(partonTopP4.pT(), weight);
      	_hComb_topPt_parton_norm->fill(partonTopP4.pT(), weight);

	if (partonTopP4.pT() >= 400.){
	  nPassParton_el += 1;
	  _hEl_cutflow->fill(1.);
	  _hEl_topY_parton->fill(partonTopP4.rapidity(), weight);
	  _hEl_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
	  _hComb_topY_parton->fill(partonTopP4.rapidity(), weight);
	  _hComb_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
	}

	if (ttbar.isParticleLep()) {
	  nParticleLep_el += 1;
	  _hEl_cutflow->fill(2.);
	}
	if (ttbar.isGenBJet()) {
	  nIsBJet_el += 1;
	  _hEl_cutflow->fill(3.);
	}

	if (ttbar.isGenTopJet()){
	  nIsTJet_el += 1;
	  _hEl_cutflow->fill(4.);
	  const FourMomentum& particleTopP4 = ttbar.particleTop().momentum();
	  
	  _hEl_topPt_particle->fill(particleTopP4.pT(), weight);
	  _hEl_topPt_particle_norm->fill(particleTopP4.pT(), weight);
	  _hComb_topPt_particle->fill(particleTopP4.pT(), weight);
	  _hComb_topPt_particle_norm->fill(particleTopP4.pT(), weight);

	  if (ttbar.passParticle()){
	    nPassParticle_el += 1;
	    _hEl_cutflow->fill(5.);
	    _hEl_topY_particle->fill(particleTopP4.rapidity(), weight);
	    _hEl_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
	    _hComb_topY_particle->fill(particleTopP4.rapidity(), weight);
	    _hComb_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
	  }
	}
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      cout << endl;
      cout << "Number of events analyzed:   " << nEvents << endl;
      cout << "Number of dileptonic events: " << nDilep << endl;
      cout << "Number of hadronic events:   " << nAllHadronic << endl;
      cout << "Number of tau events:        " << nSemilepTau << endl;
      cout << endl;
      cout << "In muon channel," << endl;
      cout << "   pass parton loose:     " << nMu << endl;
      cout << "   pass parton, pt>400:   " << nPassParton_mu << endl;
      cout << "   ==1 particle-level mu: " << nParticleLep_mu << endl;
      cout << "   >=1 gen b jet:         " << nIsBJet_mu << endl;
      cout << "   >=1 gen top jet:       " << nIsTJet_mu << endl;
      cout << "   pass particle, pt>400: " << nPassParticle_mu << endl;
      cout << endl;
      cout << "In electron channel," << endl;
      cout << "   pass parton loose:     " << nEl << endl;
      cout << "   pass parton, pt>400:   " << nPassParton_el << endl;
      cout << "   ==1 particle-level mu: " << nParticleLep_el << endl;
      cout << "   >=1 gen b jet:         " << nIsBJet_el << endl;
      cout << "   >=1 gen top jet:       " << nIsTJet_el << endl;
      cout << "   pass particle, pt>400: " << nPassParticle_el << endl;
      cout << endl;
      
      double xs_mu_parton = 252.89 * 1000. * nPassParton_mu / nMu;
      double xs_mu_particle = 252.89 * 1000. * nPassParticle_mu / nMu;
      double xs_el_parton = 252.89 * 1000. * nPassParton_el / nEl;
      double xs_el_particle = 252.89 * 1000. * nPassParticle_el / nEl;
      double xs_comb_parton = 252.89 * 1000. * (nPassParton_el + nPassParton_mu) / (nEl + nMu);
      double xs_comb_particle = 252.89 * 1000. * (nPassParticle_el + nPassParton_mu) / (nEl + nMu);
      cout << "Inclusive xsec (pt>400), mu channel, parton-level: " << xs_mu_parton << endl;
      cout << "Inclusive xsec (pt>400), mu channel, particle-level: " << xs_mu_particle << endl;
      cout << "Inclusive xsec (pt>400), el channel, parton-level: " << xs_el_parton << endl;
      cout << "Inclusive xsec (pt>400), el channel, particle-level: " << xs_el_particle << endl;
      cout << "Inclusive xsec (pt>400), comb channel, parton-level: " << xs_comb_parton << endl;
      cout << "Inclusive xsec (pt>400), comb channel, particle-level: " << xs_comb_particle << endl;

        
      normalize(_hMu_topPt_parton,xs_mu_parton,false);
      normalize(_hMu_topPt_particle,xs_mu_particle,false);
      normalize(_hMu_topY_parton,xs_mu_parton,false);
      normalize(_hMu_topY_particle,xs_mu_particle,false);
      normalize(_hEl_topPt_parton,xs_el_parton,false);
      normalize(_hEl_topPt_particle,xs_el_particle,false);
      normalize(_hEl_topY_parton,xs_el_parton,false);
      normalize(_hEl_topY_particle,xs_el_particle,false);
      normalize(_hComb_topPt_parton,xs_comb_parton,false);
      normalize(_hComb_topPt_particle,xs_comb_particle,false);
      normalize(_hComb_topY_parton,xs_comb_parton,false);
      normalize(_hComb_topY_particle,xs_comb_particle,false);

      normalize(_hMu_topPt_parton_norm,1.0,false);
      normalize(_hMu_topPt_particle_norm,1.0,false);
      normalize(_hMu_topY_parton_norm,1.0,false);
      normalize(_hMu_topY_particle_norm,1.0,false);
      normalize(_hEl_topPt_parton_norm,1.0,false);
      normalize(_hEl_topPt_particle_norm,1.0,false);
      normalize(_hEl_topY_parton_norm,1.0,false);
      normalize(_hEl_topY_particle_norm,1.0,false); 
      normalize(_hComb_topPt_parton_norm,1.0,false);
      normalize(_hComb_topPt_particle_norm,1.0,false);
      normalize(_hComb_topY_parton_norm,1.0,false);
      normalize(_hComb_topY_particle_norm,1.0,false);

    }

    //@}


  private:

    Histo1DPtr _hMu_topPt_parton;
    Histo1DPtr _hMu_topPt_particle;
    Histo1DPtr _hMu_topY_parton;
    Histo1DPtr _hMu_topY_particle;
    Histo1DPtr _hEl_topPt_parton;
    Histo1DPtr _hEl_topPt_particle;
    Histo1DPtr _hEl_topY_parton;
    Histo1DPtr _hEl_topY_particle;
    Histo1DPtr _hComb_topPt_parton;
    Histo1DPtr _hComb_topPt_particle;
    Histo1DPtr _hComb_topY_parton;
    Histo1DPtr _hComb_topY_particle;
    Histo1DPtr _hMu_topPt_parton_norm;
    Histo1DPtr _hMu_topPt_particle_norm;
    Histo1DPtr _hMu_topY_parton_norm;
    Histo1DPtr _hMu_topY_particle_norm;
    Histo1DPtr _hEl_topPt_parton_norm;
    Histo1DPtr _hEl_topPt_particle_norm;
    Histo1DPtr _hEl_topY_parton_norm;
    Histo1DPtr _hEl_topY_particle_norm;
    Histo1DPtr _hComb_topPt_parton_norm;
    Histo1DPtr _hComb_topPt_particle_norm;
    Histo1DPtr _hComb_topY_parton_norm;
    Histo1DPtr _hComb_topY_particle_norm;
    Histo1DPtr _hMu_cutflow;
    Histo1DPtr _hEl_cutflow;
    Histo1DPtr _hChannel_truth;
    Histo1DPtr _hEventsAnalyzed;

    int nEvents;
    int nDilep;
    int nSemilepTau;
    int nAllHadronic;
    int nMu;
    int nEl;
    int nPassParton_mu;
    int nParticleLep_mu;
    int nIsBJet_mu;
    int nIsTJet_mu;
    int nPassParticle_mu;
    int nPassParton_el;
    int nParticleLep_el;
    int nIsBJet_el;
    int nIsTJet_el;
    int nPassParticle_el;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1388555);


}
