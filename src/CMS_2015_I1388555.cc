// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"

#include <vector>
#include <string>


//#include "TopMonteCarlo/RivetTop/interface/PseudoBoostedTop.hh"


namespace Rivet {

  namespace {

    class PseudoBoostedTop : public FinalState {
    public:

      // Default constructor
      PseudoBoostedTop(double lepMinPt = 45., double lepMaxEta = 2.1, double jetMinPt = 30., double jetMaxEta = 2.4, double bjetR = 0.5, double tjetR = 0.8)
        : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
	  _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta), _bjetR(bjetR), _tjetR(tjetR)
      {
	setName("PseudoBoostedTop");
      }
      
      enum DecayMode {CH_HADRON, CH_MUON, CH_ELECTRON, CH_TAU};
      enum TTbarMode {CH_DILEP, CH_ALLHADRON, CH_SEMILEP_MU, CH_SEMILEP_EL, CH_SEMILEP_TAU};
      
      TTbarMode mode() const {
        if ((_topDecay == CH_MUON && _antitopDecay == CH_HADRON) || (_topDecay == CH_HADRON && _antitopDecay == CH_MUON)) return CH_SEMILEP_MU;
        else if ((_topDecay == CH_ELECTRON && _antitopDecay == CH_HADRON) || (_topDecay == CH_HADRON && _antitopDecay == CH_ELECTRON)) return CH_SEMILEP_EL;
        else if ((_topDecay == CH_TAU && _antitopDecay == CH_HADRON) || (_topDecay == CH_HADRON && _antitopDecay == CH_TAU)) return CH_SEMILEP_TAU;
        else if (_topDecay == CH_HADRON && _antitopDecay == CH_HADRON) return CH_ALLHADRON;
        else return CH_DILEP;
      }
      
      DecayMode topDecay() const { return _topDecay;}
      DecayMode antitopDecay() const {return _antitopDecay;}
      
      /// Clone on the heap.
      virtual const Projection* clone() const {
        return new PseudoBoostedTop(*this);
      }
      
      //@}
      
    public:
      Particle partonTop() const   {
        if (_antitopDecay == CH_HADRON) return _antitop;
        else if (_topDecay == CH_HADRON) return _top;
        else return NULL;
      }
      float mtt() const {return _mtt;}
      Jet particleTop() const      {return _particleTop;}
      
      bool isParticleLep() const {return _isParticleLep;}
      bool isGenBJet() const {return _isGenBJet;}
      bool isGenTopJet() const {return _isGenTopJet;}
      bool passParticle() const    {return _passParticle;}
      
    protected:
      // Apply the projection to the event
      void project(const Event& e) override;
      bool isFromTop(const GenParticle* p);
      double getDeltaR(const FourMomentum& p1, const FourMomentum& p2);
      bool fromHardMuon(const GenParticle* p);
      
    private:
      const double _lepMinPt, _lepMaxEta, _jetMinPt, _jetMaxEta;
      const double _bjetR, _tjetR;
      
      constexpr static double _halfpi = 3.1415 / 2.;
      
    private:
      bool _isValid;
      
      DecayMode _topDecay;
      DecayMode _antitopDecay;
      
      Particle _top;
      Particle _antitop;
      
      Jet _particleTop;
      
      bool _isParticleLep;
      bool _isGenBJet;
      bool _isGenTopJet;
      bool _passParticle;
      
      float _mtt;
      
    };
    
  }
  
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
        //cout << "Event " << nEvents << endl;

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

  namespace {
    struct GreaterByPt
    {
      bool operator()(const Particle& a, const Particle& b) {
	return a.pt() > b.pt();
      }
    };
    
    struct GreaterByPtJet
    {
      bool operator()(const Jet& a, const Jet& b) {
	return a.momentum().pt() > b.momentum().pt();
      }
    };
    
    double PseudoBoostedTop::getDeltaR(const FourMomentum& p1, const FourMomentum& p2){
      double deta = std::fabs(p1.eta() - p2.eta());
      double dphi = std::fabs(p1.phi() - p2.phi());
      return std::sqrt(deta * deta + dphi * dphi);
    }
    
    // Determine if lepton comes from hard top decay
    bool PseudoBoostedTop::isFromTop(const GenParticle* p){
      bool fromTop = false;
      const int pdgId = p->pdg_id();
      GenVertex* prodVtx = p->production_vertex();
      if (prodVtx != NULL) {
	foreach (const GenParticle* p1, Rivet::particles(prodVtx, HepMC::parents)) {
	  const int pdgIdParent = p1->pdg_id();
	  if (pdgIdParent == 24) { // parent is a W+ --> leptonic top
	    if (pdgId == -11) _topDecay = CH_ELECTRON;
	    if (pdgId == -13) _topDecay = CH_MUON;
	    if (pdgId == -15) _topDecay = CH_TAU;
	    fromTop = true;
	  }
	  else if (pdgIdParent == -24) { // parent is a W- --> leptonic antitop
	    if (pdgId == 11) _antitopDecay = CH_ELECTRON;
	    if (pdgId == 13) _antitopDecay = CH_MUON;
	    if (pdgId == 15) _antitopDecay = CH_TAU;
	    fromTop = true;
	  }
	}
      }
      return fromTop;  
    }
    
    bool PseudoBoostedTop::fromHardMuon(const GenParticle* p){
      GenVertex* prodVtx = p->production_vertex();
      if (prodVtx == NULL) return false;
      foreach (const GenParticle* ancestor, Rivet::particles(prodVtx, HepMC::ancestors)) {
	if (std::abs(ancestor->pdg_id()) == 13 && ancestor->status() == 3) return true;
      }
      return false;
    }
    
    void PseudoBoostedTop::project(const Event& e) {
      // Lepton : genParticle
      // (b) jets: anti-kt R=0.5 using all final-state particles excluding neutrinos
      // top jets: CA R=0.8 using all particles excluding neutrinos with 140 < mjet < 250
      
      _isValid = false;
      _theParticles.clear();
      _topDecay = CH_HADRON;
      _antitopDecay = CH_HADRON;
      _isParticleLep = false;
      _isGenBJet = false;
      _isGenTopJet = false;
      _passParticle = false;
      
      _mtt = -1.;
      
      // Get parton-level top / state
      Particles pForJet;
      
      Particle refLep;
      
      std::set<int> pFromMuon;
      
      foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
	const int status = p->status();
	const int pdgId = p->pdg_id();
	
	Particle rp(*p);
	
	// Get top quarks
	if (pdgId == 6 && (status == 3 || status == 22)){
	  _top = rp;
	}
	if (pdgId == -6 && (status == 3 || status == 22)){
	  _antitop = rp;
	}
	
	// Get electrons
	if (std::abs(pdgId) == 11){
	  if (isFromTop(p)) {
	    refLep = rp;
	  }      
	}
	
	// Get muons
	if (std::abs(pdgId) == 13){
	  if (isFromTop(p)) {
	    refLep = rp;
	  }      
	}
	
	// Get taus
	if (std::abs(pdgId) == 15){
	  if (isFromTop(p)) {
	    refLep = rp;
	  }      
	}
	
	if (status == 1) {
	  if (fromHardMuon(p)) {
	    pFromMuon.insert(p->barcode());
	  }
	}
      }
      
      const FourMomentum& ttbarP4 = _top.momentum() + _antitop.momentum();
      _mtt = ttbarP4.mass();
      
      // Now start with particle-level quantities
      
      // Get particles for jet clustering
      foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
	const int status = p->status();
	const int barcode = p->barcode();
	Particle rp(*p);
	if (status == 1 && !rp.isNeutrino() && pFromMuon.find(barcode) == pFromMuon.end()) { //consider all final state particles except neutrinos and decay products of hard muon
	  pForJet.push_back(rp);
	} 
      }
      
      // In case the clustering alters the particle collection, clone it
      Particles pForBjet = pForJet;
      Particles pForTjet = pForJet;
      
      // Then do the AK5 clustering
      FastJets ak5Jet(FastJets::ANTIKT, _bjetR);
      ak5Jet.calc(pForBjet);
      
      // Then do the CA8 clustering
      FastJets ca8Jet(FastJets::CAM, _tjetR);
      ca8Jet.calc(pForTjet);
      
      if (refLep.momentum().pt() > _lepMinPt && std::fabs(refLep.eta()) < _lepMaxEta){
	
	_isParticleLep = true;
	const FourMomentum& refLepP4 = refLep.momentum();
	
	Jets genBjets;
	Jets genTjets;
	int nGenBjets = 0;
	int nGenTjets = 0;
	
	foreach (const Jet& jet, ak5Jet.jetsByPt(_jetMinPt)) {
	  if (std::fabs(jet.eta()) > _jetMaxEta) continue;
	  if (getDeltaR(jet.momentum(),refLepP4) > _halfpi) continue;
	  if (getDeltaR(jet.momentum(),refLepP4) < 0.1) continue;
	  genBjets.push_back(jet);
	  nGenBjets += 1;
	}
	
	foreach (const Jet& jet, ca8Jet.jetsByPt(_jetMinPt)) {
	  if (std::fabs(jet.eta()) > _jetMaxEta) continue;
	  if (getDeltaR(jet.momentum(), refLepP4) < _halfpi) continue;
	  if (jet.momentum().mass() < 140.) continue;
	  if (jet.momentum().mass() > 250.) continue;
	  genTjets.push_back(jet);
	  nGenTjets += 1;
	}
	
	if (genTjets.size() > 1) std::sort(genTjets.begin(), genTjets.end(), GreaterByPtJet());
	
	if (nGenBjets >= 1) _isGenBJet = true;
	if (nGenBjets >= 1 && nGenTjets >= 1) _isGenTopJet = true;
	if (_isGenTopJet == true) _particleTop = genTjets[0];
	if (_particleTop.momentum().pT() > 400.) _passParticle = true;
      }
      _isValid = true;
    } 
  }
}
