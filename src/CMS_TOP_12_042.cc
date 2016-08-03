#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

namespace { //< only visible in this compilation unit

  // @brief Parton level top quark finder
  // 
  // Find top quark in the parton level directly tracking particle history.
  // This does not fit with the Rivet philosophy and can be generator dependent,
  // so please use this with your own risks.
  class PartonTop : public FinalState {
  public:
    enum TTbarMode { CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
    enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor.
    PartonTop() : FinalState(-MAXDOUBLE, MAXDOUBLE, 0.0*GeV)
    {
      setName("PartonTop");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new PartonTop(*this);
    }

    //@}

    TTbarMode mode() const {
      const bool isLepton1 = _mode1%3 != 0;
      const bool isLepton2 = _mode2%3 != 0;
      if      (  isLepton1 &&  isLepton2 ) return CH_FULLLEPTON;
      else if ( !isLepton1 && !isLepton2 ) return CH_FULLHADRON;
      return CH_SEMILEPTON;
    }
    DecayMode mode1() const { return _mode1; }
    DecayMode mode2() const { return _mode2; }

    Particle t1() const { return _t1; }
    Particle t2() const { return _t2; }
    Particle b1() const { return _b1; }
    Particle b2() const { return _b2; }
    ParticleVector wDecays1() const { return _wDecays1; }
    ParticleVector wDecays2() const { return _wDecays2; }
    Particle lepton1() const { return findLepton(_wDecays1); };
    Particle lepton2() const { return findLepton(_wDecays2); };

  protected:
    // Apply the projection to the event
    void project(const Event& e) {
      _theParticles.clear();
      _wDecays1.clear();
      _wDecays2.clear();
      _mode1 = _mode2 = CH_HADRON; // Set default decay mode to full-hadronic
      _t1 = _t2 = _b1 = _b2 = Particle();

      const double ptmin = 0;
      const double etamin = -MAXDOUBLE, etamax = MAXDOUBLE;

      int nTop = 0;
      bool isTau1 = false, isTau2 = false;
      foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
        const int pdgId = p->pdg_id();
        const int absId = abs(pdgId);
        if ( absId > 20 ) continue; // We are only interested in quarks and leptons
        //if ( PID::isHadron(pdgId) ) continue; // skip hadrons
        //if ( pdgId == 22 ) continue; // skip photons
        //if ( pdgId == 91 or pdgId == 92 ) continue; // Skip cluster, strings

        if ( isZero(p->momentum().perp()) || p->momentum().perp() < ptmin ) continue;
        if ( !inRange(p->momentum().eta(), etamin, etamax) ) continue;

        // Avoid double counting by skipping if particle ID == parent ID
        std::vector<GenParticle*> pps;
        if ( absId == 6 and p->end_vertex() != 0 ) {
          pps = Rivet::particles(p->end_vertex(), HepMC::children);
        }
        else if ( absId != 6 and p->production_vertex() != 0 )
        {
          pps = Rivet::particles(p->production_vertex(), HepMC::parents);
        }
        else continue;

        bool isDuplicated = false;
        foreach (GenParticle* pp, pps) {
          if ( p != pp && p->pdg_id() == pp->pdg_id() ) {
            isDuplicated = true;
            break;
          }
        }
        if ( isDuplicated ) continue;

        // Build Rivet::Particle
        Particle rp(*p);
        // Skip particles from hadronization (and keep tau decay)
        if ( rp.fromDecay() and !rp.hasAncestor(15) and !rp.hasAncestor(-15) ) continue;

        if      ( pdgId ==  6 ) { nTop++; _t1 = rp; }
        else if ( pdgId == -6 ) { nTop++; _t2 = rp; }
        else if ( pdgId ==  5 and rp.pT() > _b1.pT() ) _b1 = rp;
        else if ( pdgId == -5 and rp.pT() > _b2.pT() ) _b2 = rp;
        else if ( absId <= 16 && rp.hasAncestor( 24) ) {
          if ( pdgId == -15 ) isTau1 = true;
          else if ( pdgId == -11 ) _mode1 = CH_ELECTRON;
          else if ( pdgId == -13 ) _mode1 = CH_MUON;
          _wDecays1.push_back(rp);
        }
        else if ( absId <= 16 && rp.hasAncestor(-24) ) {
          if ( pdgId == 15 ) isTau2 = true;
          else if ( pdgId == 11 ) _mode2 = CH_ELECTRON;
          else if ( pdgId == 13 ) _mode2 = CH_MUON;
          _wDecays2.push_back(rp);
        }
      }
      if ( isTau1 ) _mode1 = static_cast<DecayMode>(_mode1+3);
      if ( isTau2 ) _mode2 = static_cast<DecayMode>(_mode2+3);

    }

    Particle findLepton(const ParticleVector& v) const {
      Particle pp = Particle();
      foreach (const Particle& p, v) {
        const int aid = std::abs(p.pdgId());
        if ( (aid == 11 or aid == 13) and  p.pT() > pp.pT() ) pp = p;
      }
      return pp;
    }

  private:
    DecayMode _mode1, _mode2;
    Particle _t1, _t2;
    Particle _b1, _b2;
    ParticleVector _wDecays1, _wDecays2;
  };

  }

  class CMS_TOP_12_042 : public Analysis {
  public:

    /// Minimal constructor
    CMS_TOP_12_042() : Analysis("CMS_TOP_12_042")
    {
    }


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Parton level top quarks
      addProjection(PartonTop(), "partonTop");

      // Complete final state
      double MAXRAPIDITY = 1e10;
      FinalState fs(-MAXRAPIDITY, MAXRAPIDITY, 0*GeV);
      
      // Projection for taus
      TauFinder taus(TauFinder::ANY);
      addProjection(taus, "Tau");
      IdentifiedFinalState nu_taus(fs);
      nu_taus.acceptIdPair(PID::NU_TAU);
      addProjection(nu_taus, "nu_tau");

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressedelectrons(photons, electrons, 0.1, Cuts::open(), true, false);
      addProjection(dressedelectrons, "DressedElectrons");
      
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressedmuons(photons, muons, 0.1, Cuts::open(), true, false);
      addProjection(dressedmuons, "DressedMuons");
      
      // Projection for jets
      VetoedFinalState fsForJets(FinalState(-MAXRAPIDITY, MAXRAPIDITY, 0*GeV));
      fsForJets.addVetoOnThisFinalState(dressedmuons);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.5), "Jets");
      
      // Projections for MET
      addProjection(MissingMomentum(), "MET");
      
      // Weight counter
      _vis_unit_weights = 0.;

      // Booking of histograms      
      const std::vector<double> binMet = {0., 27., 52., 87., 130., 172., 300.};
      _h_met = bookHisto1D("met", binMet);
      
      const std::vector<double> binHt = {120., 185., 215., 247., 283., 323., 365., 409., 458., 512., 570., 629., 691., 769., 1000.};
      _h_ht = bookHisto1D("ht", binHt);
      
      const std::vector<double> binSt = {146., 277., 319., 361., 408., 459., 514., 573., 637., 705., 774., 854., 940., 1200.};
      _h_st = bookHisto1D("st", binSt);
      
      const std::vector<double> binPtw = {0., 27., 52., 78., 105., 134., 166., 200., 237., 300.};
      _h_ptw = bookHisto1D("ptw", binPtw);
      
      _h_mode0 = bookHisto1D("mode0", 10, 0, 10);
      _h_mode1 = bookHisto1D("mode1", 10, 0, 10);
      _h_mode2 = bookHisto1D("mode2", 10, 0, 10);

    }


    void analyze(const Event& event) {
      double EPSILON = 0.1;
      
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets without taus
      const DressedLeptons& dressedelectrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressedmuons = applyProjection<DressedLeptons>(event, "DressedMuons");
      if (dressedelectrons.dressedLeptons().size() + dressedmuons.dressedLeptons().size() != 1) vetoEvent;
      
      FourMomentum lepton;
      if (dressedelectrons.dressedLeptons().size() == 1) lepton = dressedelectrons.dressedLeptons()[0].momentum();
      else lepton = dressedmuons.dressedLeptons()[0].momentum();
      
      const TauFinder& taus = applyProjection<TauFinder>(event, "Tau");
      const IdentifiedFinalState nu_taus = applyProjection<IdentifiedFinalState>(event, "nu_tau");
      foreach(const Particle& tau, taus.taus()) {
        foreach(const Particle& nu, nu_taus.particles()) {
          if (tau.pid() * nu.pid() < 0) continue;
          const FourMomentum wCandidate = tau.momentum() + nu.momentum();
          if (abs(wCandidate.mass() - 80.4) > 5.) vetoEvent;
        }
      }
      
      // plot parton level decay mode (for debugging)
      const PartonTop& partonTop = applyProjection<PartonTop>(event, "partonTop");
      _h_mode0->fill(partonTop.mode(), weight);
      _h_mode1->fill(partonTop.mode1(), weight);
      _h_mode2->fill(partonTop.mode2(), weight);

      // count weights in visible phase space
      if (weight != 0.) _vis_unit_weights += weight/std::abs(weight);
      
      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _h_met->fill(min(met.visibleMomentum().pT(), 300.-EPSILON)/GeV, weight);
      
      // HT and ST
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);
      
      double ht = 0.0;
      foreach (const Jet& j, jets) {
        if (deltaR(j.momentum(), lepton) > 0.3) {
          ht += j.pT();
        }
      }
      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _h_ht->fill(min(ht, 1000.-EPSILON)/GeV, weight);
      _h_st->fill(min(st, 1200.-EPSILON)/GeV, weight);
      
      // ptW
      FourMomentum w = lepton - met.visibleMomentum();
      _h_ptw->fill(min(w.pT(), 300.-EPSILON)/GeV, weight);
    }


    void finalize() {
      const double s = 1./_vis_unit_weights;
      scale(_h_met, s);
      scale(_h_ht, s);
      scale(_h_st, s);
      scale(_h_ptw, s);
    }

    //@}


  private:

    // @name Histogram data members
    //@{
    
    double _vis_unit_weights;
    
    Histo1DPtr _h_met, _h_ht, _h_st, _h_ptw;
    Histo1DPtr _h_mode0, _h_mode1, _h_mode2;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_TOP_12_042);

}
