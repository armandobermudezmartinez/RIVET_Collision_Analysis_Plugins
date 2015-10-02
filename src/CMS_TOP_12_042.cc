#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

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

      double MAXRAPIDITY = 1e10;
      
      vector<pair<PdgId,PdgId> > vidsW;
      vidsW.push_back(make_pair(PID::ELECTRON, PID::NU_EBAR));
      vidsW.push_back(make_pair(PID::POSITRON, PID::NU_E));
      vidsW.push_back(make_pair(PID::MUON,     PID::NU_MUBAR));
      vidsW.push_back(make_pair(PID::ANTIMUON, PID::NU_MU));
      
      vector<pair<PdgId,PdgId> > vidsNuTau;
      vidsNuTau.push_back(make_pair(PID::NU_TAU, PID::NU_TAUBAR));
      
      FinalState fs(-MAXRAPIDITY, MAXRAPIDITY, 0*GeV);
      InvMassFinalState invfsW(fs, vidsW, 75.4*GeV, 85.4*GeV);
      addProjection(invfsW, "INVFSW");
      addProjection(MissingMomentum(), "MET");
      
      VetoedFinalState fsForJets(FinalState(-MAXRAPIDITY, MAXRAPIDITY, 0*GeV));
      fsForJets.addVetoOnThisFinalState(invfsW);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.5), "Jets");

      // Booking of histograms
      _h_wprod_mult = bookHisto1D("wprod_mult", 10, 0, 10);
      
      const std::vector<double> binMet = {0., 27., 52., 87., 130., 172., 300.};
      _h_met = bookHisto1D("met", binMet);
      
      const std::vector<double> binHt = {120., 185., 215., 247., 283., 323., 365., 409., 458., 512., 570., 629., 691., 769., 1000.};
      _h_ht = bookHisto1D("ht", binHt);
      
      const std::vector<double> binSt = {146., 277., 319., 361., 408., 459., 514., 573., 637., 705., 774., 854., 940., 1200.};
      _h_st = bookHisto1D("st", binSt);
      
      const std::vector<double> binPtw = {0., 27., 52., 78., 105., 134., 166., 200., 237., 300.};
      _h_ptw = bookHisto1D("ptw", binPtw);
      
      _h_mwnutau = bookHisto1D("mwnutau", 120, 0, 120);
      _h_nu_tau_pt = bookHisto1D("nu_tau_pt", 120, 0, 120);

    }


    void analyze(const Event& event) {

      double EPSILON = 0.1;

      const double weight = event.weight();
      
      const InvMassFinalState& invMassFinalStateW = applyProjection<InvMassFinalState>(event, "INVFSW");
      const ParticleVector&  WDecayProducts =  invMassFinalStateW.particles();
      _h_wprod_mult->fill(WDecayProducts.size(), weight);
      if (WDecayProducts.size() != 2) vetoEvent; // semi-leptonic ttbar only
      
      // veto taus
      vector<HepMC::GenParticle*> allParticles = particles(event.genEvent());
      for (size_t i = 0; i < allParticles.size(); i++) {
        GenParticle* p = allParticles[i];
        if (p->status() == 1 && abs(p->pdg_id()) == PID::NU_TAU) {
          _h_nu_tau_pt->fill(p->momentum().perp(), weight);
          if (p->momentum().perp() > 10.) vetoEvent;
        }
      }
      
      // TODO: Plot m(nutau, nutaubar), find suitable cut for suppressing W>tau+nutau. Use TauFinder?
      //_h_mwnutau->fill(w.mass()/GeV, weight);
      
      int idxLep = 1;
      if ((fabs(WDecayProducts[1].pdgId()) == PID::NU_MU) || (fabs(WDecayProducts[1].pdgId()) == PID::NU_E)) {
        idxLep = 0;
      }
      
      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _h_met->fill(min(met.visibleMomentum().pT(), 300.-EPSILON)/GeV, weight);
      
      // HT and ST
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);
      
      double ht = 0.0;
      foreach (const Jet& j, jets) {
        if (deltaR(j.momentum(), WDecayProducts[idxLep].momentum()) > 0.3) {
          ht += j.pT();
        }
      }
      double st = ht + WDecayProducts[idxLep].pT() + met.visibleMomentum().pT();
      _h_ht->fill(min(ht, 1000.-EPSILON)/GeV, weight);
      _h_st->fill(min(st, 1200.-EPSILON)/GeV, weight);
      
      // ptW
      FourMomentum w = WDecayProducts[idxLep].momentum() - met.visibleMomentum();
      _h_ptw->fill(min(w.pT(), 300.-EPSILON)/GeV, weight);
    }


    void finalize() {
      //normalize(_h_met);
    }

    //@}


  private:

    // @name Histogram data members
    //@{
    
    Histo1DPtr _h_wprod_mult;
    Histo1DPtr _h_met, _h_ht, _h_st, _h_ptw;
    Histo1DPtr _h_mwnutau, _h_nu_tau_pt;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_TOP_12_042);

}
