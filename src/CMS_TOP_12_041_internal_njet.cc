#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/PartonTop.hh"

namespace Rivet {

class CMS_TOP_12_041_internal_njet : public Analysis {
public:
  /// Minimal constructor
  CMS_TOP_12_041_internal_njet() : Analysis("CMS_TOP_12_041_internal_njet") {}


  /// @name Analysis methods
  //@{

  /// Set up projections and book histograms
  void init() {
    // Parton level top quarks
    addProjection(PartonTop(), "partonTop");

    //FastJets fj(CMSGenParticle(), FastJets::ANTIKT, 0.5);
    VetoedFinalState vfs;
    vfs.addDecayProductsVeto( 24);
    vfs.addDecayProductsVeto(-24);
    FastJets fj(vfs, FastJets::ANTIKT, 0.5);
    fj.useInvisibles();
    addProjection(fj, "Jets");
    
    // Book histograms
    _hVis_nJet30 = bookHisto1D("h02_x01");
    _hVis_nJet60 = bookHisto1D("h02_x02");
    _hVis_nJet100 = bookHisto1D("h02_x03");

    _hVis_addJet1Pt  = bookHisto1D("h03_x01");
    _hVis_addJet1Eta = bookHisto1D("h03_x02");
    _hVis_addJet2Pt  = bookHisto1D("h04_x01");
    _hVis_addJet2Eta = bookHisto1D("h04_x02");
    _hVis_addJJDR   = bookHisto1D("h05_x02");
    _hVis_addJJHT   = bookHisto1D("h05_x03");

  }


  void analyze(const Event& event) {
    const double weight = event.weight();
    // The objects used in the PAPER 12-041 is defined as follows (see p.16 for details):
    // 
    //   * Leptons    : from the W boson decays after FSR
    //   * Jets       : anti-kT R=0.5 to all stable particles 
    //                               exclude W->enu, munu, taunu
    //   * B jet      : B-Ghost matched
    //   * B from top : B hadron from top->b decay
    //
    // Visible phase space definition:
    //
    //   * Leptons         : pT > 20, |eta| < 2.4
    //   * B jets from top : pT > 30, |eta| < 2.4
    //     Additional jets : pT > 20, |eta| < 2.4
    //   * 
    // Full phase space definition:
    //
    //   * Correction to dilepton BR from W boson BR
    //   * No cut on top decay products
    //   * Additional jets : pT > 20, |eta| < 2.4
    //

    const PartonTop& partonTop = applyProjection<PartonTop>(event, "partonTop");
    // Do the analysis only for the ttbar full leptonic channel, removing tau decays
    if ( partonTop.mode() != PartonTop::CH_FULLLEPTON ) vetoEvent;
    if ( partonTop.mode1() >= PartonTop::CH_TAU_HADRON or
         partonTop.mode2() >= PartonTop::CH_TAU_HADRON ) vetoEvent;

    // Apply acceptance cut on muon or electrons
    // Get the lepton of the ttbar decay
    const Particle lep1 = getLast(partonTop.lepton1());
    const Particle lep2 = getLast(partonTop.lepton2());
    if ( lep1.pT() <= 1e-9 or lep2.pT() <= 1e-9 ) vetoEvent; // Just for sanity check

    const Jets& jets = applyProjection<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);

    int nJet30 = 0, nJet60 = 0, nJet100 = 0;
    Jets topBJets, addJets, addBJets;
    foreach ( const Jet& jet, jets ) {
      //if ( deltaR(lep1.momentum(), jet.momentum()) < 0.4 ) continue;
      //if ( deltaR(lep2.momentum(), jet.momentum()) < 0.4 ) continue;

      const double pt = jet.pT();
      if ( pt >  30*GeV ) ++nJet30;
      if ( pt >  60*GeV ) ++nJet60;
      if ( pt > 100*GeV ) ++nJet100;

      bool isBtagged = false, isBFromTop = false;
      foreach(const Particle& p, jet.bTags()) {
        isBtagged = true;
        if ( isFromTop(p) ) {
          isBFromTop = true;
          break;
        }
      }
      if ( isBFromTop ) {
        if ( jet.pT() > 30*GeV ) topBJets.push_back(jet);
      } else {
        addJets.push_back(jet);
        if ( isBtagged ) addBJets.push_back(jet);
      }
    }
    const bool isVisiblePS = lep1.pT() > 20*GeV and std::abs(lep1.eta()) < 2.4 and
                             lep2.pT() > 20*GeV and std::abs(lep2.eta()) < 2.4 and
                             topBJets.size() >= 2;
    if ( isVisiblePS ) {
      fillWithOF(_hVis_nJet30, nJet30, weight);
      fillWithOF(_hVis_nJet60, nJet60, weight);
      fillWithOF(_hVis_nJet100, nJet100, weight);
    }

    // Plots with at least two additional jets
    do {
      if ( addJets.size() < 1 ) break;
      const double ht = std::accumulate(addJets.begin(), addJets.end(),
                                        0., [](double x, const Jet& jj){return x+jj.pT();});
      if ( isVisiblePS ) _hVis_addJJHT->fill(ht, weight);

      const double j1pt = addJets[0].pT(), j1aeta = std::abs(addJets[0].eta());
      if ( isVisiblePS ) {
        _hVis_addJet1Pt ->fill(j1pt  , weight);
        _hVis_addJet1Eta->fill(j1aeta, weight);
      }

      if ( addJets.size() < 2 ) break;

      const double j2pt = addJets[1].pT(), j2aeta = std::abs(addJets[1].eta());
      const double jjdR = deltaR(addJets[0], addJets[1]);

      if ( isVisiblePS ) {
        _hVis_addJet2Pt ->fill(j2pt  , weight);
        _hVis_addJet2Eta->fill(j2aeta, weight);
        _hVis_addJJDR->fill(jjdR, weight);
      }
    } while ( false );

  }

  void finalize() {
    normalize(_hVis_nJet30);
    normalize(_hVis_nJet60);
    normalize(_hVis_nJet100);
    normalize(_hVis_addJet1Pt);
    normalize(_hVis_addJet1Eta );
    normalize(_hVis_addJet2Pt);
    normalize(_hVis_addJet2Eta);
    normalize(_hVis_addJJDR);
    normalize(_hVis_addJJHT);
  }

  //@}


private:
  inline void fillWithOF(Histo1DPtr h, const double x, const double w) const {
    h->fill(std::min(x, h->xMax()-1e-9), w);
  }

  const Particle getLast(const Particle& p) const {
    if ( !p.genParticle() or !p.genParticle()->end_vertex() ) return p;

    foreach (const Particle& dp, Rivet::particles(p.genParticle()->end_vertex(), HepMC::children)) {
      if ( dp.pdgId() == p.pdgId() ) return getLast(dp);;
    }

    return p;
  }

  bool isFromTop(const Particle& p) const {
    if ( !p.genParticle() or !p.genParticle()->production_vertex() ) return false;

    foreach ( const Particle& ap, Rivet::particles(p.genParticle()->production_vertex(), HepMC::ancestors)) {
      if ( abs(ap.pdgId()) == 6 ) return true;
    }

    return false;
  }

  // @name Histogram data members
  //@{
  
  Histo1DPtr _hVis_nJet30, _hVis_nJet60, _hVis_nJet100;
  Histo1DPtr _hVis_addJet1Pt, _hVis_addJet1Eta, _hVis_addJet2Pt, _hVis_addJet2Eta;
  Histo1DPtr _hVis_addJJDR, _hVis_addJJHT;

  //@}

};



// The hook for the plugin system
DECLARE_RIVET_PLUGIN(CMS_TOP_12_041_internal_njet);

}

