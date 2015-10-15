#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/CMSGenParticle.hh"
#include "TopMonteCarlo/RivetTop/interface/PartonTop.hh"

namespace Rivet {

class CMS_TOP_12_041 : public Analysis {
public:
  /// Minimal constructor
  CMS_TOP_12_041() : Analysis("CMS_TOP_12_041") {}


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
    _hVis_addJJMass = bookHisto1D("h05_x01");
    _hVis_addJJDR   = bookHisto1D("h05_x02");
    _hVis_addJJHT   = bookHisto1D("h05_x03");

    _hFull_addJet1Pt  = bookHisto1D("h06_x01");
    _hFull_addJet1Eta = bookHisto1D("h06_x02");
    _hFull_addJet2Pt  = bookHisto1D("h07_x01");
    _hFull_addJet2Eta = bookHisto1D("h07_x02");
    _hFull_addJJMass = bookHisto1D("h08_x01");
    _hFull_addJJDR   = bookHisto1D("h08_x02");
    _hFull_addJJHT   = bookHisto1D("h08_x03");

    _hVis_addBJet1Pt  = bookHisto1D("h09_x01");
    _hVis_addBJet1Eta = bookHisto1D("h09_x02");
    _hVis_addBJet2Pt  = bookHisto1D("h09_x03");
    _hVis_addBJet2Eta = bookHisto1D("h09_x04");
    _hVis_addBBMass = bookHisto1D("h10_x01");
    _hVis_addBBDR   = bookHisto1D("h10_x02");

    _hFull_addBJet1Pt  = bookHisto1D("h11_x01");
    _hFull_addBJet1Eta = bookHisto1D("h11_x02");
    _hFull_addBJet2Pt  = bookHisto1D("h11_x03");
    _hFull_addBJet2Eta = bookHisto1D("h11_x04");
    _hFull_addBBMass = bookHisto1D("h12_x01");
    _hFull_addBBDR   = bookHisto1D("h12_x02");

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
      if ( deltaR(lep1.momentum(), jet.momentum()) < 0.4 ) continue;
      if ( deltaR(lep2.momentum(), jet.momentum()) < 0.4 ) continue;

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
      fillWithOF(_hFull_addJJHT, ht, weight);
      if ( isVisiblePS ) fillWithOF(_hVis_addJJHT, ht, weight);

      const double j1pt = addJets[0].pT(), j1aeta = std::abs(addJets[0].eta());
      fillWithOF(_hFull_addJet1Pt , j1pt  , weight);
      fillWithOF(_hFull_addJet1Eta, j1aeta, weight);
      if ( isVisiblePS ) {
        fillWithOF(_hVis_addJet1Pt , j1pt  , weight);
        fillWithOF(_hVis_addJet1Eta, j1aeta, weight);
      }

      if ( addJets.size() < 2 ) break;

      const double j2pt = addJets[1].pT(), j2aeta = std::abs(addJets[1].eta());
      const double jjmass = (addJets[0].momentum()+addJets[1].momentum()).mass();
      const double jjdR = deltaR(addJets[0], addJets[1]);

      fillWithOF(_hFull_addJet2Pt , j2pt  , weight);
      fillWithOF(_hFull_addJet2Eta, j2aeta, weight);
      if ( isVisiblePS ) {
        fillWithOF(_hVis_addJet2Pt , j2pt  , weight);
        fillWithOF(_hVis_addJet2Eta, j2aeta, weight);
      }

      fillWithOF(_hFull_addJJMass, jjmass, weight);
      fillWithOF(_hFull_addJJDR, jjdR, weight);
      if ( isVisiblePS ) {
        fillWithOF(_hVis_addJJMass, jjmass, weight);
        fillWithOF(_hVis_addJJDR, jjdR, weight);
      }
    } while ( false );

    // Same set of plots if there are additional b-jets
    do {
      if ( addBJets.size() < 1 ) break;
      const double b1pt = addBJets[0].pT(), b1aeta = std::abs(addBJets[0].eta());
      fillWithOF(_hFull_addBJet1Pt , b1pt  , weight);
      fillWithOF(_hFull_addBJet1Eta, b1aeta, weight);
      if ( isVisiblePS ) {
        fillWithOF(_hVis_addBJet1Pt , b1pt  , weight);
        fillWithOF(_hVis_addBJet1Eta, b1aeta, weight);
      }

      if ( addBJets.size() < 2 ) break;
      const double b2pt = addBJets[1].pT(), b2aeta = std::abs(addBJets[1].eta());
      const double bbmass = (addBJets[0].momentum()+addBJets[1].momentum()).mass();
      const double bbdR = deltaR(addBJets[0], addBJets[1]);

      fillWithOF(_hFull_addBJet2Pt , b2pt  , weight);
      fillWithOF(_hFull_addBJet2Eta, b2aeta, weight);
      if ( isVisiblePS ) {
        fillWithOF(_hVis_addBJet2Pt , b2pt  , weight);
        fillWithOF(_hVis_addBJet2Eta, b2aeta, weight);
      }

      fillWithOF(_hFull_addBBMass, bbmass, weight);
      fillWithOF(_hFull_addBBDR, bbdR, weight);

      // Fill plots in visible phase space
      if ( isVisiblePS ) {
        fillWithOF(_hVis_addBBMass, bbmass, weight);
        fillWithOF(_hVis_addBBDR, bbdR, weight);
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
    normalize(_hVis_addJJMass);
    normalize(_hVis_addJJDR);
    normalize(_hVis_addJJHT);
    normalize(_hFull_addJet1Pt);
    normalize(_hFull_addJet1Eta);
    normalize(_hFull_addJet2Pt);
    normalize(_hFull_addJet2Eta);
    normalize(_hFull_addJJMass);
    normalize(_hFull_addJJDR);
    normalize(_hFull_addJJHT);
    normalize(_hVis_addBJet1Pt);
    normalize(_hVis_addBJet1Eta);
    normalize(_hVis_addBJet2Pt);
    normalize(_hVis_addBJet2Eta);
    normalize(_hVis_addBBMass);
    normalize(_hVis_addBBDR);
    normalize(_hFull_addBJet1Pt);
    normalize(_hFull_addBJet1Eta);
    normalize(_hFull_addBJet2Pt);
    normalize(_hFull_addBJet2Eta);
    normalize(_hFull_addBBMass);
    normalize(_hFull_addBBDR);
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
  Histo1DPtr _hVis_addJJMass, _hVis_addJJDR, _hVis_addJJHT;
  Histo1DPtr _hFull_addJet1Pt, _hFull_addJet1Eta, _hFull_addJet2Pt, _hFull_addJet2Eta;
  Histo1DPtr _hFull_addJJMass, _hFull_addJJDR, _hFull_addJJHT;
  Histo1DPtr _hVis_addBJet1Pt, _hVis_addBJet1Eta, _hVis_addBJet2Pt, _hVis_addBJet2Eta;
  Histo1DPtr _hVis_addBBMass, _hVis_addBBDR;
  Histo1DPtr _hFull_addBJet1Pt, _hFull_addBJet1Eta, _hFull_addBJet2Pt, _hFull_addBJet2Eta;
  Histo1DPtr _hFull_addBBMass, _hFull_addBBDR;

  //@}

};



// The hook for the plugin system
DECLARE_RIVET_PLUGIN(CMS_TOP_12_041);

}

