// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

#include <cstdlib>
#include <cstring>

#include <iostream>

namespace Rivet {

  class CMS_2022_PAS_SMP_21_003 : public Analysis {
  public:

    /// Constructor
    CMS_2022_PAS_SMP_21_003()
      : Analysis("CMS_2022_PAS_SMP_21_003")
    {  }


    /// Book histograms and initialise projections before the run
    void init() {

      // Get options 
      //_mode = 1;
      //if ( getOption("LMODE") == "EL" ) _mode = 0;
      //if ( getOption("LMODE") == "MU" ) _mode = 1;

      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);

      ZFinder zFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, (_mode ? PID::MUON : PID::ELECTRON), 76.0*GeV, 106.0*GeV);
      declare(zFinder, "ZFinder");

      VetoedFinalState jetConstits(visfs);
      jetConstits.addVetoOnThisFinalState(zFinder);

      FastJets jets(jetConstits, FastJets::ANTIKT, 0.4);
      declare(jets, "jets");

      book(_h_deltaPhi_exJMult_Zptbin1,1, 1, 1);
      book(_h_deltaPhi_exJMult_Zptbin2,2, 1, 1);
      book(_h_deltaPhi_exJMult_Zptbin3,3, 1, 1);
      book(_h_deltaPhi_exJMult_Zptbin4,4, 1, 1);
      book(_h_deltaPhi_exJMult_Zptbin5,5, 1, 1);
      book(_h_deltaPhi_ZJ_Zptbin1,6, 1, 1);
      book(_h_deltaPhi_ZJ_Zptbin2,7, 1, 1);
      book(_h_deltaPhi_ZJ_Zptbin3,8, 1, 1);
      book(_h_deltaPhi_ZJ_Zptbin4,9, 1, 1);
      book(_h_deltaPhi_ZJ_Zptbin5,10, 1, 1);
      book(_h_deltaPhi_JJ_Zptbin1,11, 1, 1);
      book(_h_deltaPhi_JJ_Zptbin2,12, 1, 1);
      book(_h_deltaPhi_JJ_Zptbin3,13, 1, 1);
      book(_h_deltaPhi_JJ_Zptbin4,14, 1, 1);
      book(_h_deltaPhi_JJ_Zptbin5,15, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {;

      const ZFinder& zFS = applyProjection<ZFinder>(event, "ZFinder");
      const Particles& zs = zFS.bosons();

      // We did not find exactly one Z. No good.
      if (zs.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      // Find the (dressed!) leptons
      //const Particles& dressedLeptons =  zFS.constituents(cmpMomByPt);
      const Particles& dressedLeptons =  zFS.constituents();
      if (dressedLeptons.size() != 2 || dressedLeptons[0].charge() * dressedLeptons[1].charge() > 0) vetoEvent;

      // leading lepton pt > 25 GeV
      if (dressedLeptons[0].pT() < 25.) vetoEvent;

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = applyProjection<FastJets>(event, "jets");
      const Jets& jets = fj.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);

      // Perform lepton-jet overlap
      Jets goodjets;
      for (const Jet& j : jets) {
        // Decide if this jet is "good", i.e. isolated from the leptons
        bool overlap = false;
        for (const Particle& l : dressedLeptons) {
          if (Rivet::deltaR(j, l) < 0.4) {
            overlap = true;
            break;
          }
        }

        // Fill HT and good-jets collection
        if (overlap) continue;
        goodjets.push_back(j);
      }

      // Weight to be used for histo filling
      const Particle& z = zs[0];

      // Fill jet number histograms
      if (z.pT()<10.) _h_deltaPhi_exJMult_Zptbin1->fill(goodjets.size());
      if (10.<=z.pT() && z.pT()<30.) _h_deltaPhi_exJMult_Zptbin2->fill(goodjets.size());
      if (30.<=z.pT() && z.pT()<50.) _h_deltaPhi_exJMult_Zptbin3->fill(goodjets.size());
      if (50.<=z.pT() && z.pT()<100.) _h_deltaPhi_exJMult_Zptbin4->fill(goodjets.size());
      if (100.<=z.pT()) _h_deltaPhi_exJMult_Zptbin5->fill(goodjets.size());

      // Fill leading jet histograms
      if (goodjets.size() < 1) return;      
      const Jet& j1 = goodjets[0];
      const double dphiZJ = deltaPhi(j1.phi(), z.phi());
      if (z.pT()<10.) _h_deltaPhi_ZJ_Zptbin1->fill(dphiZJ);
      if (10.<=z.pT() && z.pT()<30.) _h_deltaPhi_ZJ_Zptbin2->fill(dphiZJ);
      if (30.<=z.pT() && z.pT()<50.) _h_deltaPhi_ZJ_Zptbin3->fill(dphiZJ);
      if (50.<=z.pT() && z.pT()<100.) _h_deltaPhi_ZJ_Zptbin4->fill(dphiZJ);
      if (100.<=z.pT()) _h_deltaPhi_ZJ_Zptbin5->fill(dphiZJ);      

      // Fill 2nd jet histograms
      if (goodjets.size() < 2) return;
      const Jet& j2 = goodjets[1];
      const double dphiJJ = deltaPhi(j1.phi(), j2.phi());      
      if (z.pT()<10.) _h_deltaPhi_JJ_Zptbin1->fill(dphiJJ);
      if (10.<=z.pT() && z.pT()<30.) _h_deltaPhi_JJ_Zptbin2->fill(dphiJJ);
      if (30.<=z.pT() && z.pT()<50.) _h_deltaPhi_JJ_Zptbin3->fill(dphiJJ);
      if (50.<=z.pT() && z.pT()<100.) _h_deltaPhi_JJ_Zptbin4->fill(dphiJJ);
      if (100.<=z.pT()) _h_deltaPhi_JJ_Zptbin5->fill(dphiJJ);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      scale(_h_deltaPhi_exJMult_Zptbin1, norm);
      scale(_h_deltaPhi_exJMult_Zptbin2, norm);
      scale(_h_deltaPhi_exJMult_Zptbin3, norm);
      scale(_h_deltaPhi_exJMult_Zptbin4, norm);
      scale(_h_deltaPhi_exJMult_Zptbin5, norm);
      scale(_h_deltaPhi_ZJ_Zptbin1, norm);
      scale(_h_deltaPhi_ZJ_Zptbin2, norm);
      scale(_h_deltaPhi_ZJ_Zptbin3, norm);
      scale(_h_deltaPhi_ZJ_Zptbin4, norm);
      scale(_h_deltaPhi_ZJ_Zptbin5, norm);
      scale(_h_deltaPhi_JJ_Zptbin1, norm);
      scale(_h_deltaPhi_JJ_Zptbin2, norm);
      scale(_h_deltaPhi_JJ_Zptbin3, norm);
      scale(_h_deltaPhi_JJ_Zptbin4, norm);
      scale(_h_deltaPhi_JJ_Zptbin5, norm);

    }

  protected:

    size_t _mode;

  private:

    /// @name Histograms

    Histo1DPtr _h_deltaPhi_exJMult_Zptbin1;
    Histo1DPtr _h_deltaPhi_exJMult_Zptbin2;
    Histo1DPtr _h_deltaPhi_exJMult_Zptbin3;
    Histo1DPtr _h_deltaPhi_exJMult_Zptbin4;
    Histo1DPtr _h_deltaPhi_exJMult_Zptbin5;
    Histo1DPtr _h_deltaPhi_ZJ_Zptbin1;
    Histo1DPtr _h_deltaPhi_ZJ_Zptbin2;
    Histo1DPtr _h_deltaPhi_ZJ_Zptbin3;
    Histo1DPtr _h_deltaPhi_ZJ_Zptbin4;
    Histo1DPtr _h_deltaPhi_ZJ_Zptbin5;
    Histo1DPtr _h_deltaPhi_JJ_Zptbin1;
    Histo1DPtr _h_deltaPhi_JJ_Zptbin2;
    Histo1DPtr _h_deltaPhi_JJ_Zptbin3;
    Histo1DPtr _h_deltaPhi_JJ_Zptbin4;
    Histo1DPtr _h_deltaPhi_JJ_Zptbin5;

  };


  DECLARE_RIVET_PLUGIN(CMS_2022_PAS_SMP_21_003);


}
