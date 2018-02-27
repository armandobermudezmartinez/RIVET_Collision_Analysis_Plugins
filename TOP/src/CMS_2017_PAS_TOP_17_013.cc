#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
using namespace fastjet;
using namespace fastjet::contrib;

#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/MatrixDiag.hh"
using Rivet::Matrix;
using Rivet::EigenSystem;

namespace Rivet {
  class CMS_2017_PAS_TOP_17_013 : public Analysis {
  public:

    /// Minimal constructor
    CMS_2017_PAS_TOP_17_013() : Analysis("CMS_2017_PAS_TOP_17_013")
    {
    }
  
  private:

    enum Reconstruction { CHARGED=0, ALL=1 };
    enum Observable { MULT=0, PTDS=1, GA_LHA=2, GA_WIDTH=3, GA_THRUST=4, ECC=5, ZG=6, ZGDR=7, NSD=8, TAU21=9, TAU32=10, TAU43=11, C1_00=12, C1_02=13, C1_05=14, C1_10=15, C1_20=16, C2_00=17, C2_02=18, C2_05=19, C2_10=20, C2_20=21, C3_00=22, C3_02=23, C3_05=24, C3_10=25, C3_20=26, M2_B1=27, N2_B1=28, N3_B1=29, M2_B2=30, N2_B2=31, N3_B2=32 };
    enum Flavor { INCL=0, BOTTOM=1, QUARK=2, GLUON=3 };
  
  public:

    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Cuts
      particle_cut = (Cuts::abseta < 5.0) and (Cuts::pT >  0.*GeV);
      lepton_cut   = (Cuts::abseta < 2.4) and (Cuts::pT > 15.*GeV);
      jet_cut      = (Cuts::abseta < 2.4) and (Cuts::pT > 30.*GeV);
      
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
      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, 
                     lepton_cut, /*cluster*/ true, /*useDecayPhotons*/ true);
      addProjection(dressed_leptons, "DressedLeptons");
      
      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.4,
                             JetAlg::ALL_MUONS, JetAlg::NO_INVISIBLES), "Jets");

      // Booking of histograms
      
      for (int r = 0; r < 2; ++r) { // reconstruction (charged, all)
        for (int o = 0; o < 33; ++o) { // observable
          for (int f = 0; f < 4; ++f) { // flavor
            char buffer [11];
            sprintf(buffer, "d%02d-x%02d-y%02d", r+1, o+1, f+1);
            _h[r][o][f] = bookHisto1D(buffer);
          }
        }
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets
      const std::vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      int nsel_leptons = 0;
      for (const DressedLepton& lepton : leptons) {
        if (lepton.pt() > 26.)
          nsel_leptons += 1;
        else
          vetoEvent; // found veto lepton
      }
      if (nsel_leptons != 1)
        vetoEvent;
      
      const Jets all_jets = applyProjection<FastJets>(event, "Jets").jetsByPt(jet_cut);
      if (all_jets.size() < 4)
        vetoEvent;
      
      // categorize jets
      int nsel_bjets = 0;
      int nsel_wjets = 0;
      Jets jets[4];
      for (const Jet& jet : all_jets) {
        // check for jet-lepton overlap -> do not consider for selection
        if (deltaR(jet, leptons[0]) < 0.4)
          continue;
        
        bool overlap = false;
        bool w_jet   = false;
        for (const Jet& jet2 : all_jets) {
          if (jet.momentum() == jet2.momentum())
            continue;
          // check for jet-jet overlap -> do not consider for analysis
          if (deltaR(jet, jet2) < 0.8)
            overlap = true;
          // check for W candidate
          if (jet.bTagged() or jet2.bTagged())
            continue;
          FourMomentum w_cand = jet.momentum() + jet2.momentum();
          if (abs(w_cand.mass() - 80.4) < 15.)
            w_jet = true;
        }

        // count jets for event selection
        if (jet.bTagged())
          nsel_bjets += 1;
        if (w_jet)
          nsel_wjets += 1;
        
        // jets for analysis
        if (jet.abseta() > 2. or overlap) continue;
        
        jets[INCL].push_back(jet);
        if (jet.bTagged())
          jets[BOTTOM].push_back(jet);
        else if (w_jet)
          jets[QUARK].push_back(jet);
        else
          jets[GLUON].push_back(jet);
      }
      
      if (nsel_bjets != 2)
        vetoEvent;
      if (nsel_wjets < 2)
        vetoEvent;
      
      // substructure analysis
      for (int f = 0; f < 4; ++f) {
        for (const Jet& jet : jets[f]) { //FIXME: this calculates everything 2 times :s
          // apply cuts on constituents
          std::vector<PseudoJet> particles[2];
          for (auto p : jet.particles(Cuts::pT > 1.*GeV)) {
            particles[ALL].push_back( PseudoJet(p.px(), p.py(), p.pz(), p.energy()) );
            if (p.charge3() != 0)
              particles[CHARGED].push_back( PseudoJet(p.px(), p.py(), p.pz(), p.energy()) );
          }
          
          if (particles[CHARGED].size() == 0)
            continue;
          
          // recluster with C/A and anti-kt+WTA
          PseudoJet ca_jet[2];
          JetDefinition ca_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
          ClusterSequence ca_charged(particles[CHARGED], ca_def);
          ClusterSequence ca_all(particles[ALL], ca_def);
          ca_jet[CHARGED] = ca_charged.exclusive_jets(1)[0];
          ca_jet[ALL] = ca_all.exclusive_jets(1)[0];
          
          PseudoJet akwta_jet[2];
          JetDefinition akwta_def(fastjet::antikt_algorithm, fastjet::JetDefinition::max_allowable_R, fastjet::RecombinationScheme::WTA_pt_scheme);
          ClusterSequence akwta_charged(particles[CHARGED], akwta_def);
          ClusterSequence akwta_all(particles[ALL], akwta_def);
          akwta_jet[CHARGED] = akwta_charged.exclusive_jets(1)[0];
          akwta_jet[ALL]     = akwta_all.exclusive_jets(1)[0];
          
          // calculate observables
          for (int r = 0; r < 2; ++r) {
            int mult = akwta_jet[r].constituents().size();
            // generalized angularities
            _h[r][MULT][f]->fill(mult, weight);
            if (mult > 1) {
              _h[r][PTDS][f]->fill(getPtDs(akwta_jet[r]), weight);
              _h[r][GA_LHA][f]->fill(calcGA(0.5, 1., akwta_jet[r]), weight);
              _h[r][GA_WIDTH][f]->fill(calcGA(1., 1., akwta_jet[r]), weight);
              _h[r][GA_THRUST][f]->fill(calcGA(2., 1., akwta_jet[r]), weight);
            }
            // eccentricity
            if (mult > 3) {
              _h[r][ECC][f]->fill(getEcc(akwta_jet[r]), weight);
            }
            // N-subjettiness
            if (mult > 2)
              _h[r][TAU21][f]->fill(getTau(2, 1, ca_jet[r]), weight);
            if (mult > 3)
              _h[r][TAU32][f]->fill(getTau(3, 2, ca_jet[r]), weight);
            if (mult > 4)
              _h[r][TAU43][f]->fill(getTau(4, 3, ca_jet[r]), weight);
            // soft drop
            if (mult > 1) {
              std::vector<double> sd_results = getZg(ca_jet[r]);
              if (sd_results[0] > 0.) {
                _h[r][ZG][f]->fill(sd_results[0], weight);
                _h[r][ZGDR][f]->fill(sd_results[1], weight);
              }
            }
            _h[r][NSD][f]->fill(getNSD(0.007, -1., ca_jet[r]), weight);
            // C-series energy correlation ratios
            if (mult > 1) {
              _h[r][C1_00][f]->fill(getC(1, 0.0, ca_jet[r]), weight);
              _h[r][C1_02][f]->fill(getC(1, 0.2, ca_jet[r]), weight);
              _h[r][C1_05][f]->fill(getC(1, 0.5, ca_jet[r]), weight);
              _h[r][C1_10][f]->fill(getC(1, 1.0, ca_jet[r]), weight);
              _h[r][C1_20][f]->fill(getC(1, 2.0, ca_jet[r]), weight);
            }
            if (mult > 2) {
              _h[r][C2_00][f]->fill(getC(2, 0.0, ca_jet[r]), weight);
              _h[r][C2_02][f]->fill(getC(2, 0.2, ca_jet[r]), weight);
              _h[r][C2_05][f]->fill(getC(2, 0.5, ca_jet[r]), weight);
              _h[r][C2_10][f]->fill(getC(2, 1.0, ca_jet[r]), weight);
              _h[r][C2_20][f]->fill(getC(2, 2.0, ca_jet[r]), weight);
            }
            if (mult > 3) {
              _h[r][C3_00][f]->fill(getC(3, 0.0, ca_jet[r]), weight);
              _h[r][C3_02][f]->fill(getC(3, 0.2, ca_jet[r]), weight);
              _h[r][C3_05][f]->fill(getC(3, 0.5, ca_jet[r]), weight);
              _h[r][C3_10][f]->fill(getC(3, 1.0, ca_jet[r]), weight);
              _h[r][C3_20][f]->fill(getC(3, 2.0, ca_jet[r]), weight);
            }
            // M/N-series energy correlation ratios
            if (mult > 2) {
              _h[r][M2_B1][f]->fill(getM(2, 1., ca_jet[r]), weight);
              _h[r][M2_B2][f]->fill(getM(2, 2., ca_jet[r]), weight);
              _h[r][N2_B1][f]->fill(getN(2, 1., ca_jet[r]), weight);
              _h[r][N2_B2][f]->fill(getN(2, 2., ca_jet[r]), weight);
            }
            if (mult > 3) {
              _h[r][N3_B1][f]->fill(getN(3, 1., ca_jet[r]), weight);
              _h[r][N3_B2][f]->fill(getN(3, 2., ca_jet[r]), weight);
            }
          }
        }
      }
    }


    void finalize() {
      for (int r = 0; r < 2; ++r) { // reconstruction (charged, all)
        for (int o = 0; o < 33; ++o) { // observable
          for (int f = 0; f < 4; ++f) { // flavor
            normalize(_h[r][o][f], 1.0, false);
          }
        }
      }
    }

    //@}


  private:
    
    double deltaR(PseudoJet j1, PseudoJet j2) {
      double deta = j1.eta() - j2.eta();
      double dphi = j1.delta_phi_to(j2);
      return sqrt(deta*deta + dphi*dphi);
    }
    
    double getPtDs(PseudoJet jet) {
      double mult   = jet.constituents().size();
      double sumpt  = 0.; // would be jet.pt() in WTA scheme but better keep it generic
      double sumpt2 = 0.;
      for (auto p : jet.constituents()) {
        sumpt  += p.pt();
        sumpt2 += pow(p.pt(), 2);
      }
      double ptd = sumpt2/pow(sumpt,2);
      return max(0., sqrt((ptd-1./mult) * mult/(mult-1.)));
    }
    
    double calcGA(double beta, double kappa, PseudoJet jet) {
      double sumpt = 0.;
      for (auto p : jet.constituents()) {
        sumpt += p.pt();
      }
      double ga = 0.;
      for (auto p : jet.constituents()) {
        ga += pow(p.pt()/sumpt, kappa) * pow(deltaR(jet, p)/0.4, beta);
      }
      return ga;
    }
    
    double getEcc(PseudoJet jet) {
      // Covariance matrix
      Matrix<2> M;
      foreach (auto p, jet.constituents()) {
        Matrix<2> MPart;
        MPart.set(0, 0, (p.eta() - jet.eta()) * (p.eta() - jet.eta()));
        MPart.set(0, 1, (p.eta() - jet.eta()) * mapAngleMPiToPi(p.phi() - jet.phi()));
        MPart.set(1, 0, mapAngleMPiToPi(p.phi() - jet.phi()) * (p.eta() - jet.eta()));
        MPart.set(1, 1, mapAngleMPiToPi(p.phi() - jet.phi()) * mapAngleMPiToPi(p.phi() - jet.phi()));
        M += MPart * p.e();
      }
      // Calculate eccentricity from eigenvalues
      const EigenSystem<2> eigen = diagonalize(M);
      return 1. - eigen.getEigenValues()[1]/eigen.getEigenValues()[0];
    }
    
    double getTau(int N, int M, PseudoJet jet) {
      NsubjettinessRatio tau_ratio(N, M, OnePass_WTA_KT_Axes(), NormalizedMeasure(1.0, 0.4));
      return tau_ratio(jet);
    }
    
    std::vector<double> getZg(PseudoJet jet) {
      PseudoJet jet0 = jet;
      PseudoJet jet1, jet2;
      double zg = 0.;
      while (zg < 0.1 and jet0.has_parents(jet1, jet2)) {
        zg   = jet2.pt()/jet0.pt();
        jet0 = jet1;
      }
      if (zg < 0.1) return {-1., -1.};
      std::vector<double> results;
      results.push_back(zg);
      results.push_back(jet1.delta_R(jet2));
      return results;
    }
    
    int getNSD(double zcut, double beta, PseudoJet jet) {
      PseudoJet jet0 = jet;
      PseudoJet jet1, jet2;
      int nsd = 0.;
      double zg = 0.;
      while (jet0.has_parents(jet1, jet2)) {
        zg = jet2.pt()/jet0.pt();
        if (zg > zcut * pow(jet1.delta_R(jet2)/0.4, beta))
          nsd += 1;
        jet0 = jet1;
      }
      return nsd;
    }
    
    double getC(int N, double beta, PseudoJet jet) {
      EnergyCorrelatorDoubleRatio C(N, beta);
      return C(jet);
    }
    
    double getM(int N, double beta, PseudoJet jet) {
      EnergyCorrelatorMseries CM(N, beta);
      return CM(jet);
    }
    
    double getN(int N, double beta, PseudoJet jet) {
      EnergyCorrelatorNseries CN(N, beta);
      return CN(jet);
    }

    // @name Histogram data members
    //@{
    
    Cut particle_cut, lepton_cut, jet_cut;
    Histo1DPtr _h[2][33][4];

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_PAS_TOP_17_013);

}
