#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {
  namespace { //< only visible in this compilation unit

    /// @brief Special dressed lepton finder
    ///
    /// Find dressed leptons by clustering all leptons and photons
    class SpecialDressedLeptons : public FinalState {
      public:
        /// The default constructor. May specify cuts
        SpecialDressedLeptons(const FinalState& fs, const Cut& cut)
          : FinalState(cut)
        {
          setName("SpecialDressedLeptons");
          IdentifiedFinalState ifs(fs);
          ifs.acceptIdPair(PID::PHOTON);
          ifs.acceptIdPair(PID::ELECTRON);
          ifs.acceptIdPair(PID::MUON);
          addProjection(ifs, "IFS");
          addProjection(FastJets(ifs, FastJets::ANTIKT, 0.1), "LeptonJets");
        }
        
        /// Retrieve the dressed leptons
        const vector<DressedLepton>& dressedLeptons() const { return _clusteredLeptons; }
        
      private:
        /// Container which stores the clustered lepton objects
        vector<DressedLepton> _clusteredLeptons;
        
      public:
        void project(const Event& e) {
          _theParticles.clear();
          _clusteredLeptons.clear();
          
          vector<DressedLepton> allClusteredLeptons;
          
          const Jets jets = applyProjection<FastJets>(e, "LeptonJets").jetsByPt(5*GeV);
          foreach (const Jet& jet, jets) {
            Particle lepCand;
            for (const Particle& cand : jet.particles()) {
              const int absPdgId = abs(cand.pdgId());
              if (absPdgId == PID::ELECTRON || absPdgId == PID::MUON) {
                if (cand.pt() > lepCand.pt()) lepCand = cand;
              }
            }
            //Central lepton must be the major component
            if ((lepCand.pt() < jet.pt()/2.) || (lepCand.pdgId() == 0)) continue;
            
            DressedLepton lepton = DressedLepton(lepCand);
            
            for (const Particle& cand : jet.particles()) {
              if (cand == lepCand) continue;
              lepton.addPhoton(cand, true);
            }
            allClusteredLeptons.push_back(lepton);
          }
          
          for (const DressedLepton& lepton : allClusteredLeptons) {
            if (accept(lepton)) {
              _clusteredLeptons.push_back(lepton);
              _theParticles.push_back(lepton.constituentLepton());
              _theParticles += lepton.constituentPhotons();
            }
          }
        }
    };
  }

  class CMS_TOP_15_006 : public Analysis {
  public:

    /// Minimal constructor
    CMS_TOP_15_006() : Analysis("CMS_TOP_15_006")
    {
    }
  
  public:

    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);
      
      Cut superLooseLeptonCuts = Cuts::pt > 5*GeV;
      SpecialDressedLeptons dressedleptons(fs, superLooseLeptonCuts);
      addProjection(dressedleptons, "DressedLeptons");
      
      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressedleptons);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.5), "Jets");

      // Booking of histograms
      _normedElectronMuonHisto = bookHisto1D("normedElectronMuonHisto", 7, 3.5, 10.5, "Normalized Differential Cross Section in Lepton+Jets Channel", "Jet Multiplicity", "Normed units");
      _absXSElectronMuonHisto = bookHisto1D("absXSElectronMuonHisto", 7, 3.5, 10.5, "Differential Cross Section in Lepton+Jets Channel", "Jet Multiplicity", "pb");
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets
      const SpecialDressedLeptons& dressedleptons = applyProjection<SpecialDressedLeptons>(event, "DressedLeptons");
      
      std::vector<FourMomentum> selleptons;
      
      foreach (const DressedLepton& dressedlepton, dressedleptons.dressedLeptons()) {
        // select good leptons
        if      (dressedlepton.pt() > 30 && dressedlepton.abseta() < 2.4)
          selleptons.push_back(dressedlepton.momentum());
        // veto loose leptons
        else if (dressedlepton.pt() > 15 && dressedlepton.abseta() < 2.5)
          vetoEvent;
      }
      if (selleptons.size() != 1) vetoEvent;
      
      const FourMomentum lepton = selleptons[0];
      
      // jets
      const FastJets& jets   = applyProjection<FastJets>(event, "Jets");
      const Jets      jets30 = jets.jetsByPt(30*GeV);
      int nJets  = 0;
      int nBJets = 0;
      foreach (const Jet& jet, jets30) {
        if (jet.abseta() > 2.5) continue;
        if (deltaR(jet.momentum(), lepton) < 0.5) continue;
        ++nJets;
        if (jet.bTagged(Cuts::pt > 5*GeV))
          ++nBJets;
      }
      
      if (nJets < 4 ||  nBJets < 2) vetoEvent;
      
      // fill histograms
      _normedElectronMuonHisto->fill(min(nJets, 10), weight);
      _absXSElectronMuonHisto ->fill(min(nJets, 10), weight);
    }


    void finalize() {
      const double ttbarXS = !std::isnan(crossSectionPerEvent()) ? crossSection() : 252.89*picobarn;
      if (std::isnan(crossSectionPerEvent()))
        MSG_INFO("No valid cross-section given, using NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV): " << ttbarXS/picobarn << " pb");
      
      const double xsPerWeight = ttbarXS/picobarn / sumOfWeights();
      scale(_absXSElectronMuonHisto, xsPerWeight);
      
      normalize(_normedElectronMuonHisto);
    }

    //@}


  private:

    // @name Histogram data members
    //@{
    
    Histo1DPtr _normedElectronMuonHisto, _absXSElectronMuonHisto;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_TOP_15_006);

}
