#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  class CMS_TOP_15_006_2 : public Analysis {
  public:

    /// Minimal constructor
    CMS_TOP_15_006_2() : Analysis("CMS_TOP_15_006_2")
    {
    }
  
  public:

    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for dressed electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      IdentifiedFinalState leptons(fs);
      leptons.acceptIdPair(PID::ELECTRON);
      leptons.acceptIdPair(PID::MUON);
      addProjection(leptons, "Leptons");
      Cut looseLeptonCuts = Cuts::abseta < 2.5 && Cuts::pt > 15*GeV;
      //Cut superLooseLeptonCuts = Cuts::pt > 5*GeV;
      DressedLeptons dressedleptons(photons, leptons, 0.1, looseLeptonCuts, true, false);
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
      ////std::cout << "analyze()" << std::endl;
      const double weight = event.weight();
      
      // select ttbar -> lepton+jets
      const DressedLeptons& dressedleptons = applyProjection<DressedLeptons>(event, "DressedLeptons");
      
      //std::cout << "-- Markus --" << std::endl;
      //for(unsigned int i=0 ; i<dressedleptons.dressedLeptons().size() ; i++) cout<<"lepton pT: "<<dressedleptons.dressedLeptons()[i].momentum().pT()<<" eta: "<<dressedleptons.dressedLeptons()[i].momentum().eta()<<" id: "<<dressedleptons.dressedLeptons()[i].pdgId()<<" const lepton pT: "<<dressedleptons.dressedLeptons()[i].constituentLepton().pT()<<endl;
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
      const double ttbarXS = !isnan(crossSectionPerEvent()) ? crossSection() : 252.89*picobarn;
      if (isnan(crossSectionPerEvent()))
        MSG_INFO("No valid cross-section given, using NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV): " << ttbarXS/picobarn << " pb");
      
      normalize(_normedElectronMuonHisto);
      
      const double xsPerWeight = ttbarXS/picobarn / sumOfWeights();
      scale(_absXSElectronMuonHisto, xsPerWeight);
    }

    //@}


  private:

    // @name Histogram data members
    //@{
    
    Histo1DPtr _normedElectronMuonHisto, _absXSElectronMuonHisto;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_TOP_15_006_2);

}
