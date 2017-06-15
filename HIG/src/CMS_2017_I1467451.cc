// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetFastJet.hh"
#include "Rivet/Projections/JetAlg.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2017_I1467451 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1467451);

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs(-2.5,2.5,0.0*GeV);
      FinalState fsm(-5,5,0.0*GeV);
      addProjection(fs, "FS");
      addProjection(fsm, "FSM");

      IdentifiedFinalState e(fs,11);
      e.acceptId(-11);
      e.acceptId(13);
      e.acceptId(-13);
      addProjection(e,"E");

      MissingMomentum Met(fsm);
      addProjection(Met, "MET");

      FastJets jets(fsm, FastJets::ANTIKT, 0.5);
      addProjection(jets, "JETS");


      // Book histograms
      histoPtH=bookHisto1D(1,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      Particles e = applyProjection<IdentifiedFinalState>(event, "E").particlesByPt(10.0*GeV);

      if(e.size()<2) vetoEvent;
      if(e[0].momentum().pT()<20*GeV || e[1].momentum().pT()<10*GeV) vetoEvent;
      if(e[0].charge()==e[1].charge()) vetoEvent;
      if(abs(e[0].pdgId())==abs(e[1].pdgId())) vetoEvent;

      FourMomentum LL=(e[0].momentum()+e[1].momentum());

      if(LL.mass()<12*GeV) vetoEvent;
      if(LL.pT()<30*GeV) vetoEvent;

      FourMomentum EtMiss = applyProjection<MissingMomentum>(event,"MET").missingMomentum();
      FourMomentum PtH = LL+EtMiss;

      double phiLL = LL.phi(); 
      double phiEtMiss = EtMiss.phi();
      double phi;
      if (phiLL<phiEtMiss) phi = phiEtMiss-phiLL;
      else phi = phiLL-phiEtMiss;

      double mT = sqrt(2*LL.pT()*EtMiss.pT()*(1-cos(phi)));
      if (mT<50*GeV) vetoEvent;

      const FastJets& jetfs = applyProjection<FastJets>(event, "JETS");
      const Jets& jets = jetfs.jetsByPt(30.*GeV);
      std::vector<Jet> realJets;

      for (size_t i=0; i<jets.size(); ++i) {
        if (deltaR(e[0], jets[i])>0.5 && deltaR(e[1], jets[i])>0.5) {
          realJets.push_back(jets[i]);
        }
      }

      histoPtH->fill(PtH.pT(),weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      histoPtH->normalize(crossSection());

    }


  private:

    Histo1DPtr histoPtH;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1467451);

}
