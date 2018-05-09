// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class CMS_2012_S9383875 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2012_S9383875()
      : Analysis("CMS_2012_S9383875")
    {        

    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      const FinalState cnfs(-4, 4);
      addProjection(cnfs, "FS");
      addProjection(FastJets(cnfs, FastJets::ANTIKT, 0.5), "Jets");
      
      std::vector<std::pair<double, double> > eta_m;
      eta_m.push_back(make_pair(-2.4,2.4));
      
      IdentifiedFinalState mufs(Cuts::abseta < 2.4 && Cuts::pT > 9.0*GeV);
      mufs.acceptId(PID::MUON);
      mufs.acceptId(PID::ANTIMUON);
      addProjection(mufs, "Muons");
      
      /// @todo Book histograms here, e.g.:
      _h_dsigdpty05 = bookHisto1D(4, 1, 1);
      _h_dsigdpty10 = bookHisto1D(5, 1, 1);
      _h_dsigdpty15 = bookHisto1D(6, 1, 1);
      _h_dsigdpty20 = bookHisto1D(7, 1, 1);
      _h_dsigdpty22 = bookHisto1D(8, 1, 1);
      _h_dsigdpt    = bookHisto1D(9, 1, 1);
      _h_dsigdy     = bookHisto1D(11, 1, 1);
       
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      /// @todo Do the event by event analysis here

      const FastJets& fastjets = applyProjection<FastJets>(event, "Jets"); 
      const Jets jets = fastjets.jetsByPt(10.);

      const FinalState& muons      = applyProjection<FinalState>(event, "Muons");

      bool onebtag=false;

      foreach (const Jet& j, jets) {
	
	const double ptB= j.momentum().pT();
	const double yB= j.momentum().y();

	bool btag=false;

	foreach (const GenParticle* p, particles(event.genEvent())) {
	  
	  const PdgId pid = p->pdg_id();
	  if (abs(pid) == 5) { 
	    double difference=deltaR(j.momentum().eta(),j.momentum().phi(),p->momentum().eta(),p->momentum().phi());
	    if(sqrt(difference)<0.3){
	      btag=true;
	      onebtag=true;
	    }
	  }
	}
	
	if(btag){
	
       	  if( fabs(yB) < 0.5) { _h_dsigdpty05->fill( ptB, weight );}
	  else if( fabs(yB) > 0.5 && fabs(yB) < 1.0) { _h_dsigdpty10->fill( ptB, weight );}
	  else if( fabs(yB) > 1.0 && fabs(yB) < 1.5) { _h_dsigdpty15->fill( ptB, weight );}
	  else if( fabs(yB) > 1.5 && fabs(yB) < 2.0) { _h_dsigdpty20->fill( ptB, weight );}
	  else if( fabs(yB) > 2.0 && fabs(yB) < 2.2) { _h_dsigdpty22->fill( ptB, weight );}
	}

	if(ptB>30 && onebtag){

	  if(muons.size()>1){	  
	    _h_dsigdpt->fill( ptB, weight );
	    _h_dsigdy->fill( fabs(yB), weight );	
	  }
	}
      }
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      /// @todo Normalise, scale and otherwise manipulate histograms here
      
      double invlumi = crossSection()/picobarn/sumOfWeights();

      scale(_h_dsigdpty05, invlumi); 
      scale(_h_dsigdpty10, invlumi); 
      scale(_h_dsigdpty15, invlumi); 
      scale(_h_dsigdpty20, invlumi); 
      scale(_h_dsigdpty22, invlumi/0.4); 

      double invlumiNano = crossSection()/nanobarn/sumOfWeights();

      scale(_h_dsigdpt, invlumiNano); 
      scale(_h_dsigdy, invlumiNano); 

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_dsigdpty05;
    Histo1DPtr _h_dsigdpty10;
    Histo1DPtr _h_dsigdpty15;
    Histo1DPtr _h_dsigdpty20;
    Histo1DPtr _h_dsigdpty22;
    Histo1DPtr _h_dsigdpt;
    Histo1DPtr _h_dsigdy;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_S9383875);

}
