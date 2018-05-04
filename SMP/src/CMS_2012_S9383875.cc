// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

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
       double Bjet18 ;
       double Bjet32 ;
       double Bjetmuon ;

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
       
       IdentifiedFinalState mufs(cnfs);
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
       
       Bjet18 = 0.;
       Bjet32 = 0.;
       Bjetmuon = 0.;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      /// @todo Do the event by event analysis here
      unsigned int num_b_jets = 0;
      const FastJets& fastjets = applyProjection<FastJets>(event, "Jets"); 
      const Jets jets = fastjets.jetsByPt(10.);

      const FinalState& muons      = applyProjection<FinalState>(event, "Muons");

      bool muongoodevent=false;

      bool onebtag=false;

      foreach (const Jet& j, jets) {
	
	bool btag=false;

	foreach (const GenParticle* p, particles(event.genEvent())) {
	  
	  const PdgId pid = p->pdg_id();
	  if (abs(pid) == 5) { 
	    double difference=(j.momentum().eta()-p->momentum().eta())*(j.momentum().eta()-p->momentum().eta())+deltaPhi(j.momentum().phi(),p->momentum().phi())*deltaPhi(j.momentum().phi(),p->momentum().phi());
	    if(sqrt(difference)<0.3){
	      btag=true;
	      onebtag=true;
	    }
	  }
	}
	
	if(btag){
	
	  num_b_jets += 1;
	  const double ptB= j.momentum().pT();
	  const double yB= j.momentum().rapidity();
	  if( fabs(yB) < 0.5) { _h_dsigdpty05->fill( ptB, weight );}
	  else if( fabs(yB) > 0.5 && fabs(yB) < 1.0) { _h_dsigdpty10->fill( ptB, weight );}
	    else if( fabs(yB) > 1.0 && fabs(yB) < 1.5) { _h_dsigdpty15->fill( ptB, weight );}
	    else if( fabs(yB) > 1.5 && fabs(yB) < 2.0) { _h_dsigdpty20->fill( ptB, weight );}
	    else if( fabs(yB) > 2.0 && fabs(yB) < 2.2) { _h_dsigdpty22->fill( ptB, weight );}
	  //          cout << " pt = " << ptB << " y = " << yP << endl;
	  if ( fabs(yB) < 2.2 && ptB > 18. ) { Bjet18 = Bjet18 + weight; } 
	  if ( fabs(yB) < 2.2 && ptB > 32. ) { Bjet32 = Bjet32 + weight; }       
	}
     
	if(j.momentum().pT()>30 && fabs(j.momentum().rapidity())<2.4){
	  
	  foreach (const Particle& muon, muons.particles()) {
	  
	    if ( fabs(muon.momentum().eta()) < 2.4 && muon.momentum().pT() > 9 ) muongoodevent=true;
	  }
	  
	  const double ptB= j.momentum().pT();
          const double yB= j.momentum().rapidity();
	
	  if(muongoodevent && onebtag){
	    Bjetmuon+=weight;
	    _h_dsigdpt->fill( ptB, weight );
	    _h_dsigdy->fill( fabs(yB), weight );	
	  }
	}
      }
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      /// @todo Normalise, scale and otherwise manipulate histograms here
      
      // scale(_h_YYYY, crossSection()/sumOfWeights()); # norm to cross section
      
      cout << " CMS_2012_S9383875 [pb] gen xsec = " << crossSection() << endl;
      cout << " CMS_2012_S9383875  Bjet > 18 GeV  xsec = " << Bjet18/sumOfWeights()*crossSection() << endl;
      cout << " CMS_2012_S9383875  Bjet > 32 GeV  xsec = " << Bjet32/sumOfWeights()*crossSection() << endl;

      double invlumi = crossSection()/picobarn/sumOfWeights();

      scale(_h_dsigdpty05, invlumi); 
      scale(_h_dsigdpty10, invlumi); 
      scale(_h_dsigdpty15, invlumi); 
      scale(_h_dsigdpty20, invlumi); 
      scale(_h_dsigdpty22, invlumi/0.4); 

      double invlumiNano = crossSection()/nanobarn/sumOfWeights();

      cout << " CMS_2012_S9383875  Bjet > 30 GeV with muon  xsec = " << Bjetmuon/sumOfWeights()*crossSection()/nanobarn << endl;

      scale(_h_dsigdpt, invlumiNano); 
      //scale(_h_dsigdy, invlumiNano/0.5); 
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
