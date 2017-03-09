// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {

  using namespace Cuts;

  class CMS_2015_I1380605 : public Analysis {
  public:
  
    /// Constructor
    CMS_2015_I1380605()
      : Analysis("CMS_2015_I1380605")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      ntracks  = 0;

      /// @todo Initialise and register projections here
      const ChargedFinalState cfs(-7., 7., 0.0*GeV);
      addProjection(cfs, "CFS");
      addProjection(FastJets(cfs, FastJets::ANTIKT, 0.5),"Jets");

      /// @todo Book histograms here, e.g.:
      _h_tracks = bookHisto1D(1, 1, 1);
      _h_jets = bookHisto1D(2, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      // do track analysis here
      const ChargedFinalState & cfs = applyProjection<ChargedFinalState>(event, "CFS");
      int count_plus = 0, count_minus = 0;

      foreach (const Particle& p, cfs.particles()) {
        if (inRange(p.eta(),  5.3,  6.5)) count_plus++;
        if (inRange(p.eta(), -6.5, -5.3)) count_minus++;
      }
      const bool cutsor  = (count_plus > 0 || count_minus > 0);
      
      if (cutsor)  { 
        ntracks = ntracks + weight; 
        //find pttrackmax
        double track_ptmax = -99.9;

        if (cfs.size() > 0) { 
          foreach (const Particle& j, cfs.particles()) {
            if (inRange(j.eta(),  -2.4,  2.4)) { 
              if (j.momentum().pT() > track_ptmax) track_ptmax = j.momentum().pT();
            } 
          }
        }

        for (size_t i = 0; i < _h_tracks->numBins(); ++i) {
          double binlimitlow_t = _h_tracks->bin(i).xMin();
          double weightbw_t = weight * _h_tracks->bin(i).xWidth();
          double xbin_t = _h_tracks->bin(i).xMid();
        
          if (track_ptmax > binlimitlow_t) { _h_tracks -> fill(xbin_t, weightbw_t);}   
        }

        // do jet analysis here
        const FastJets &fj = applyProjection<FastJets>(event,"Jets");
        const Jets jetsdeta = fj.jetsByPt(1.00*GeV);

        //find ptjetmax
        double jet_ptmax = -99.9;
        if (jetsdeta.size() > 0)   { 
          foreach(const Jet &j, jetsdeta) {
            if (inRange(j.eta(),  -1.9,  1.9) and inRange(j.pT(),1.0,60.0)) { 
              if (j.momentum().pT() > jet_ptmax) jet_ptmax = j.momentum().pT();
            }
          }
        }

        for (size_t i = 0; i < _h_jets->numBins(); ++i) {
          double binlimitlow_j = _h_jets->bin(i).xMin();
          double weightbw_j = weight * _h_jets->bin(i).xWidth();
          double xbin_j = _h_jets->bin(i).xMid();

          if (jet_ptmax > binlimitlow_j) { _h_jets -> fill(xbin_j, weightbw_j);}   
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      double norm_t0 =  _h_tracks->bin(7).height()/2.056170e-03;
      double norm_t1 =  _h_tracks->bin(7).sumW()/2.056170e-03;
      double norm_j0 =  _h_jets->bin(13).height()/3.575290e-03;
      double norm_j1 =  _h_jets->bin(13).sumW()/3.575290e-03;
      cout << " norm track " << norm_t0 << " " << norm_t1 << endl;
      cout << " norm  jets " << norm_j0 << " " << norm_j1 << endl;
      if (norm_t0 > 0 ) scale(_h_tracks, 1./ norm_t0);
      if (norm_j0 > 0 ) scale(_h_jets, 1./ norm_j0);
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here

    /// @name Histograms
    //@{
    Histo1DPtr _h_tracks, _h_jets;
    double ntracks;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1380605);

}
