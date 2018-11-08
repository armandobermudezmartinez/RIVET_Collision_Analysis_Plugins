// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "Rivet/Math/MathUtils.hh"

namespace Rivet {


  class CMS_2018_I1682495 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2018_I1682495()
      : Analysis("CMS_2018_I1682495"),
        _softdrop(fastjet::contrib::SoftDrop(0, 0.1, 0.8) ) // parameters are beta, zcut, R0
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // define a projection that keeps all the particles up to |eta|=5
      const FinalState fs(Cuts::abseta < 5.);

      // use FastJet, anti-kt(R=0.8) to do the clustering
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.8), "JetsAK8");

      // Histograms
      for (size_t i = 0; i < N_PT_BINS_dj; ++i ) {
        _h_ungroomedJetMass_dj[i][0] = bookHisto1D(i+1+0*N_PT_BINS_dj, 1, 1); // Ungroomed mass, absolute
        _h_sdJetMass_dj[i][0]        = bookHisto1D(i+1+1*N_PT_BINS_dj, 1, 1); // Groomed mass, absolute
        _h_ungroomedJetMass_dj[i][1] = bookHisto1D(i+1+2*N_PT_BINS_dj, 1, 1); // Ungroomed mass, normalized
        _h_sdJetMass_dj[i][1]        = bookHisto1D(i+1+3*N_PT_BINS_dj, 1, 1); // Groomed mass, normalized
      }

    }


    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin(double ptJ) {
      const double ptBins_dj[N_PT_BINS_dj+1] = { 200., 260., 350., 460., 550., 650., 760., 900, 1000, 1100, 1200, 1300, 13000};
      for (size_t ibin = 0; ibin < N_PT_BINS_dj; ++ibin) {
        if (inRange(ptJ, ptBins_dj[ibin], ptBins_dj[ibin+1])) return ibin;
      }
      return N_PT_BINS_dj;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Look at events with >= 2 jets
      auto jetsAK8 = applyProjection<FastJets>(event, "JetsAK8").jetsByPt(Cuts::pT > 200*GeV and Cuts::abseta < 2.4);
      if (jetsAK8.size() < 2) vetoEvent;

      // Get the leading two jets
      const fastjet::PseudoJet& j0 = jetsAK8[0].pseudojet();
      const fastjet::PseudoJet& j1 = jetsAK8[1].pseudojet();

      // Calculate delta phi and the pt asymmetry
      double deltaPhi = Rivet::deltaPhi( j0.phi(), j1.phi() );
      double ptasym = (j0.pt() - j1.pt()) / (j0.pt() + j1.pt());
      if (deltaPhi < 2.0 ) vetoEvent;
      if (ptasym > 0.3) vetoEvent;

      // Find the appropriate pT bins and fill the histogram
      const size_t njetBin0 = findPtBin(j0.pt()/GeV);
      const size_t njetBin1 = findPtBin(j1.pt()/GeV);
      if (njetBin0 < N_PT_BINS_dj && njetBin1 < N_PT_BINS_dj) {
        for ( size_t jbin = 0; jbin < N_CATEGORIES; jbin++ ){
          _h_ungroomedJetMass_dj[njetBin0][jbin]->fill(j0.m()/GeV, weight);
          _h_ungroomedJetMass_dj[njetBin1][jbin]->fill(j1.m()/GeV, weight);
        }
      }

      // Now run the substructure algs...
      fastjet::PseudoJet sd0 = _softdrop(j0);
      fastjet::PseudoJet sd1 = _softdrop(j1);
      // ... and repeat
      if (njetBin0 < N_PT_BINS_dj && njetBin1 < N_PT_BINS_dj) {
        for ( size_t jbin = 0; jbin < N_CATEGORIES; jbin++ ){
          _h_sdJetMass_dj[njetBin0][jbin]->fill(sd0.m()/GeV, weight);
          _h_sdJetMass_dj[njetBin1][jbin]->fill(sd1.m()/GeV, weight);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalize the normalized cross section histograms to unity,
      for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
        normalize(_h_ungroomedJetMass_dj[i][1]);
        normalize(_h_sdJetMass_dj[i][1]);
      }
      // Normalize the absolute cross section histograms to xs * lumi.
      for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
        scale(_h_ungroomedJetMass_dj[i][0],   crossSection()/picobarn / sumOfWeights());
        scale(_h_sdJetMass_dj[i][0],          crossSection()/picobarn / sumOfWeights());
      }
    }
    //@}


  private:

    /// @name FastJet grooming tools (configured in constructor init list)
    //@{
    const fastjet::contrib::SoftDrop _softdrop;
    //@}


    /// @name Histograms
    //@{
    enum { PT_200_260_dj=0, 
           PT_260_350_dj,
           PT_350_460_dj,
           PT_460_550_dj,
           PT_550_650_dj,
           PT_650_760_dj,
           PT_760_900_dj,
           PT_900_1000_dj,
           PT_1000_1100_dj,
           PT_1100_1200_dj,
           PT_1200_1300_dj,
           PT_1300_Inf_dj,
           N_PT_BINS_dj } BINS_dj;
    static const int N_CATEGORIES=2;
    Histo1DPtr _h_ungroomedJet0pt, _h_ungroomedJet1pt;
    Histo1DPtr _h_sdJet0pt, _h_sdJet1pt;
    // Here, store both the absolute (index 0) and normalized (index 1) cross sections. 
    Histo1DPtr _h_ungroomedJetMass_dj[N_PT_BINS_dj][2];
    Histo1DPtr _h_sdJetMass_dj[N_PT_BINS_dj][2];
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_I1682495);


}