// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
// @Author: Sercan Sen <ssen@cern.ch>

namespace Rivet {

  class CMS_FSQ_15_005 : public Analysis {
  public:

    CMS_FSQ_15_005()
      : Analysis("CMS_FSQ_15_005"), _xi_sd_cut(1E-6), _xi_hf_cut(1E-6), _xi_castor_cut(1E-7)
    {    }

  public:

    void init() {
        addProjection(FinalState(),"FS");
        addProjection(FinalState(-6.6, -5.2, 0.0*GeV), "CASTORM");
        addProjection(FinalState(-5.205, -3.152, 0.0*GeV), "HFM");
        addProjection(FinalState(3.152, 5.205, 0.0*GeV), "HFP");

        //_h_xsec = bookHisto1D("xsec",3,0.,3.);
        _h_xsec = bookHisto1D(1, 1, 1); // TODO: Use this auto-booking system when the data points in yoda 

    }

    void analyze(const Event& event) {

        const double weight = event.weight();

        const FinalState& fs = applyProjection<FinalState>(event, "FS");
        // Calculate gap sizes and midpoints
        const ParticleVector particlesByRapidity = fs.particlesByPt();
        const size_t num_particles = particlesByRapidity.size();

        vector<double> gaps;
        vector<double> midpoints;
        for (size_t ip = 1; ip < num_particles; ++ip) {
          const Particle& p1 = particlesByRapidity[ip-1];
          const Particle& p2 = particlesByRapidity[ip];
          const double gap = p2.momentum().rapidity() - p1.momentum().rapidity();
          const double mid = (p2.momentum().rapidity() + p1.momentum().rapidity()) / 2.;
          gaps.push_back(gap);
          midpoints.push_back(mid);
        }

        // Separate X and Y systems
        int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
        double gapcenter = midpoints[imid];

        FourMomentum MxFourVector(0.,0.,0.,0.);
        FourMomentum MyFourVector(0.,0.,0.,0.);

        foreach(const Particle& p, fs.particlesByPt()) {
            if (p.momentum().rapidity() < gapcenter) { // X system is at minus side.
                MxFourVector += p.momentum();
            } else {
                MyFourVector += p.momentum(); // Y system is at plus side.
            }
        }

        double Mx = MxFourVector.mass();
        double My = MyFourVector.mass();

        double xix = (Mx * Mx) / (sqrtS()/GeV * sqrtS()/GeV);
        double xiy = (My * My) / (sqrtS()/GeV * sqrtS()/GeV);
        double xi_sd = std::max(xix, xiy);

        //sigmainel // Total inelastic cross section. Eq.8 in the paper. This is 71.26 +-... mb.
        _h_xsec->fill(0.5, weight);
        //sigmaHF // Eg. 6 in the paper. 65.77 +- ... mb. 
        if (xi_sd > _xi_sd_cut) {  // This means either xix or xiy greater than 1e-6. "if (xix > _xi_hf_cut || xiy > _xi_hf_cut) {"
            _h_xsec->fill(1.5, weight);
        }
        //sigmaHFCASTOR // Eq. 7 in the paper. 66.85 +- ... mb.  p.s. X(Y) system is the minus(plus) side. CASTOR is only at minus side. Therefore, X is limited by the CASTOR acceptance (1e-7) and Y is limited by HF acceptance (1e-6). 
        if (xix > _xi_castor_cut || xiy > _xi_hf_cut) { 
            _h_xsec->fill(2.5, weight);
        }

    } // end of events loop
        
    void finalize() {

        scale(_h_xsec, crossSection()/millibarn/sumOfWeights());

    }

  private:

    Histo1DPtr _h_xsec;
    double _xi_sd_cut;
    double _xi_hf_cut;
    double _xi_castor_cut;

  };

  DECLARE_RIVET_PLUGIN(CMS_FSQ_15_005);

}
