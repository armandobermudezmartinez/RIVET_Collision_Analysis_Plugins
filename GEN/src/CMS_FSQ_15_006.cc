// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class CMS_FSQ_15_006 : public Analysis {
  public:

    CMS_FSQ_15_006()
      : Analysis("CMS_FSQ_15_006"), _noe_inel(0.), _noe_nsd(0.), _noe_bsc(0.)
    {    }

  public:

    void init() {

      addProjection(FinalState(),"FS");

      const ChargedFinalState cfsBSCplus(3.9, 4.4, 0*MeV);
      addProjection(cfsBSCplus, "cfsBSCplus");

      const ChargedFinalState cfsBSCminus(-4.4, -3.9, 0*MeV);
      addProjection(cfsBSCminus, "cfsBSCminus");

      _h_inel = bookHisto1D(1, 1, 1);
      _h_nsd  = bookHisto1D(2, 1, 1);
      _h_et   = bookHisto1D(3, 1, 1);

    }

    void analyze(const Event& event) {

        const double weight = event.weight();

        double ybeam = 9.54;

        bool bscplus  = true;
        bool bscminus = true;

        const ChargedFinalState& cfsBSCplus = applyProjection<ChargedFinalState>(event, "cfsBSCplus");
        if (cfsBSCplus.empty()) bscplus = false;
        const ChargedFinalState& cfsBSCminus = applyProjection<ChargedFinalState>(event, "cfsBSCminus");
        if (cfsBSCminus.empty()) bscminus = false;

        const FinalState& fs = applyProjection<FinalState>(event, "FS");

        const ParticleVector particlesByPt = fs.particlesByPt();
        const size_t num_particles = particlesByPt.size();

        vector<double> gaps;
        vector<double> midpoints;
        for (size_t ip = 1; ip < num_particles; ++ip) {
          const Particle& p1 = particlesByPt[ip-1];
          const Particle& p2 = particlesByPt[ip];
          const double gap = p2.momentum().rapidity() - p1.momentum().rapidity();
          const double mid = (p2.momentum().rapidity() + p1.momentum().rapidity()) / 2.;
          gaps.push_back(gap);
          midpoints.push_back(mid);
        }

        int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
        double gapcenter = midpoints[imid];

        FourMomentum MxFourVector(0.,0.,0.,0.);
        FourMomentum MyFourVector(0.,0.,0.,0.);

        foreach(const Particle& p, fs.particlesByPt()) {
            if (p.momentum().rapidity() > gapcenter) {
                MxFourVector += p.momentum();
            } else {
                MyFourVector += p.momentum();
            }
        }

        double Mx = MxFourVector.mass();
        double My = MyFourVector.mass();

        double xix = (Mx * Mx) / (sqrtS()/GeV * sqrtS()/GeV);
        double xiy = (My * My) / (sqrtS()/GeV * sqrtS()/GeV);
        double xi_sd = std::max(xix, xiy);

        bool inel = false;

        if (xi_sd > 1E-6) { 
            inel = true;
            ++_noe_inel;    
        }

        bool nsd = false;
        bool bsc = false;
        int nplus  = 0;    
        int nminus = 0;    

        if (bscplus && bscminus) { 
            ++_noe_bsc;
            bsc = true;
        }

        foreach(const Particle& p, fs.particlesByPt()) {
            double eta = p.momentum().eta();
            double tenergy = p.momentum().Et();
            if (abs(p.pid()) >= 12 && abs(p.pid()) <= 16) continue;
            if (eta > 2.866 && eta < 5.205) ++nplus;
            if (eta < -2.866 && eta > -5.205) ++nminus;

            if (bsc) { 
                _h_et->fill(eta-ybeam, tenergy*weight);
            }
        }

        if (nminus > 0 && nplus > 0) { 
            nsd = true;
            ++_noe_nsd;
        }

        foreach(const Particle& p, fs.particlesByPt()) {
            double eta = p.momentum().eta();
            double energy = p.momentum().E();
            if (inel) { 
                _h_inel->fill(eta, energy*weight);
            }
            if (nsd) { 
                _h_nsd->fill(eta, energy*weight);
            }
        }

    }

    void finalize() {

       scale(_h_inel, 1./_noe_inel);
       scale(_h_nsd, 1./_noe_nsd);
       scale(_h_et, 1./_noe_bsc);

    }

  private:

    Histo1DPtr _h_inel;
    Histo1DPtr _h_nsd;
    Histo1DPtr _h_et;
    double _noe_inel;
    double _noe_nsd;
    double _noe_bsc;

  };


  DECLARE_RIVET_PLUGIN(CMS_FSQ_15_006);

}
