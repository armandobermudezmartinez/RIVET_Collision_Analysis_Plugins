#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  class CMS_2017_I1511284 : public Analysis {
  public:

    /// Constructor
    CMS_2017_I1511284()
      : Analysis("CMS_2017_I1511284")
    {    }

    /// Book histograms and initialise projections before the run
    void init() {

    addProjection(FinalState(), "FS");
    _h_totEnergy = bookHisto1D(1, 1, 1);
    _h_emEnergy = bookHisto1D(2, 1, 1);
    _h_hadEnergy = bookHisto1D(3, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        const double weight = event.weight();

        const FinalState& fs = applyProjection<FinalState>(event, "FS");
        if (fs.size() < 2) vetoEvent; // need at least two particles to calculate gaps

        double gapCenter = 0.;
        double largestGap = 0.;
        double previousRapidity = 0.;
        bool first = true;

        foreach(const Particle& p, fs.particles(cmpMomByRap)) {
            if (first) { // First particle
                first = false;
                previousRapidity = p.rapidity();
            } else {
                double gap = fabs(p.rapidity()-previousRapidity);
                if (gap > largestGap) {
                   largestGap = gap; // largest gap
                   gapCenter = (p.rapidity()+previousRapidity)/2.; // find the center of the gap to separate the X and Y systems.
                }
                previousRapidity = p.rapidity();
            }
        }

        FourMomentum mxFourVector, myFourVector;
        foreach(const Particle& p, fs.particles(cmpMomByRap)) {
            ((p.rapidity() > gapCenter) ? mxFourVector : myFourVector) += p.momentum();
        }
        const double xiX = mxFourVector.mass2()/sqr(sqrtS());
        const double xiY = myFourVector.mass2()/sqr(sqrtS());
        const double xi = max(xiX,xiY);
        if (xi < 1e-6) vetoEvent;

        double totEnergy = 0.;
        double emEnergy = 0.;
        double hadEnergy = 0.;
        foreach(const Particle& p, fs.particles(cmpMomByRap)) {
            if (p.eta()>-6.6 && p.eta()<-5.2){
                if (p.isVisible and p.absPid() != 13){
                    totEnergy += p.energy();
                    if ( p.abspid() == 11 || p.abspid() == 22 || p.abspid() == 111){
                        emEnergy += p.energy();
                    }
                    if ( p.abspid() != 11 && p.abspid() != 22 && p.abspid() != 111){
                        hadEnergy += p.energy();
                    }
                }
            }
        }

        _h_totEnergy->fill(totEnergy, weight);
        _h_emEnergy->fill(emEnergy, weight);
        _h_hadEnergy->fill(hadEnergy, weight);
    }

    void finalize() {

        scale(_h_totEnergy, crossSection()/microbarn/sumOfWeights());
        scale(_h_emEnergy, crossSection()/microbarn/sumOfWeights());
        scale(_h_hadEnergy, crossSection()/microbarn/sumOfWeights());
    }

  private:

    Histo1DPtr _h_totEnergy;
    Histo1DPtr _h_emEnergy;
    Histo1DPtr _h_hadEnergy;

  };

  DECLARE_RIVET_PLUGIN(CMS_2017_I1511284);
}
