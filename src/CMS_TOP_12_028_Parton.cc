#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "GeneratorInterface/RivetTop/interface/PartonTop.hh"

namespace Rivet {

class CMS_TOP_12_028_Parton : public Analysis {
public:
  CMS_TOP_12_028_Parton() : Analysis("CMS_TOP_12_028_Parton") {
  }

  void init() {
    // Parton level top quarks
    PartonTop ttbarState;
    addProjection(ttbarState, "ttbar");

    _h14_diffXSecTopSemiLepPartontopPt         = bookHisto1D("h14_diffXSecTopSemiLepPartontopPt"        );
    _h15_diffXSecTopSemiLepPartontopPtTtbarSys = bookHisto1D("h15_diffXSecTopSemiLepPartontopPtTtbarSys");
    _h16_diffXSecTopSemiLepPartontopY          = bookHisto1D("h16_diffXSecTopSemiLepPartontopY"         );
    _h17_diffXSecTopSemiLepPartonttbarDelPhi   = bookHisto1D("h17_diffXSecTopSemiLepPartonttbarDelPhi"  );
    _h18_diffXSecTopSemiLepPartontopPtLead     = bookHisto1D("h18_diffXSecTopSemiLepPartontopPtLead"    );
    _h19_diffXSecTopSemiLepPartontopPtSubLead  = bookHisto1D("h19_diffXSecTopSemiLepPartontopPtSubLead" );
    _h20_diffXSecTopSemiLepPartonttbarPt       = bookHisto1D("h20_diffXSecTopSemiLepPartonttbarPt"      );
    _h21_diffXSecTopSemiLepPartonttbarY        = bookHisto1D("h21_diffXSecTopSemiLepPartonttbarY"       );
    _h22_diffXSecTopSemiLepPartonttbarMass     = bookHisto1D("h22_diffXSecTopSemiLepPartonttbarMass"    );

    _h23_diffXSecTopDiLepPartontopPt           = bookHisto1D("h23_diffXSecTopDiLepPartontopPt"          );
    _h24_diffXSecTopDiLepPartontopPtTtbarSys   = bookHisto1D("h24_diffXSecTopDiLepPartontopPtTtbarSys"  );
    _h25_diffXSecTopDiLepPartontopY            = bookHisto1D("h25_diffXSecTopDiLepPartontopY"           );
    _h26_diffXSecTopDiLepPartonttbarDelPhi     = bookHisto1D("h26_diffXSecTopDiLepPartonttbarDelPhi"    );
    _h27_diffXSecTopDiLepPartontopPtLead       = bookHisto1D("h27_diffXSecTopDiLepPartontopPtLead"      );
    _h28_diffXSecTopDiLepPartontopPtSubLead    = bookHisto1D("h28_diffXSecTopDiLepPartontopPtSubLead"   );
    _h29_diffXSecTopDiLepPartonttbarPt         = bookHisto1D("h29_diffXSecTopDiLepPartonttbarPt"        );
    _h30_diffXSecTopDiLepPartonttbarY          = bookHisto1D("h30_diffXSecTopDiLepPartonttbarY"         );
    _h31_diffXSecTopDiLepPartonttbarMass       = bookHisto1D("h31_diffXSecTopDiLepPartonttbarMass"      );
  };

  void analyze(const Event& event) {
    const double weight = event.weight();

    // Get the parton level ttbar candidate
    const PartonTop& ttbarState = applyProjection<PartonTop>(event, "ttbar");
    const Particle& tCand1 = ttbarState.t1();
    const Particle& tCand2 = ttbarState.t2();

    // Do the anlaysis only for semileptonic and full leptonic channels.
    // Veto tau decays
    if ( ttbarState.mode() != PartonTop::CH_SEMILEPTON and
         ttbarState.mode() != PartonTop::CH_FULLLEPTON ) vetoEvent;
    if ( ttbarState.mode1() >= PartonTop::CH_TAU_HADRON ||
         ttbarState.mode2() >= PartonTop::CH_TAU_HADRON ) vetoEvent;

    // Fill top quarks they are defined in the parton level, full phase space
    const FourMomentum& t1P4 = tCand1.momentum();
    const FourMomentum& t2P4 = tCand2.momentum();
    const double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
    const FourMomentum ttbarP4 = t1P4+t2P4;
    const FourMomentum t1P4AtCM = LorentzTransform(-ttbarP4.boostVector()).transform(t1P4);
    const double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

    if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
      _h14_diffXSecTopSemiLepPartontopPt->fill(t1Pt, weight); 
      _h14_diffXSecTopSemiLepPartontopPt->fill(t2Pt, weight); 
      _h15_diffXSecTopSemiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight); 
      _h16_diffXSecTopSemiLepPartontopY->fill(t1P4.rapidity(), weight); 
      _h16_diffXSecTopSemiLepPartontopY->fill(t2P4.rapidity(), weight); 
      _h17_diffXSecTopSemiLepPartonttbarDelPhi->fill(dPhi, weight); 
      _h18_diffXSecTopSemiLepPartontopPtLead->fill(std::max(t1Pt, t2Pt), weight); 
      _h19_diffXSecTopSemiLepPartontopPtSubLead->fill(std::min(t1Pt, t2Pt), weight); 
      _h20_diffXSecTopSemiLepPartonttbarPt->fill(ttbarP4.pT(), weight); 
      _h21_diffXSecTopSemiLepPartonttbarY->fill(ttbarP4.rapidity(), weight); 
      _h22_diffXSecTopSemiLepPartonttbarMass->fill(ttbarP4.mass(), weight); 
    }
    else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
      _h23_diffXSecTopDiLepPartontopPt->fill(t1Pt, weight);
      _h23_diffXSecTopDiLepPartontopPt->fill(t2Pt, weight);
      _h24_diffXSecTopDiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t1P4.rapidity(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t2P4.rapidity(), weight);
      _h26_diffXSecTopDiLepPartonttbarDelPhi->fill(dPhi, weight);
      _h27_diffXSecTopDiLepPartontopPtLead->fill(std::max(t1Pt, t2Pt), weight);
      _h28_diffXSecTopDiLepPartontopPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
      _h29_diffXSecTopDiLepPartonttbarPt->fill(ttbarP4.pT(), weight);
      _h30_diffXSecTopDiLepPartonttbarY->fill(ttbarP4.rapidity(), weight);
      _h31_diffXSecTopDiLepPartonttbarMass->fill(ttbarP4.mass(), weight);
    }
  };

  void finalize() {
    normalize(_h14_diffXSecTopSemiLepPartontopPt        );
    normalize(_h15_diffXSecTopSemiLepPartontopPtTtbarSys);
    normalize(_h16_diffXSecTopSemiLepPartontopY         );
    normalize(_h17_diffXSecTopSemiLepPartonttbarDelPhi  );
    normalize(_h18_diffXSecTopSemiLepPartontopPtLead    );
    normalize(_h19_diffXSecTopSemiLepPartontopPtSubLead );
    normalize(_h20_diffXSecTopSemiLepPartonttbarPt      );
    normalize(_h21_diffXSecTopSemiLepPartonttbarY       );
    normalize(_h22_diffXSecTopSemiLepPartonttbarMass    );

    normalize(_h23_diffXSecTopDiLepPartontopPt        );
    normalize(_h24_diffXSecTopDiLepPartontopPtTtbarSys);
    normalize(_h25_diffXSecTopDiLepPartontopY         );
    normalize(_h26_diffXSecTopDiLepPartonttbarDelPhi  );
    normalize(_h27_diffXSecTopDiLepPartontopPtLead    );
    normalize(_h28_diffXSecTopDiLepPartontopPtSubLead );
    normalize(_h29_diffXSecTopDiLepPartonttbarPt      );
    normalize(_h30_diffXSecTopDiLepPartonttbarY       );
    normalize(_h31_diffXSecTopDiLepPartonttbarMass    );
  };

private:
  Histo1DPtr _h14_diffXSecTopSemiLepPartontopPt        ;
  Histo1DPtr _h15_diffXSecTopSemiLepPartontopPtTtbarSys;
  Histo1DPtr _h16_diffXSecTopSemiLepPartontopY         ;
  Histo1DPtr _h17_diffXSecTopSemiLepPartonttbarDelPhi  ;
  Histo1DPtr _h18_diffXSecTopSemiLepPartontopPtLead    ;
  Histo1DPtr _h19_diffXSecTopSemiLepPartontopPtSubLead ;
  Histo1DPtr _h20_diffXSecTopSemiLepPartonttbarPt      ;
  Histo1DPtr _h21_diffXSecTopSemiLepPartonttbarY       ;
  Histo1DPtr _h22_diffXSecTopSemiLepPartonttbarMass    ;

  Histo1DPtr _h23_diffXSecTopDiLepPartontopPt        ;
  Histo1DPtr _h24_diffXSecTopDiLepPartontopPtTtbarSys;
  Histo1DPtr _h25_diffXSecTopDiLepPartontopY         ;
  Histo1DPtr _h26_diffXSecTopDiLepPartonttbarDelPhi  ;
  Histo1DPtr _h27_diffXSecTopDiLepPartontopPtLead    ;
  Histo1DPtr _h28_diffXSecTopDiLepPartontopPtSubLead ;
  Histo1DPtr _h29_diffXSecTopDiLepPartonttbarPt      ;
  Histo1DPtr _h30_diffXSecTopDiLepPartonttbarY       ;
  Histo1DPtr _h31_diffXSecTopDiLepPartonttbarMass    ;

};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_TOP_12_028_Parton> plugin_CMS_TOP_12_028_Parton;

}
