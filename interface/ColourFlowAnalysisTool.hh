#ifndef _ColourFlowAnalysisTool_hh_
#define _ColourFlowAnalysisTool_hh_
#ifndef _CMS_2016_Viesturs_hh_
#define _CMS_2016_Viesturs_hh_
#endif
#include "TLorentzVector.h"
#include "TH1.h"
#include "TString.h"
#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Analysis.hh"
 //#include "TopMonteCarlo/RivetTop/interface/CMS_2016_Viesturs.hh"

class TCanvas;
class TH2F;
#include <map>
namespace Rivet
{
  class CMS_2016_Viesturs;
};
class Event;

using namespace std;
using namespace Rivet;
class ColourFlowAnalysisTool
{

  friend class Rivet::CMS_2016_Viesturs;//::AssignHistograms();
  struct PullVector: public TVector2
  {
    PullVector(Double_t, Double_t);
    Double_t & phi_component = fX;
    Double_t & eta_component = fY;
  };
  static const unsigned char N_jet_types_;
  static const char* tag_channel_;
  static const char * tag_charge_types_[2];
  static const char* tag_jet_types_[];
  static const unsigned char N_DeltaR_types_;
  static const char * tag_DeltaR_types_[];
  static const unsigned char _N_levels;
  static const char * tag_levels_[];
  static const TLorentzVector beam_;
  TCanvas * canvas_;
  TH2F    * PtRadProf_;
public:
  unsigned char                work_mode_;
  unsigned char                event_display_mode_;
  unsigned char                PtRadiation_mode_;
  unsigned long                event_number_;
  ColourFlowAnalysisTool();
  const Rivet::Event                  * event_ptr_;
  const vector<const Jet *> * _light_jets_ptr;
  const vector<unsigned char>  * light_jets_indices_ptr_;
  unsigned char jet_indices_[2] = {255, 255}; 
  const vector<const Jet *> * _b_jets_ptr;
  const vector<unsigned char>  * b_jets_indices_ptr_;
  const TLorentzVector         * neutrino_ptr_;
  const Particle         * lepton_ptr_;
  const Jet         * leading_light_jet_ptr_;
  vector<const Jet *> vect_jets_;
  const Jets * _jets;
  float leading_light_jet_index_;
  float second_leading_light_jet_index_;
  const Jet         * second_leading_light_jet_ptr_;
  FourMomentum had_t_;
  FourMomentum lept_t_;
  
  map<TString, Histo1DPtr>     * plots_ptr_;
  //map<TString, TH1*>     & plots_ = *plots_ptr_;
  float weight_;
  FourMomentum GetChargedJet(const Jet &) const;
  vector<const Jet *> IdentifyJets() ;
  
  //  void AssignHistograms();
  void Work();
  float PullAngle(const PullVector & pull_vector, const TVector2 & jet_difference) const;
  PullVector CalculatePullVector(const Jet & jet, /*unsigned char index,*/ bool OnlyChargedConstituents) const; 
  void PlotAngleBetweenJets() const;
  void Do();
  void AnalyseParticleFlow() const;
  void EventDisplay(const PullVector &, float, const TVector2 &, bool) const;
  void PtRadiationProfile() const;
  void ConfigurePtRadiationProfileTool();
  void EndPtRadiationProfile();
};

#endif
