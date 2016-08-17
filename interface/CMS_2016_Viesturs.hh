#ifndef _CMS_2016_Viesturs_hh_
#define _CMS_2016_Viesturs_hh_
#include "TopMonteCarlo/RivetTop/interface/ColourFlowAnalysisTool.hh"
#include <map>
using namespace std;
class Event;
namespace Rivet
{
  class CMS_2016_Viesturs:
    public Analysis 
  {
    double _vis_unit_weights;
    double _full_unit_weights;
    map<TString, Histo1DPtr> _plots_1D_ptr;
    ColourFlowAnalysisTool _cfat;
    public:
    CMS_2016_Viesturs() : Analysis("CMS_2016_Viesturs") {}
    void init() ;
    void analyze(const Event& event) ;
    void finalize();
    void AssignHistograms();
  };
}
#endif
