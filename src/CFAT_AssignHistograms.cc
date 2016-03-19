#include "TMath.h"
#include "TopMonteCarlo/RivetTop/interface/CMS_2016_Viesturs.hh"

/*namespace Rivet
{
*/
void  CMS_2016_Viesturs::AssignHistograms() 
  {
    static const unsigned char ncategories_PA         = 2;
    static const char * cat_title_PA[ncategories_PA]     = {"pull_angle",     "cos_pull_angle" };
    static const unsigned short nbins_PA[ncategories_PA] = {25,               100              };
    static const double min_PA[ncategories_PA]           = {-1.2*TMath::Pi(), -1.2             };
    static const double max_PA[ncategories_PA]           = {- min_PA[0],         -min_PA[1]          };
    //static const char* axes_PA[ncategories_PA]           = {"; pull angle [rad]; Events",  "; cos; Events"  };
    for (unsigned char cat_index = 0; cat_index < ncategories_PA; cat_index++)
      {
	for (unsigned char level_index = 0; level_index < ColourFlowAnalysisTool::_N_levels; level_index ++)
	  {
	    for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	      {
		for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
		  {
		    for (unsigned char jet2_index = 0; jet2_index < ColourFlowAnalysisTool::N_jet_types_; jet2_index ++)
		      {
			for (unsigned char DeltaR_index = 0; DeltaR_index < ColourFlowAnalysisTool::N_DeltaR_types_; DeltaR_index ++)
			  {
	
			    if (jet1_index == jet2_index)
			      continue;
			    const TString hash_key = TString(cat_title_PA[cat_index]) + "_" + 
			      ColourFlowAnalysisTool::tag_charge_types_[charge_index] + "_" + 
			      ColourFlowAnalysisTool::tag_levels_[level_index] + "_" +
			      ColourFlowAnalysisTool::tag_jet_types_[jet1_index] + "_:_" + 
			      ColourFlowAnalysisTool::tag_jet_types_[jet2_index] + "_" + 
			      ColourFlowAnalysisTool::tag_DeltaR_types_[DeltaR_index] + "_" + 
			      ColourFlowAnalysisTool::tag_channel_;
			    _plots_1D_ptr . operator[](hash_key) = bookHisto1D(hash_key.Data(), 
									     /*hash_key + axes_PA[cat_index], */
									     nbins_PA[cat_index], 
									     min_PA[cat_index], 
									       max_PA[cat_index]);
			    printf("opening %s\n", hash_key.Data());
			  }
		      }
		  }
    
	      }
	  }
      }
    
    static const unsigned char ncategories_PV = 6;
    static const char * cat_title_PV[ncategories_PV]     = {"phi_PV",         "eta_PV",         "mag_PV",   "phi_jet_dif",  "eta_jet_dif", "mag_jet_dif" };
    static const unsigned short nbins_PV[ncategories_PV] = {100,              100,              100,        100,    100, 100    };
    static const double min_PV[ncategories_PV]           = {-0.01*TMath::Pi(), -0.05,               0,      -1.1*TMath::Pi(), -5.5, 0       };
    static const double max_PV[ncategories_PV]           = {- min_PV[0],         -min_PV[1],    0.05,       1.1*TMath::Pi(), 5.5, 10      };
    //static const char* axes_PV[ncategories_PV]           = {"; #phi [rad]; Events",  "; #eta [a.u.]; Events", "; magn. [a.u.]; Events", "; #phi [rad]; Events",  "; #eta [a.u.]; Events", "; magn. [a.u.]; Events"};
    for (unsigned char cat_index = 0; cat_index < ncategories_PV; cat_index++)
      {
	for (unsigned char level_index = 0; level_index < ColourFlowAnalysisTool::_N_levels; level_index ++)
	  { 
	    for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	      {
		for (unsigned char jet1_index = 0; jet1_index < 2; jet1_index ++)
		  {
		    const TString hash_key = TString(cat_title_PV[cat_index]) + "_" + 
		      ColourFlowAnalysisTool::tag_charge_types_[charge_index] + "_" + 
		      ColourFlowAnalysisTool::tag_levels_[level_index] + "_" +
		      ColourFlowAnalysisTool::tag_jet_types_[jet1_index] + "_" + 
		      ColourFlowAnalysisTool::tag_channel_;
		    _plots_1D_ptr . operator[](hash_key) = bookHisto1D(hash_key.Data(), 
								     /*hash_key + axes_PV[cat_index], */
								     nbins_PV[cat_index], 
								     min_PV[cat_index], 
								     max_PV[cat_index]);
    
		    printf("opening %s\n", hash_key.Data());
		  }
	      }
	  }
      }
  
    /*    static const unsigned char ncategories_CTRL              = 2;
    static const char * branch[2]                            = {"had", "lept"};
    static const char * cat_title_CTRL[ncategories_CTRL]     = {"W_mass",         "t_mass"    };
    static const unsigned short nbins_CTRL[ncategories_CTRL] = {100,              100         };
    static const double min_CTRL[ncategories_CTRL]           = {0,                50          };
    static const double max_CTRL[ncategories_CTRL]           = {150,              400         };
    static const char* axes_CTRL[ncategories_CTRL]           = {"; GeV; Events",  "; GeV; Events"};
    for (unsigned char branch_index = 0; branch_index < 2; branch_index ++)
      {
	for (unsigned char level_index = 0; level_index < 2; level_index ++)
	  {
	    for (unsigned char cat_index = 0; cat_index < ncategories_CTRL; cat_index ++)
	      { 
		const TString hash_key = TString(branch[branch_index]) + "_" + tag_levels_[level_index] + "_" + cat_title_CTRL[cat_index]; 
		plots_ptr_ -> operator[](hash_key) = bookHisto1D(hash_key, /*hash_key + axes_CTRL[cat_index],*//* nbins_CTRL[cat_index], min_CTRL[cat_index], max_CTRL[cat_index]);
	      }
	  }
      }
      													       */
    /*const char * particles[2] = {"all", "jet"};
  
    for (unsigned char DeltaR_index = 0; DeltaR_index < N_DeltaR_types_; DeltaR_index ++)
      {
	for (unsigned char level_index = 0; level_index < 2; level_index ++)
	  { 
	    for (unsigned char jet1_index = 0; jet1_index < N_jet_types_; jet1_index ++)
	      {
		for (unsigned char jet2_index = jet1_index + 1; jet2_index < N_jet_types_; jet2_index ++)
		  {
		    const TString hash_key = TString("angle_") + tag_jet_types_[jet1_index] + "_:_" + tag_jet_types_[jet2_index] + "_" + tag_levels_[level_index] + "_" + tag_DeltaR_types_[DeltaR_index];
		    plots_ptr_ -> operator[](hash_key) = bookHisto1D(hash_key, hash_key + "; angle [rad]; Events", 50, -0.1, TMath::Pi() + 0.1);
		  }
	      }
	    for (unsigned char charge_index = 0; charge_index < 2; charge_index ++)
	      {
		for (unsigned char particle_index = 0; particle_index < 2; particle_index ++)
		  {
		    const TString hash_key = TString("chi_") + tag_charge_types_[charge_index] + "_" + tag_levels_[level_index] + "_" + tag_DeltaR_types_[DeltaR_index] + "_" + particles[particle_index];
		    plots_ptr_ -> operator[](hash_key) = bookHisto1D(hash_key, hash_key + "; #chi; Events", 50, -0.1, 2.1);
		  }
	      }
	  }
      }
    plots_ptr_ -> operator[]("had_b_flavour")              = bookHisto1D("had_b_flavour", "had_b_flavour; flavour; Events", 21, -10.5, 10.5);
    plots_ptr_ -> operator[]("lept_b_flavour")             = bookHisto1D("lept_b_flavour", "lept_b_flavour; flavour; Events", 21, -10.5, 10.5);

    plots_ptr_ -> operator[]("leading_jet_flavour")        = bookHisto1D("leading_jet_flavour", "leading_jet_flavour; flavour; Events", 21, -10.5, 10.5);
    plots_ptr_ -> operator[]("second_leading_jet_flavour") = bookHisto1D("second_leading_jet_flavour", "second_leadiing_jet_flavour; flavour; Events", 21, -10.5, 10.5);
    //  ConfigurePtRadiationProfileTool();

    }*/
}
