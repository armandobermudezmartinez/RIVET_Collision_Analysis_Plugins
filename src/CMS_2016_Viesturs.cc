#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "TopMonteCarlo/RivetTop/interface/CMSGenParticle.hh"
#include "TopMonteCarlo/RivetTop/interface/PartonTop.hh"
#include "TopMonteCarlo/RivetTop/interface/CMS_2016_Viesturs.hh"


#include "TMath.h"
#include "TString.h"
namespace Rivet 
{
  DECLARE_RIVET_PLUGIN(CMS_2016_Viesturs);

}
void CMS_2016_Viesturs::init()
{
  printf("*** running analysis CMS_2016_Viesturs ***\n");
  //      getchar();
  /*      FourMomentum mom(10, 12, 14, 100);
          Jet j(mom);
          Jet *j_ptr = &j;
          FourMomentum *mom_ptr = (FourMomentum *)j_ptr;
          printf("jet px %f py %f pz %f E %f\n", j_ptr -> px(), j_ptr -> py(), j_ptr -> pz(), j_ptr -> E());
          printf("mom px %f py %f pz %f E %f\n", mom_ptr -> px(), mom_ptr -> py(), mom_ptr -> pz(), mom_ptr -> E());
          getchar();*/
  _vis_unit_weights = 0;
  _full_unit_weights = 0;
  _plots_1D_ptr["btag_jets_pt"]        = bookHisto1D("btag_jets_pt", 50, 20, 300);
  _plots_1D_ptr["light_jets_pt"]       = bookHisto1D("light_jets_pt", 50, 20, 300);

  _plots_1D_ptr["electrons_pt"]        = bookHisto1D("electrons_pt", 50, 20, 300);
  addProjection(FinalState(), "FS");
  addProjection(FastJets(FinalState(-2.5, 2.5, 0*GeV), FastJets::ANTIKT, 0.5), "Jets");
  AssignHistograms();

  _cfat.plots_ptr_ = &_plots_1D_ptr;
}

void CMS_2016_Viesturs::analyze(const Event& event) 
{
  vector<const Particle *> electrons;
  vector<const Particle *> muons;
  const FinalState& fs = applyProjection<FinalState>(event, "FS");
  const Particles & particles = fs.particles();
  const double weight = event.weight();
  if (weight != 0.) _vis_unit_weights += weight/std::abs(weight);
  if (weight != 0.) _full_unit_weights += weight/std::abs(weight);
  _cfat.weight_ = weight;
  for (Particles::const_iterator it = particles.cbegin(); it != particles.cend(); it++)
  {
    const PdgId absid = it -> abspid();
    if (absid == 11)
      if (it -> pt() > 15*GeV and it -> abseta() < 2.5)
        electrons.push_back(&*it);
    if (absid == 13)
      if (it -> pt() > 10*GeV and it -> abseta() < 2.4)
        muons.push_back(&*it);

  }
  if (electrons.size() != 1 or muons.size() != 0)
    return;

  _plots_1D_ptr["electrons_pt"] -> fill(electrons[0] ->  pt(), weight);
  const Jets& jets = applyProjection<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);

  vector<const Jet *> btag_jets;
  vector<const Jet *> light_jets;

  foreach( const Jet& jet, jets ) 
  {
    if (jet.bTagged())
      btag_jets.push_back(&jet);
    else
      light_jets.push_back(&jet);

  }
  if (btag_jets.size() != 2 or light_jets.size() !=2)
    return;
  _cfat._b_jets_ptr = &btag_jets;
  _cfat._light_jets_ptr = &light_jets;
  _cfat._jets = &jets;
  _cfat.Work();
  for (vector<const Jet *>::const_iterator cit = btag_jets.cbegin(); cit != btag_jets.cend(); cit ++)
  {
    _plots_1D_ptr["btag_jets_pt"] -> fill((*cit) -> pt(), weight);
  }
  for (vector<const Jet *>::const_iterator cit = light_jets.cbegin(); cit != light_jets.cend(); cit ++)
  {
    _plots_1D_ptr["light_jets_pt"] -> fill((*cit) -> pt(), weight);
  }
}

void CMS_2016_Viesturs::finalize() 
{
  double ttbarXS = 0.;
  if (!std::isnan(crossSectionPerEvent())) 
  {
    std::cout << "Using generator cross section: " << crossSection() << " pb" << std::endl;
    ttbarXS = crossSection();
  }
  else 
  {
    std::cout << "No valid cross section given, using NNLO value: 252.89 pb" << std::endl;
    ttbarXS = 252.89; // NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV)
    // see also https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
  }

  std::cout << "Sum vis unit weights: " << _vis_unit_weights
    << ", sum full unit weights: " << _full_unit_weights << std::endl;
  for (map<TString, Histo1DPtr>::const_iterator cit = _plots_1D_ptr.cbegin(); cit != _plots_1D_ptr.cend(); cit ++)
  {
    printf("number of entries %lu\n", cit -> second -> numEntries());
    scale(cit -> second,       ttbarXS / _vis_unit_weights);
    normalize(cit -> second);
  }
}

//    void AssignHistograms();
/*  };
    DECLARE_RIVET_PLUGIN(CMS_2016_Viesturs);



    }*/

