#include "TopMonteCarlo/RivetTop/interface/CMSGenParticle.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

void CMSGenParticle::project(const Event& e) {
  _theParticles.clear();

  // Collect particles from resonance
  const std::vector<GenParticle*> allParticles = Rivet::particles(e.genEvent());
  std::vector<GenParticle*> pFromReso;
  pFromReso.reserve(allParticles.size());
  foreach (GenParticle* p, allParticles) {
    if ( p->status() == 1 or p->status() == 4 ) continue;
    if ( !isResonance(std::abs(p->pdg_id())) ) continue;
    if ( !p->production_vertex() or !p->end_vertex() ) continue; // For sately

    // Skip if this resonance is decay product of other resonance
    bool isDecayFromResonance = false;
    foreach (GenParticle* gp, Rivet::particles(p->production_vertex(), HepMC::parents) ) {
      if ( isResonance(std::abs(gp->pdg_id())) ) {
        isDecayFromResonance = true;
        break;
      }
    }
    if ( isDecayFromResonance ) continue;

    // Collect all stable particles from this resonance particle
    foreach (GenParticle* sp, Rivet::particles(p->end_vertex(), HepMC::descendants)) {
      if ( sp->status() != 1 ) continue;
      pFromReso.push_back(sp);
    }
  }
  std::sort(pFromReso.begin(), pFromReso.end());

  // Collect stable particles vetoing uninteresting ones
  foreach (GenParticle* p, allParticles) {
    if ( p->status() != 1 ) continue;

    const int aid = std::abs(p->pdg_id());
    // Skip exotic particles
    if ( std::binary_search(_vetoIds.begin(), _vetoIds.end(), aid) ) continue;
    // Skip muons and neutrinos from resonance
    if ( std::binary_search(_vetoIdsFromResonances.begin(), _vetoIdsFromResonances.end(), aid) &&
         std::binary_search(pFromReso.begin(), pFromReso.end(), p) ) continue;

    _theParticles.push_back(Particle(*p));
  }
}
