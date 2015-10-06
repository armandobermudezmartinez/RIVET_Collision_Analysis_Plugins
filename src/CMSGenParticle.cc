#include "TopMonteCarlo/RivetTop/interface/CMSGenParticle.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

void CMSGenParticle::project(const Event& e) {
  _theParticles.clear();

  const std::vector<GenParticle*> particles = Rivet::particles(e.genEvent());
  std::set<GenParticle*> selected, invalid;

  foreach (GenParticle* p, particles ) {
    if ( invalid.find(p) != invalid.end() ) continue;
    if ( p->status() == 1 ) selected.insert(p);
  }

  foreach (GenParticle* p, selected ) {
    if ( invalid.find(p) != invalid.end() ) continue;

    if ( fromResonance(invalid, particles, p) ) {
      invalid.insert(p);
      continue;
    }

    if ( isIgnored(p->pdg_id()) ) continue;

    _theParticles.push_back(p);
  }
}

int CMSGenParticle::fromResonance(std::set<GenParticle*>& invalid, const std::vector<GenParticle*>& pv, GenParticle* p) const {
  const int id = p->pdg_id();
  const unsigned int aid = std::abs(id);

  if ( invalid.find(p) != invalid.end() ) return 2;
  if ( isResonance(aid) && p->status() == 3 ) return 1;
  if ( !isIgnored(aid) && isParton(aid) ) return 0;

  GenVertex* vtx = p->production_vertex();
  if ( !vtx ) return 0;

  const std::vector<GenParticle*> mothers = Rivet::particles(vtx, HepMC::parents);
  if ( mothers.empty() ) return 0;

  foreach (GenParticle* mother, mothers) {
    const int result = fromResonance(invalid, pv, mother);
    switch ( result ) {
      case 0: break;
      case 1:
        if ( mother->pdg_id() == id or isResonance(id) ) return 1;
        if ( !isExcludedFromResonance(aid) ) break;
      case 2: return 2;
    }
  }
  return 0;
}

