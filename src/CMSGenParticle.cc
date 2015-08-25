#include "TopMonteCarlo/RivetTop/interface/CMSGenParticle.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

bool CMSGenParticle::isFromResonance(GenParticle* p) const {
  const int aid = std::abs(p->pdg_id());
  if ( isResonance(aid) ) return true;

  if ( p->production_vertex() == 0 ) return false;
  foreach ( GenParticle* pp,  Rivet::particles(p->production_vertex(), HepMC::parents) ) {
    if ( _vetoIds.find(aid) != _vetoIds.end() ) continue;
    if ( isFromResonance(pp) ) return true;
  }

  return false;
}

void CMSGenParticle::project(const Event& e) {
  _theParticles.clear();

  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    if ( status != 1 ) continue;

    const int aid = std::abs(p->pdg_id());
    if ( _vetoIds.find(aid) != _vetoIds.end() ) continue;
    if ( _vetoIdsFromResonances.find(aid) != _vetoIdsFromResonances.end() && isFromResonance(p) ) continue;

    _theParticles.push_back(Particle(*p));
  }
}
