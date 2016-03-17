#include "TopMonteCarlo/RivetTop/interface/PseudoBoostedTop.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

struct GreaterByPt
{
  bool operator()(const Particle& a, const Particle& b) {
    return a.pt() > b.pt();
  }
};

struct GreaterByPtJet
{
  bool operator()(const Jet& a, const Jet& b) {
    return a.momentum().pt() > b.momentum().pt();
  }
};

double getDeltaR(const FourMomentum& p1, const FourMomentum& p2){
  double deta = std::fabs(p1.eta() - p2.eta());
  double dphi = std::fabs(p1.phi() - p2.phi());
  return std::sqrt(deta * deta + dphi * dphi);
}

// Determine if lepton comes from hard top decay
bool PseudoBoostedTop::isFromTop(const GenParticle* p, const Event&e){
  bool fromTop = false;
  const int pdgId = p->pdg_id();
  const GenVertex* prodVtx = p->production_vertex();
  if (prodVtx != NULL) {
    const int prodVtxId = prodVtx->barcode();

    // Find parent particle
    foreach (GenParticle* p1, Rivet::particles(e.genEvent())) {
      const GenVertex* decayVtx = p1->end_vertex();
      if (decayVtx == NULL) continue;
      const int decayVtxId = decayVtx->barcode();
      if (prodVtxId != decayVtxId) continue;
      
      // Now we have found parent particle...
      const int pdgIdParent = p1->pdg_id();
      if (pdgIdParent == 24) { // parent is a W+ --> leptonic top
	if (pdgId == -11) _topDecay = CH_ELECTRON;
	if (pdgId == -13) _topDecay = CH_MUON;
	if (pdgId == -15) _topDecay = CH_TAU;
	fromTop = true;
      }
      else if (pdgIdParent == -24) { // parent is a W- --> leptonic antitop
	if (pdgId == 11) _antitopDecay = CH_ELECTRON;
	if (pdgId == 13) _antitopDecay = CH_MUON;
	if (pdgId == 15) _antitopDecay = CH_TAU;
	fromTop = true;
      }
    }
  }
  return fromTop;  
}

void PseudoBoostedTop::project(const Event& e) {
  // Lepton : genParticle
  // (b) jets: anti-kt R=0.5 using all final-state particles excluding neutrinos
  // top jets: CA R=0.8 using all particles excluding neutrinos with 140 < mjet < 250

  _isValid = false;
  _theParticles.clear();
  _topDecay = CH_HADRON;
  _antitopDecay = CH_HADRON;
  _isParticleLep = false;
  _isGenBJet = false;
  _isGenTopJet = false;
  _passParticle = false;

  _mtt = -1.;

  // Get parton-level top / state
  Particles pGenLep;
  Particles pForJet;

  Particle refLep;
  
  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    const int pdgId = p->pdg_id();

    Particle rp(*p);

    // Get lepton
    if (std::abs(pdgId) == 11 || std::abs(pdgId) == 13 || std::abs(pdgId) == 15) {
      if (isFromTop(p,e)) refLep = rp;
    }

    // Get parton-level tops
    if (std::abs(status) > 30) continue; //not from hadronization
    if (pdgId == 6){
      _top = rp;
    }
    if (pdgId == -6){
      _antitop = rp;
    }
  }

  const FourMomentum& ttbarP4 = _top.momentum() + _antitop.momentum();
  _mtt = ttbarP4.mass();

  // Now start with particle-level quantities

  // Get particles for jet clustering
  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    Particle rp(*p);
    if (status == 1 && !rp.isNeutrino()) { //consider all final state particles except neutrinos
      pForJet.push_back(rp);
    } 
  }

  // In case the clustering alters the particle collection, clone it
  Particles pForBjet = pForJet;
  Particles pForTjet = pForJet;
    
  // Then do the AK5 clustering
  FastJets ak5Jet(FastJets::ANTIKT, _bjetR);
  ak5Jet.calc(pForBjet);

  // Then do the CA8 clustering
  FastJets ca8Jet(FastJets::CAM, _tjetR);
  ca8Jet.calc(pForTjet);

  if (refLep.momentum().pt() > _lepMinPt && std::fabs(refLep.eta()) < _lepMaxEta){

    _isParticleLep = true;
    const FourMomentum& refLepP4 = refLep.momentum();

    Jets genBjets;
    Jets genTjets;
    int nGenBjets = 0;
    int nGenTjets = 0;
  
    foreach (const Jet& jet, ak5Jet.jetsByPt(_jetMinPt)) {
      if (std::fabs(jet.eta()) > _jetMaxEta) continue;
      if (getDeltaR(jet.momentum(),refLepP4) > _halfpi) continue;
      if (getDeltaR(jet.momentum(),refLepP4) < 0.1) continue;
      genBjets.push_back(jet);
      nGenBjets += 1;
    }
    
    foreach (const Jet& jet, ca8Jet.jetsByPt(_jetMinPt)) {
      if (std::fabs(jet.eta()) > _jetMaxEta) continue;
      if (getDeltaR(jet.momentum(), refLepP4) < _halfpi) continue;
      if (jet.momentum().mass() < 140.) continue;
      if (jet.momentum().mass() > 250.) continue;
      genTjets.push_back(jet);
      nGenTjets += 1;
    }
  
    if (genTjets.size() > 1) std::sort(genTjets.begin(), genTjets.end(), GreaterByPtJet());

    if (nGenBjets >= 1) _isGenBJet = true;
    if (nGenBjets >= 1 && nGenTjets >= 1) _isGenTopJet = true;
    if (_isGenTopJet == true) _particleTop = genTjets[0];
    if (_particleTop.momentum().pT() > 400.) _passParticle = true;
  }
  _isValid = true;
}
