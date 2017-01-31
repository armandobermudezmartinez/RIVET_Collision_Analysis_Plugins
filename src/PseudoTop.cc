#include "TopMonteCarlo/RivetTop/interface/PseudoTop.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

struct GreaterByPt
{
  bool operator()(const Particle& a, const Particle& b) {
    return a.pt() > b.pt();
  }
};

void PseudoTop::cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup) const
{
  std::vector<std::map<double, std::pair<size_t, size_t> >::const_iterator> toErase;
  std::set<size_t> usedLeg1, usedLeg2;
  if ( !doCrossCleanup ) {
    for (auto key = v.begin(); key != v.end(); ++key) {
      const size_t leg1 = key->second.first;
      const size_t leg2 = key->second.second;
      if (usedLeg1.find(leg1) == usedLeg1.end() and
          usedLeg2.find(leg2) == usedLeg2.end()) {
        usedLeg1.insert(leg1);
        usedLeg2.insert(leg2);
      } else {
        toErase.push_back(key);
      }
    }
  }
  else {
    for (auto key = v.begin(); key != v.end(); ++key) {
      const size_t leg1 = key->second.first;
      const size_t leg2 = key->second.second;
      if (usedLeg1.find(leg1) == usedLeg1.end() and
          usedLeg1.find(leg2) == usedLeg1.end()) {
        usedLeg1.insert(leg1);
        usedLeg1.insert(leg2);
      } else {
        toErase.push_back(key);
      }
    }
  }
  for (auto& key : toErase) v.erase(key);
}

void PseudoTop::project(const Event& event) {
  // Leptons : do the lepton clustering anti-kt R=0.1 using stable photons and leptons not from hadron decay
  // Neutrinos : neutrinos not from hadron decay
  // MET : vector sum of all invisible particles in x-y plane 
  // Jets : anti-kt clustering using all particles excluding neutrinos and particles used in lepton clustering
  //        add ghost B hadrons during the jet clustering to identify B jets.

  // W->lv : dressed lepton and neutrino pairs
  // W->jj : light flavored dijet
  // W candidate : select lv or jj pairs which minimise |mW1-80.4|+|mW2-80.4|
  //               lepton-neutrino pair will be selected with higher priority

  // t->Wb : W candidate + b jet
  // t candidate : select Wb pairs which minimise |mtop1-172.5|+|mtop2-172.5|

  _isValid = false;
  _theParticles.clear();
  _wDecays1.clear();
  _wDecays2.clear();
  _jets.clear();
  _bjets.clear();
  _ljets.clear();
  _mode1 = _mode2 = CH_HADRON;
  
  // Get analysis objects from projections
  
  const vector<DressedLepton> leptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();

  Cut jet_cut = (Cuts::abseta < _jetMaxEta) and (Cuts::pT > _jetMinPt*GeV);
  _jets = apply<FastJets>(event, "Jets").jetsByPt(jet_cut);
  for (const Jet& jet : _jets) {
    if (jet.bTagged()) _bjets.push_back(jet);
    else               _ljets.push_back(jet);
  }
  
  const Particles neutrinos = apply<IdentifiedFinalState>(event, "Neutrinos").particlesByPt();
  
  _met = -apply<MissingMomentum>(event, "MET").vectorEt();

  // All building blocks are ready. Continue to pseudo-W and pseudo-top combination

  if (_bjets.size() < 2) return; // Ignore single top for now
  std::map<double, std::pair<size_t, size_t> > wLepCandIdxs;
  std::map<double, std::pair<size_t, size_t> > wHadCandIdxs;

  // Collect leptonic-decaying W's
  for (size_t iLep = 0, nLep = leptons.size(); iLep < nLep; ++iLep) {
    const DressedLepton& lep = leptons.at(iLep);
    for (size_t iNu = 0, nNu = neutrinos.size(); iNu < nNu; ++iNu) {
      const Particle& nu = neutrinos.at(iNu);
      if (nu.fromHadron())
        continue;
      const double m = (lep.momentum()+nu.momentum()).mass();
      const double dm = std::abs(m-_wMass);
      wLepCandIdxs[dm] = make_pair(iLep, iNu);
    }
  }

  // Continue to hadronic decaying W's
  for (size_t i = 0, nLjet = _ljets.size(); i < nLjet; ++i) {
    const Jet& ljet1 = _ljets[i];
    for (size_t j = i+1; j < nLjet; ++j) {
      const Jet& ljet2 = _ljets[j];
      const double m = (ljet1.momentum()+ljet2.momentum()).mass();
      const double dm = std::abs(m-_wMass);
      wHadCandIdxs[dm] = make_pair(i, j);
    }
  }

  // Cleanup W candidate, choose pairs with minimum dm if they share decay products
  cleanup(wLepCandIdxs);
  cleanup(wHadCandIdxs, true);
  const size_t nWLepCand = wLepCandIdxs.size();
  const size_t nWHadCand = wHadCandIdxs.size();

  if (nWLepCand + nWHadCand < 2) return; // We skip single top

  int w1Q = 1, w2Q = -1;
  int w1dau1Id = 1, w2dau1Id = -1;
  FourMomentum w1dau1LVec, w1dau2LVec;
  FourMomentum w2dau1LVec, w2dau2LVec;
  if (nWLepCand == 0) { // Full hadronic case
    const auto& idPair1 = wHadCandIdxs.begin()->second;
    const auto& idPair2 = std::next(wHadCandIdxs.begin())->second;
    const auto& w1dau1 = _ljets[idPair1.first];
    const auto& w1dau2 = _ljets[idPair1.second];
    const auto& w2dau1 = _ljets[idPair2.first];
    const auto& w2dau2 = _ljets[idPair2.second];

    w1dau1LVec = w1dau1.momentum();
    w1dau2LVec = w1dau2.momentum();
    w2dau1LVec = w2dau1.momentum();
    w2dau2LVec = w2dau2.momentum();
  } else if (nWLepCand == 1) { // Semi-leptonic case
    const auto& idPair1 = wLepCandIdxs.begin()->second;
    const auto& idPair2 = wHadCandIdxs.begin()->second;
    const auto& w1dau1 = leptons[idPair1.first];
    const auto& w1dau2 = neutrinos[idPair1.second];
    const auto& w2dau1 = _ljets[idPair2.first];
    const auto& w2dau2 = _ljets[idPair2.second];

    w1dau1LVec = w1dau1.momentum();
    w1dau2LVec = w1dau2.momentum();
    w2dau1LVec = w2dau1.momentum();
    w2dau2LVec = w2dau2.momentum();
    w1dau1Id = leptons[idPair1.first].pid();
    w1Q = w1dau1Id > 0 ? -1 : 1;
    w2Q = -w1Q;

    switch (w1dau1Id) {
      case 13: case -13: _mode1 = CH_MUON; break;
      case 11: case -11: _mode1 = CH_ELECTRON; break;
    }
  } else { // Full leptonic case
    const auto& idPair1 = wLepCandIdxs.begin()->second;
    const auto& idPair2 = std::next(wLepCandIdxs.begin())->second;
    const auto& w1dau1 = leptons[idPair1.first];
    const auto& w1dau2 = neutrinos[idPair1.second];
    const auto& w2dau1 = leptons[idPair2.first];
    const auto& w2dau2 = neutrinos[idPair2.second];

    w1dau1LVec = w1dau1.momentum();
    w1dau2LVec = w1dau2.momentum();
    w2dau1LVec = w2dau1.momentum();
    w2dau2LVec = w2dau2.momentum();
    w1dau1Id = leptons[idPair1.first].pid();
    w2dau1Id = leptons[idPair2.first].pid();
    w1Q = w1dau1Id > 0 ? -1 : 1;
    w2Q = w2dau1Id > 0 ? -1 : 1;

    switch (w1dau1Id) {
      case 13: case -13: _mode1 = CH_MUON; break;
      case 11: case -11: _mode1 = CH_ELECTRON; break;
    }
    switch (w2dau1Id) {
      case 13: case -13: _mode2 = CH_MUON; break;
      case 11: case -11: _mode2 = CH_ELECTRON; break;
    }
  }
  const auto w1LVec = w1dau1LVec+w1dau2LVec;
  const auto w2LVec = w2dau1LVec+w2dau2LVec;

  // Combine b jets
  double sumDm = 1e9;
  FourMomentum b1LVec, b2LVec;
  for (size_t i = 0, n = _bjets.size(); i < n; ++i) {
    const Jet& bjet1 = _bjets[i];
    const double mtop1 = (w1LVec+bjet1.momentum()).mass();
    const double dmtop1 = std::abs(mtop1-_tMass);
    for (size_t j=0; j<n; ++j) {
      if (i == j) continue;
      const Jet& bjet2 = _bjets[j];
      const double mtop2 = (w2LVec+bjet2.momentum()).mass();
      const double dmtop2 = std::abs(mtop2-_tMass);

      if (sumDm <= dmtop1+dmtop2) continue;

      sumDm = dmtop1+dmtop2;
      b1LVec = bjet1.momentum();
      b2LVec = bjet2.momentum();
    }
  }
  if (sumDm >= 1e9) return; // Failed to make top, but this should not happen.

  const auto t1LVec = w1LVec + b1LVec;
  const auto t2LVec = w2LVec + b2LVec;

  // Put all of them into candidate collection
  _t1 = Particle(w1Q*6, t1LVec);
  _b1 = Particle(w1Q*5, b1LVec);
  _w1 = Particle(w1Q*24, w1LVec);
  _wDecays1.push_back(Particle(w1dau1Id, w1dau1LVec));
  _wDecays1.push_back(Particle(-w1dau1Id+w1Q, w1dau2LVec));

  _t2 = Particle(w2Q*6, t2LVec);
  _b2 = Particle(w2Q*5, b2LVec);
  _w2 = Particle(w2Q*24, w2LVec);
  _wDecays2.push_back(Particle(w2dau1Id, w2dau1LVec));
  _wDecays2.push_back(Particle(-w2dau1Id+w2Q, w2dau2LVec));

  _isValid = true;
}

