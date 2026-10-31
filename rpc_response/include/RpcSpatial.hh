#ifndef PKMUON_RPC_SPATIAL_HH
#define PKMUON_RPC_SPATIAL_HH
#include "RpcResponse.hh"
#include <set>

namespace rpc {
struct SpatialPosition {
  Hit hit;
  int tracks=0;
};
struct SpatialWindow {
  int layer=0, tracks=0;
  double t_ns=0;
  std::vector<SpatialPosition> positions;
};
// Ideal spatial-response approximation authorized by the user, NOT waveform fitting.
// Build one ionizing-deposit centroid per simulated track in each fixed time group.
// Connected centroids at distance <= separation_mm form one unresolved component.
// Track identities are simulation information; no claim of measured pulse pairing.
inline std::vector<SpatialWindow> SpatialReadout(std::vector<Deposit> deposits,
                                               double window_ns,double separation_mm) {
  Require(std::isfinite(window_ns)&&window_ns>=0 &&
          std::isfinite(separation_mm)&&separation_mm>0,"Invalid spatial readout parameters");
  for(const auto& d:deposits)
    Require(d.energy_ev>0&&std::isfinite(d.energy_ev+d.x_mm+d.y_mm+d.t_ns),"Invalid spatial deposit");
  std::sort(deposits.begin(),deposits.end(),[](const Deposit& a,const Deposit& b) {
    if(a.layer!=b.layer) return a.layer<b.layer;
    if(a.t_ns!=b.t_ns) return a.t_ns<b.t_ns;
    return a.track<b.track;
  });
  std::vector<SpatialWindow> result;
  for(std::size_t first=0;first<deposits.size();) {
    const auto& start=deposits[first]; std::size_t last=first;
    std::map<int,Hit> tracks;
    while(last<deposits.size()&&deposits[last].layer==start.layer&&
          deposits[last].t_ns-start.t_ns<=window_ns) {
      const auto& d=deposits[last++]; auto& h=tracks[d.track];
      const double sum=h.weight+d.energy_ev;
      h.x_mm+=(d.x_mm-h.x_mm)*d.energy_ev/sum;
      h.y_mm+=(d.y_mm-h.y_mm)*d.energy_ev/sum;
      if(h.contributors==0) h.t_ns=d.t_ns;
      h.layer=d.layer; h.weight=sum; ++h.contributors;
    }
    std::vector<Hit> centers; for(const auto& item:tracks) centers.push_back(item.second);
    std::vector<std::size_t> component(centers.size());
    for(std::size_t i=0;i<component.size();++i) component[i]=i;
    auto root=[&](std::size_t i) { while(component[i]!=i) i=component[i]; return i; };
    for(std::size_t i=0;i<centers.size();++i) for(std::size_t j=0;j<i;++j)
      if(std::hypot(centers[i].x_mm-centers[j].x_mm,centers[i].y_mm-centers[j].y_mm)<=separation_mm)
        component[root(i)]=root(j);
    SpatialWindow group; group.layer=start.layer; group.t_ns=start.t_ns;
    group.tracks=static_cast<int>(centers.size()); std::map<std::size_t,SpatialPosition> merged;
    for(std::size_t i=0;i<centers.size();++i) {
      auto& p=merged[root(i)]; const auto& h=centers[i]; const double sum=p.hit.weight+h.weight;
      p.hit.x_mm+=(h.x_mm-p.hit.x_mm)*h.weight/sum;
      p.hit.y_mm+=(h.y_mm-p.hit.y_mm)*h.weight/sum;
      if(p.tracks==0) p.hit.t_ns=h.t_ns; else p.hit.t_ns=std::min(p.hit.t_ns,h.t_ns);
      p.hit.layer=start.layer; p.hit.weight=sum; p.hit.contributors+=h.contributors; ++p.tracks;
    }
    for(const auto& item:merged) group.positions.push_back(item.second);
    // Require emitted centers to remain mutually separated as well. Averaging a
    // long connected component can move its centroid toward another component.
    bool changed=true;
    while(changed) {
      changed=false;
      for(std::size_t i=0;i<group.positions.size()&&!changed;++i)
        for(std::size_t j=0;j<i;++j) {
          auto& a=group.positions[j]; const auto& b=group.positions[i];
          if(std::hypot(a.hit.x_mm-b.hit.x_mm,a.hit.y_mm-b.hit.y_mm)>separation_mm) continue;
          const double sum=a.hit.weight+b.hit.weight;
          a.hit.x_mm+=(b.hit.x_mm-a.hit.x_mm)*b.hit.weight/sum;
          a.hit.y_mm+=(b.hit.y_mm-a.hit.y_mm)*b.hit.weight/sum;
          a.hit.weight=sum; a.hit.contributors+=b.hit.contributors; a.tracks+=b.tracks;
          a.hit.t_ns=std::min(a.hit.t_ns,b.hit.t_ns);
          group.positions.erase(group.positions.begin()+i); changed=true; break;
        }
    }
    result.push_back(group); first=last;
  }
  return result;
}
}
#endif
