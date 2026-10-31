#include "RpcResponse.hh"
#include "RpcSpatial.hh"
#include <chrono>
#include <iostream>

int main() try {
  std::mt19937_64 rng(20261018);
  const int trials = 120000;
  auto checkMean = [&](double a, double b, double ds, double tolerance) {
    double sum = 0;
    for (int i = 0; i < trials; ++i) sum += rpc::Advance(4, a, b, ds, rng);
    const double expected = 4 * std::exp((a-b)*ds);
    rpc::Require(std::abs(sum/trials - expected) < tolerance, "Birth/death mean validation failed");
  };
  checkMean(2, 0, .3, .04); checkMean(0, 2, .3, .025);
  checkMean(2, 2, .3, .04); checkMean(1, 3, .3, .03);
  rpc::Transport t;
  t.gap_mm=2; t.voltage_v=1000; t.temperature_k=298.15; t.pressure_torr=760;
  t.velocity_mm_ns=.1; t.alpha_per_mm=0; t.eta_per_mm=0;
  auto av=rpc::Grow(100, 1., t, .03, 100000, rng);
  rpc::Require(std::abs(av.electron_induced_e-50)<1e-9 && std::abs(av.duration_ns-10)<1e-9,
               "Planar Ramo charge/time validation failed");
  // Validate integrated mean against analytic unsaturated gain, including step convergence.
  t.alpha_per_mm=2;
  double charge=0;
  for(int i=0;i<20000;++i) charge+=rpc::Grow(10,.5,t,.01,100000000,rng).electron_induced_e;
  const double expected=10*std::expm1(2*.5)/(2*t.gap_mm);
  rpc::Require(std::abs(charge/20000-expected)<.07,"Integrated avalanche charge validation failed");
  rpc::DelayCalibration dc{2,3,1,-1};
  rpc::Require(dc.x(5,2)==7 && dc.y(5,2)==8,"Delay signs/offsets validation failed");
  rpc::NominalDelayLine delay;
  const auto ends=delay.Encode(12,-9,4,0,rng);
  rpc::Require(std::abs(delay.calibration.x(ends.xl,ends.xr)-12)<1e-10 &&
               std::abs(delay.calibration.y(ends.yd,ends.yu)+9)<1e-10 &&
               std::abs(ends.xl+ends.xr-328)<1e-10,"Nominal delay round-trip failed");
  std::vector<rpc::Deposit> d(3);
  d[0].energy_ev=1; d[1].energy_ev=3; d[1].x_mm=4; d[1].t_ns=1;
  d[2].energy_ev=1; d[2].t_ns=20;
  auto hits=rpc::Empirical(d,5,0,rng);
  rpc::Require(hits.size()==2 && hits[0].x_mm==3,"Time grouping/barycenter validation failed");
  // Same-time distinct tracks far apart produce two spatial positions; nearby
  // secondaries merge, preserving ionizing energy and weighted centroid.
  std::vector<rpc::Deposit> spatial(4);
  for(int i=0;i<4;++i) { spatial[i].energy_ev=1; spatial[i].track=i+1; }
  spatial[0].x_mm=0; spatial[1].x_mm=1; spatial[1].energy_ev=3;
  spatial[2].x_mm=20; spatial[3].x_mm=40; spatial[3].t_ns=50.001;
  auto groups=rpc::SpatialReadout(spatial,50,3.745);
  rpc::Require(groups.size()==2&&groups[0].positions.size()==2&&groups[1].positions.size()==1,
    "Spatial/time separation validation failed");
  rpc::Require(groups[0].positions[0].hit.x_mm==.75&&groups[0].positions[0].hit.weight==4&&
    groups[0].positions[0].tracks==2,"Nearby-track merge failed");
  spatial[3].t_ns=50; spatial[3].layer=1;
  rpc::Require(rpc::SpatialReadout(spatial,50,3.745).size()==2,"Spatial layer isolation failed");
  spatial[0].x_mm=0; spatial[1].x_mm=3; spatial[2].x_mm=6;
  rpc::Require(rpc::SpatialReadout(spatial,50,3.745)[0].positions.size()==1,
    "Connected unresolved chain must not be split into artificial hits");
  double sx=0,sxx=0,sy=0,syy=0,sxy=0;
  const auto start=std::chrono::steady_clock::now();
  for(int i=0;i<trials;++i) {
    auto h=rpc::Empirical({d[0]},0,.749,rng)[0];
    sx+=h.x_mm; sxx+=h.x_mm*h.x_mm; sy+=h.y_mm; syy+=h.y_mm*h.y_mm; sxy+=h.x_mm*h.y_mm;
  }
  const double mx=sx/trials,my=sy/trials;
  const double sigx=std::sqrt(sxx/trials-mx*mx),sigy=std::sqrt(syy/trials-my*my);
  rpc::Require(std::abs(sigx-.749)<.008 && std::abs(sigy-.749)<.008 &&
               std::abs(sxy/trials-mx*my)<.01,"Measured spatial sigma/independence validation failed");
  t.Save("validation_transport.rpc");
  rpc::Require(rpc::Transport::Load("validation_transport.rpc").gap_mm==2,"Response round-trip failed");
  bool rejected=false;
  try { rpc::Grow(100,1,t,.01,1,rng); } catch(const std::exception&) { rejected=true; }
  rpc::Require(rejected,"Uncalibrated saturation must reject, not silently clip");
  const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::cout << "PASS birth/death, Ramo integral, delay calibration, time grouping, response round-trip, validity guard\n"
            << "Gaussian sigma X=" << sigx << " mm Y=" << sigy << " mm\n"
            << "Isolated hits=" << trials << " elapsed_s=" << elapsed << '\n';
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
