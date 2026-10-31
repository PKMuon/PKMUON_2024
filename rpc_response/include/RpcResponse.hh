#ifndef PKMUON_RPC_RESPONSE_HH
#define PKMUON_RPC_RESPONSE_HH

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// All runtime quantities: mm, ns, eV, elementary-charge counts; fields V/mm.
namespace rpc {
inline void Require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error(message);
}
struct Transport {
  double gap_mm = 0, voltage_v = 0, temperature_k = 0, pressure_torr = 0;
  double velocity_mm_ns = 0, alpha_per_mm = 0, eta_per_mm = 0;
  double diffusion_l_sqrt_mm = 0, diffusion_t_sqrt_mm = 0;
  std::string gas_file;
  void Validate() const {
    Require(std::isfinite(gap_mm) && gap_mm > 0 && std::isfinite(voltage_v) && voltage_v > 0,
            "Missing/invalid actual gas gap or effective gas voltage");
    Require(std::isfinite(velocity_mm_ns) && velocity_mm_ns > 0 &&
            std::isfinite(alpha_per_mm) && alpha_per_mm >= 0 &&
            std::isfinite(eta_per_mm) && eta_per_mm >= 0, "Invalid transport coefficients");
    Require(std::isfinite(temperature_k) && temperature_k > 0 &&
            std::isfinite(pressure_torr) && pressure_torr > 0, "Invalid gas conditions");
  }
  void Save(const std::string& path) const {
    Validate();
    std::ofstream f(path);
    Require(bool(f), "Cannot create response file: " + path);
    f.precision(17);
    f << "schema=1\nmodel=uniform_field_transport\nlength_unit=mm\ntime_unit=ns\n"
      << "gas_file=" << gas_file << '\n'
      << "gap_mm=" << gap_mm << "\nvoltage_v=" << voltage_v
      << "\ntemperature_k=" << temperature_k << "\npressure_torr=" << pressure_torr
      << "\nvelocity_mm_ns=" << velocity_mm_ns << "\nalpha_per_mm=" << alpha_per_mm
      << "\neta_per_mm=" << eta_per_mm << "\ndiffusion_l_sqrt_mm=" << diffusion_l_sqrt_mm
      << "\ndiffusion_t_sqrt_mm=" << diffusion_t_sqrt_mm << '\n';
    Require(bool(f), "Response write failed");
  }
  static Transport Load(const std::string& path) {
    std::ifstream f(path); Require(bool(f), "Cannot open response: " + path);
    std::map<std::string, std::string> kv;
    for (std::string line; std::getline(f, line); ) {
      const auto p = line.find('=');
      if (p != std::string::npos) kv[line.substr(0, p)] = line.substr(p + 1);
    }
    Require(kv["schema"] == "1" && kv["model"] == "uniform_field_transport" &&
            kv["length_unit"] == "mm" && kv["time_unit"] == "ns", "Unsupported response schema/units");
    Transport t; t.gas_file = kv.at("gas_file");
    auto get = [&](const char* key) { return std::stod(kv.at(key)); };
    t.gap_mm=get("gap_mm"); t.voltage_v=get("voltage_v");
    t.temperature_k=get("temperature_k"); t.pressure_torr=get("pressure_torr");
    t.velocity_mm_ns=get("velocity_mm_ns"); t.alpha_per_mm=get("alpha_per_mm");
    t.eta_per_mm=get("eta_per_mm"); t.diffusion_l_sqrt_mm=get("diffusion_l_sqrt_mm");
    t.diffusion_t_sqrt_mm=get("diffusion_t_sqrt_mm"); t.Validate(); return t;
  }
};

// Exact linear birth/death transition, rather than deterministic mean gain.
// Sum of surviving geometric families = survivors + negative-binomial draw.
inline std::uint64_t Advance(std::uint64_t n, double alpha, double eta, double ds,
                             std::mt19937_64& rng) {
  Require(alpha >= 0 && eta >= 0 && ds >= 0 && std::isfinite(alpha + eta + ds), "Invalid avalanche step");
  if (!n || ds == 0 || (alpha == 0 && eta == 0)) return n;
  const double r = alpha - eta;
  double p0 = 0, q = 0;
  if (std::abs(r) < 1.e-12 * std::max(1., alpha + eta)) {
    p0 = q = alpha * ds / (1. + alpha * ds);
  } else if (r > 0) {
    const double h = -std::expm1(-r * ds);
    const double den = r + eta * h;
    p0 = eta * h / den; q = alpha * h / den;
  } else {
    const double h = -std::expm1(r * ds);
    const double den = -r + alpha * h;
    p0 = eta * h / den; q = alpha * h / den;
  }
  p0 = std::clamp(p0, 0., 1.); q = std::clamp(q, 0., 1.);
  Require(q < 1., "Unbounded linear gain; reduce step or use a calibrated saturation model");
  const auto survivors = std::binomial_distribution<std::uint64_t>(n, 1. - p0)(rng);
  if (!survivors || q == 0) return survivors;
  const double expected = double(survivors) * q / (1. - q);
  Require(expected < 1.e15, "Linear avalanche outside supported range");
  return survivors + std::negative_binomial_distribution<std::uint64_t>(survivors, 1. - q)(rng);
}
struct Avalanche {
  double electron_induced_e = 0, duration_ns = 0;
  std::uint64_t final_electrons = 0;
};
inline Avalanche Grow(std::uint64_t n, double distance_mm, const Transport& t,
                      double step_mm, std::uint64_t validity_limit, std::mt19937_64& rng) {
  t.Validate();
  Require(std::isfinite(distance_mm) && distance_mm >= 0 && distance_mm <= t.gap_mm &&
          std::isfinite(step_mm) && step_mm > 0 && validity_limit > 0, "Invalid drift geometry/limit");
  Avalanche out;
  for (double s = 0; s < distance_mm && n; ) {
    const double ds = std::min(step_mm, distance_mm - s);
    const auto next = Advance(n, t.alpha_per_mm, t.eta_per_mm, ds, rng);
    Require(n <= validity_limit && next <= validity_limit,
            "Linear model validity limit exceeded: saturation requires calibration");
    // Planar electron-only Ramo integral. Not a segmented-electrode weighting map.
    out.electron_induced_e += 0.5 * (double(n) + double(next)) * ds / t.gap_mm;
    s += ds; out.duration_ns = s / t.velocity_mm_ns; n = next;
  }
  out.final_electrons = n; return out;
}

struct Deposit {
  int layer = 0, track = 0, parent = 0, pid = 0;
  double x_mm = 0, y_mm = 0, z_mm = 0, t_ns = 0, energy_ev = 0;
};
struct Hit {
  int layer = 0; std::size_t contributors = 0;
  double x_mm = 0, y_mm = 0, t_ns = 0, weight = 0;
};
// An explicit phenomenological coincidence window, NOT a two-hit resolution model.
// No spatial merge is inferred from the single-hit sigma.
inline std::vector<Hit> Empirical(std::vector<Deposit> deposits, double window_ns,
                                  double sigma_mm, std::mt19937_64& rng) {
  Require(std::isfinite(window_ns) && window_ns >= 0 &&
          std::isfinite(sigma_mm) && sigma_mm >= 0, "Invalid readout parameters");
  for (const auto& d : deposits)
    Require(d.energy_ev > 0 && std::isfinite(d.energy_ev + d.x_mm + d.y_mm + d.t_ns), "Invalid deposit");
  std::sort(deposits.begin(), deposits.end(), [](const Deposit& a, const Deposit& b) {
    return a.layer != b.layer ? a.layer < b.layer : a.t_ns < b.t_ns;
  });
  std::vector<Hit> result;
  for (const auto& d : deposits) {
    if (result.empty() || result.back().layer != d.layer || d.t_ns - result.back().t_ns > window_ns)
      result.push_back({d.layer, 0, d.x_mm, d.y_mm, d.t_ns, 0});
    auto& h = result.back(); const double sum = h.weight + d.energy_ev;
    h.x_mm += (d.x_mm - h.x_mm) * d.energy_ev / sum;
    h.y_mm += (d.y_mm - h.y_mm) * d.energy_ev / sum;
    h.weight = sum; ++h.contributors;
  }
  if (sigma_mm > 0) for (auto& h : result) {
    h.x_mm += std::normal_distribution<double>(0, sigma_mm)(rng);
    h.y_mm += std::normal_distribution<double>(0, sigma_mm)(rng);
  }
  return result;
}
struct DelayCalibration {
  double ax_mm_ns = 0, ay_mm_ns = 0, bx_mm = 0, by_mm = 0;
  double x(double left_ns, double right_ns) const {
    Require(std::isfinite(ax_mm_ns) && ax_mm_ns != 0, "Missing X delay calibration");
    return ax_mm_ns * (left_ns - right_ns) + bx_mm;
  }
  double y(double down_ns, double up_ns) const {
    Require(std::isfinite(ay_mm_ns) && ay_mm_ns != 0, "Missing Y delay calibration");
    return ay_mm_ns * (down_ns - up_ns) + by_mm;
  }
};
struct EndpointTimes { double xl=0,xr=0,yd=0,yu=0; };
struct NominalDelayLine {
  // 8 complete 1507-40A blocks per axis: 80 intervals, 81 readout positions.
  // Pitches/offset derived from current newrpc_readout.yaml, not scan calibration.
  double total_ns=320.;
  DelayCalibration calibration{3.457/8., 3.463/8., 1.062, 0.};
  EndpointTimes Encode(double x_mm,double y_mm,double origin_ns,double sigma_mm,
                       std::mt19937_64& rng) const {
    Require(total_ns>0 && sigma_mm>=0 && calibration.ax_mm_ns!=0 && calibration.ay_mm_ns!=0,
             "Invalid delay-line model");
    // Equal independent endpoint jitter is ONE realization of measured spatial sigma.
    // Common timing jitter and real two-pulse performance remain uncalibrated.
    const double sx=sigma_mm/(std::sqrt(2.)*std::abs(calibration.ax_mm_ns));
    const double sy=sigma_mm/(std::sqrt(2.)*std::abs(calibration.ay_mm_ns));
    auto jitter=[&](double s) { return s>0?std::normal_distribution<double>(0,s)(rng):0.; };
    const double dx=(x_mm-calibration.bx_mm)/(2*calibration.ax_mm_ns);
    const double dy=(y_mm-calibration.by_mm)/(2*calibration.ay_mm_ns);
    return {origin_ns+total_ns/2+dx+jitter(sx),origin_ns+total_ns/2-dx+jitter(sx),
            origin_ns+total_ns/2+dy+jitter(sy),origin_ns+total_ns/2-dy+jitter(sy)};
  }
};
}  // namespace rpc
#endif
