#include "CLHEP/Random/Randomize.h" 
//#include "CLHEP/Random/Ranlux64Engine.h" 

unsigned long  GetSeedFromTime(){
  auto now = std::chrono::system_clock::now();
  auto generator = CLHEP::HepRandom::getTheGenerator();
  auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count() % static_cast<int64_t>(1e16 * generator->flat());
  return seed;
}
