#ifndef __END_End__
#define __END_End__

#include <Config.hh>
#include <RAT/Rat.hh>
#include <RAT/AnyParse.hh>
#include <RAT/ProcBlockManager.hh>
#include <RAT/ProcAllocator.hh>
#include <RAT/GLG4Gen.hh>
#include <RAT/Factory.hh>
#include <GeoEndFactory.hh>
//#include <DichroiconArrayFactory.hh>
#include <HitmanProc.hh>
#include <NtupleProc.hh>
#include <LaserballGenerator.hh>
#include <string>

namespace END {

class End : public RAT::Rat {
 public:
  End(RAT::AnyParse* parser, int argc, char** argv);
};

}  // namespace END

#endif
