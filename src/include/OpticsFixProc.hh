#ifndef __END_OpticsFixProc__
#define __END_OpticsFixProc__

#include <RAT/Processor.hh>

namespace END {

class OpticsFixProc : public RAT::Processor {
 public:
  OpticsFixProc();
  virtual ~OpticsFixProc() {}
  virtual void BeginOfRun(RAT::DS::Run *run);
  virtual RAT::Processor::Result Event(RAT::DS::Root *ds, RAT::DS::EV *ev);
};

}  // namespace END

#endif
