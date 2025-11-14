#ifndef __END_MyDAQProc_hh__
#define __END_MyDAQProc_hh__

#include <RAT/Processor.hh>

namespace END {

class MyDAQProc : public RAT::Processor {
public:
    MyDAQProc();
    virtual ~MyDAQProc();
    virtual RAT::Processor::Result DSEvent(RAT::DS::Root *ds);

private:
    int fEventCounter;
};

} // namespace END

#endif // __END_MyDAQProc_hh__
