#ifndef __RAT_GeoArrayFactory__
#define __RAT_GeoArrayFactory__

#include <RAT/GeoFactory.hh>
#include <G4VPhysicalVolume.hh>

namespace RAT {

class GeoArrayFactory : public GeoFactory {
public:
  GeoArrayFactory();
  virtual ~GeoArrayFactory() {}
  virtual G4VPhysicalVolume *Construct(DBLinkPtr table);
};

} // namespace RAT

#endif
