#ifndef __RAT_GeoGdmlFactory__
#define __RAT_GeoGdmlFactory__

#include <RAT/GeoSolidFactory.hh>

namespace RAT {

class GeoGdmlFactory : public GeoSolidFactory {
 public:
  GeoGdmlFactory();
  virtual ~GeoGdmlFactory() {};
  virtual G4VPhysicalVolume *Construct(DBLinkPtr table);
  virtual G4VSolid *ConstructSolid(DBLinkPtr table);
};

}  // namespace RAT

#endif // __RAT_GeoGdmlFactory__
