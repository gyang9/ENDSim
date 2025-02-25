#ifndef __RAT_GeoEndFactory__
#define __RAT_GeoEndFactory__

#include <RAT/GeoSolidFactory.hh>

namespace END {
class GeoEndFactory : public RAT::GeoSolidFactory {
 public:
  GeoEndFactory() : GeoSolidFactory("end"){};
  virtual G4VSolid *ConstructSolid(RAT::DBLinkPtr table);
};

}  // namespace END

#endif
