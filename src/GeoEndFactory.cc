#include <CLHEP/Units/SystemOfUnits.h>

#include <G4Ellipsoid.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>
#include <G4SubtractionSolid.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <GeoEndFactory.hh>
#include <RAT/GeoTubeFactory.hh>

namespace END {

G4VSolid *GeoEndFactory::ConstructSolid(RAT::DBLinkPtr table) {
  // Inner radius of cylindrical part of End
  G4double r_min = table->GetD("r_min");
  // Outer radius of cylindrical part of End
  G4double r_max = table->GetD("r_max");
  // Half-height of cylindrical part of End
  G4double size_z = table->GetD("size_z");
  // Radius of ellipical bottom/top cap
  G4double top_radius = table->GetD("top_radius");
  // Height of ellipical bottom/top cap
  G4double top_height = table->GetD("top_height");
  // Offset for cap
  G4double offset = table->GetD("offset");

  // Solids for the cylindrical body and ellipical caps
  G4Tubs *body = new G4Tubs("body", r_min * CLHEP::mm, r_max * CLHEP::mm, size_z * CLHEP::mm, 0., CLHEP::twopi);

  G4Ellipsoid *head = new G4Ellipsoid("head", top_radius * CLHEP::mm, top_radius * CLHEP::mm, top_height * CLHEP::mm,
                                      0., top_height * CLHEP::mm);

  G4Ellipsoid *bot = new G4Ellipsoid("bot", top_radius * CLHEP::mm, top_radius * CLHEP::mm, top_height * CLHEP::mm,
                                     -top_height * CLHEP::mm, 0.0);

  // Location and rotation of the top cap
  G4ThreeVector *trans = new G4ThreeVector(0., 0., (size_z - offset) * CLHEP::mm);
  G4RotationMatrix *rotation = new G4RotationMatrix();
  G4Transform3D *transf = new G4Transform3D(*rotation, *trans);

  // Add the top cap to the cylinder for both the vessel and cavity
  G4VSolid *EndVolumeP1 = new G4UnionSolid("end_v1", body, head, *transf);

  // Location of the bottom cap (same rotation as the top)
  G4ThreeVector *neg_trans = new G4ThreeVector(0., 0., -(size_z - offset) * CLHEP::mm);
  G4Transform3D *neg_transf = new G4Transform3D(*rotation, *neg_trans);

  // Add the bottom cap to the (top cap + cylinder) for both the vessel and cavity
  G4UnionSolid *EndVolume = new G4UnionSolid("end", EndVolumeP1, bot, *neg_transf);

  return EndVolume;
}

}  // namespace END
