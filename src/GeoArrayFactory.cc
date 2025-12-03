#include "GeoArrayFactory.hh"
#include <RAT/DB.hh>
#include <RAT/GeoBuilder.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <iostream>

namespace RAT {

GeoArrayFactory::GeoArrayFactory() : GeoFactory("geoarray") {
  std::cout << "DEBUG: GeoArrayFactory instantiated and registered." << std::endl;
}

G4VPhysicalVolume *GeoArrayFactory::Construct(DBLinkPtr table) {
  std::string volume_name = table->GetS("volume"); // The name of the logical volume to place
  std::string index = table->GetIndex();
  
  // Find the logical volume
  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* logical = nullptr;
  
  // Search for the volume by name. 
  // Note: RAT might prefix names or use specific naming conventions.
  // Usually the logical volume name matches the index in the GEO table if created by standard factories.
  for (auto vol : *store) {
      if (vol->GetName() == volume_name) {
          logical = vol;
          break;
      }
  }

  if (!logical) {
      // Try to find it via GeoBuilder if it hasn't been built yet?
      // But GeoBuilder builds in order. If 'volume' is listed before this array in dependency, it should be built.
      // However, RATDB doesn't guarantee order unless 'mother' is set.
      // If the template volume has no mother (or world), it might be built.
      // If we can't find it, we can't place it.
      std::cerr << "GeoArrayFactory Error: Could not find logical volume '" << volume_name << "'" << std::endl;
      return nullptr;
  }

  // Get position and rotation
  std::vector<double> pos = table->GetDArray("position");
  std::vector<double> rot = table->GetDArray("rotation");
  
  G4ThreeVector position(0, 0, 0);
  if (pos.size() == 3) {
      position.set(pos[0], pos[1], pos[2]);
  }
  
  G4RotationMatrix* rotation = new G4RotationMatrix();
  if (rot.size() == 3) {
      rotation->rotateX(rot[0] * CLHEP::deg);
      rotation->rotateY(rot[1] * CLHEP::deg);
      rotation->rotateZ(rot[2] * CLHEP::deg);
  }

  // Create the placement
  // Note: The mother volume is handled by GeoBuilder, which calls this factory.
  // This factory returns the physical volume, and GeoBuilder adds it to the mother.
  // BUT, G4PVPlacement constructor *requires* the mother logical volume if we want to place it immediately.
  // GeoFactory::Construct signature doesn't provide the mother.
  // Standard RAT factories return a G4VPhysicalVolume*. 
  // Wait, if we return a G4PVPlacement, we must provide a mother?
  // Or can we pass nullptr as mother and let GeoBuilder handle it?
  // Let's check how GeoTubeFactory does it.
  // Usually: new G4PVPlacement(rotation, position, logical, name, mother_logical, false, 0);
  
  // If GeoBuilder handles the mother, it must be taking the returned volume and... wait.
  // G4VPhysicalVolume cannot be reparented easily.
  // RAT::GeoBuilder::Construct() calls factory->Construct(table).
  // Inside GeoBuilder, it finds the mother volume.
  // Does it pass the mother to the factory? No.
  
  // Let's look at GeoEndFactory.cc again.
  // It returns G4VSolid* ?? No, it returns G4VPhysicalVolume* in the header?
  // Wait, I viewed GeoEndFactory.cc and it had:
  // G4VSolid *GeoEndFactory::ConstructSolid(RAT::DBLinkPtr table)
  // It didn't have Construct().
  
  // Ah! GeoEndFactory inherits from GeoFactory but implements ConstructSolid?
  // Standard GeoFactory has Construct().
  // Maybe GeoEndFactory inherits from something else?
  // Or maybe it implements Construct() and calls ConstructSolid()?
  
  // If I look at RAT::GeoTubeFactory (standard), it likely implements Construct().
  // And Construct() usually does:
  // 1. ConstructSolid()
  // 2. Create Logical Volume
  // 3. Return new G4PVPlacement(..., mother=0) ??
  
  // If mother is 0, it's a world volume?
  // But these are daughters.
  
  // Actually, in Geant4, if you pass 0 as mother, it's not placed in the hierarchy yet.
  // But G4PVPlacement *requires* a mother logical volume unless it's the very top world.
  
  // Let's check RAT::GeoFactory header if possible.
  // I can't.
  
  // Let's assume I can pass nullptr as mother.
  // But wait, if I pass nullptr, it's a floating volume.
  // Does RAT GeoBuilder add it?
  // If I look at `GeoEndFactory.cc` again...
  // It has `ConstructSolid`. It does NOT have `Construct`.
  // This implies `GeoEndFactory` might be a `GeoSolidFactory`?
  // Or `End` class registers it differently?
  
  // Let's check `End.cc`.
  // `new GeoEndFactory();`
  
  // If `GeoEndFactory` only implements `ConstructSolid`, then it must be used by a generic factory that handles the placement?
  // But `END.geo` has `type: "box"` for detector.
  // `GeoEndFactory` is likely NOT used for "box".
  // `GeoEndFactory` is used for... what?
  // `END.geo` doesn't use `type: "end"`.
  
  // Wait, `END.geo` has:
  // `type: "box"` for world.
  // `type: "box"` for detector.
  // `type: "box"` for rock.
  // `type: "pmtarray"` for km3net_doms.
  
  // So `GeoEndFactory` is NOT USED in the current `END.geo`!
  // It might be legacy code or for a specific `type: "end"` that isn't currently used.
  
  // So I should look at how standard factories work.
  // Since I can't see RAT source, I have to guess or use `GeoGdmlFactory` as reference.
  // `GeoGdmlFactory` implements `Construct`.
  // `G4VPhysicalVolume *GeoGdmlFactory::Construct(DBLinkPtr table)`
  // It returns a `G4VPhysicalVolume*`.
  
  // In `GeoGdmlFactory.cc`:
  // return new G4PVPlacement(0, G4ThreeVector(0,0,0), world_log, filename, 0, false, 0);
  // It passes `0` (nullptr) as mother!
  // So yes, I can pass nullptr as mother.
  
  return new G4PVPlacement(rotation, position, logical, index, nullptr, false, 0);
}

} // namespace RAT
