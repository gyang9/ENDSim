#include "GeoGdmlFactory.hh"
#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <G4GDMLParser.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4Material.hh>
#include <G4VisAttributes.hh>
#include <iostream>

namespace RAT {

GeoGdmlFactory::GeoGdmlFactory() : GeoSolidFactory("gdml") {
}

G4VSolid *GeoGdmlFactory::ConstructSolid(DBLinkPtr table) {
    return nullptr;
}

//static GeoGdmlFactory gGeoGdmlFactory;

// Helper function to recursively extract aluminum volumes and their global transforms
void ExtractAluminumVolumes(G4VPhysicalVolume* physVol, 
                            G4ThreeVector globalPos, 
                            G4RotationMatrix globalRot,
                            G4LogicalVolume* targetMother,
                            const std::string& namePrefix,
                            int& copyNumber) {
    
    if (!physVol) return;
    
    G4LogicalVolume* logVol = physVol->GetLogicalVolume();
    if (!logVol) return;
    
    // Get this volume's local transform
    G4ThreeVector localPos = physVol->GetTranslation();
    G4RotationMatrix* localRot = physVol->GetRotation();
    
    // Compute global transform
    G4ThreeVector newGlobalPos = globalPos + globalRot * localPos;
    G4RotationMatrix newGlobalRot = globalRot;
    if (localRot) {
        newGlobalRot = globalRot * (*localRot);
    }
    
    // Check if this volume is aluminum
    G4Material* mat = logVol->GetMaterial();
    bool isAluminum = (mat && mat->GetName() == "G4_Al");
    
    if (isAluminum) {
        // Create a copy of this aluminum volume in the target mother
        G4RotationMatrix* newRot = new G4RotationMatrix(newGlobalRot);
        std::string newName = namePrefix + "_" + logVol->GetName() + "_" + std::to_string(copyNumber);
        
        new G4PVPlacement(newRot, newGlobalPos, logVol, newName, targetMother, false, copyNumber);
        copyNumber++;
        
        RAT::info << "GeoGdmlFactory: Extracted aluminum volume " << newName 
                  << " at position (" << newGlobalPos.x()/CLHEP::mm << ", " 
                  << newGlobalPos.y()/CLHEP::mm << ", " 
                  << newGlobalPos.z()/CLHEP::mm << ") mm" << newline;
    }
    
    // Recursively process daughter volumes
    for (int i = 0; i < logVol->GetNoDaughters(); i++) {
        G4VPhysicalVolume* daughter = logVol->GetDaughter(i);
        ExtractAluminumVolumes(daughter, newGlobalPos, newGlobalRot, targetMother, namePrefix, copyNumber);
    }
}

G4VPhysicalVolume *GeoGdmlFactory::Construct(DBLinkPtr table) {
  
  std::string filename;
  try {
    filename = table->GetS("gdml_file");
  } catch (DBNotFoundError &e) {
    Log::Die("GeoGdmlFactory: gdml_file parameter missing from table");
  }

  // Check if we should flatten the hierarchy (extract only aluminum)
  bool flattenHierarchy = false;
  try {
    flattenHierarchy = table->GetI("flatten_aluminum");
  } catch (DBNotFoundError &e) {}

  G4GDMLParser parser;
  parser.Read(filename, false);
  
  G4VPhysicalVolume* worldPhys = parser.GetWorldVolume();
  if (!worldPhys) {
    Log::Die("GeoGdmlFactory: Failed to load world volume from " + filename);
  }
  
  G4LogicalVolume* worldLog = worldPhys->GetLogicalVolume();

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

  std::string motherName = "";
  try {
      motherName = table->GetS("mother");
  } catch (DBNotFoundError &e) {}

  G4LogicalVolume* motherLog = nullptr;
  if (!motherName.empty()) {
      G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
      for (auto vol : *store) {
          if (vol->GetName() == motherName) {
              motherLog = vol;
              break;
          }
      }
      if (!motherLog) {
          Log::Die("GeoGdmlFactory: Could not find mother volume '" + motherName + "'");
      }
  }

  G4VPhysicalVolume* phys = nullptr;
  
  if (flattenHierarchy && motherLog) {
      // Extract only aluminum volumes and place them directly in mother
      RAT::info << "GeoGdmlFactory: Flattening hierarchy - extracting aluminum volumes only" << newline;
      
      int copyNumber = 0;
      G4ThreeVector initialPos = position;
      G4RotationMatrix initialRot = *rotation;
      
      ExtractAluminumVolumes(worldPhys, initialPos, initialRot, motherLog, table->GetIndex(), copyNumber);
      
      RAT::info << "GeoGdmlFactory: Extracted " << copyNumber << " aluminum volumes" << newline;
      
      // Do NOT place the GDML structure - we've already extracted what we need
      // Return nullptr since we're not placing a physical volume
      phys = nullptr;
      
  } else {
      // Normal mode: place the entire GDML structure as-is
      phys = new G4PVPlacement(rotation, position, worldLog, table->GetIndex(), motherLog, false, 0);
  }

  // Apply aluminum surface properties
  static G4OpticalSurface* alSurface = nullptr;
  if (!alSurface) {
      try {
          RAT::DBLinkPtr loptics = RAT::DB::Get()->GetLink("OPTICS_CUSTOM", "AluminumSurface");
          if (loptics) {
              alSurface = new G4OpticalSurface("AluminumSurface");
              alSurface->SetType(dielectric_metal);
              alSurface->SetModel(unified);
              alSurface->SetFinish(polished);

              G4MaterialPropertiesTable* alMPT = new G4MaterialPropertiesTable();
              
              std::string refl_opt = "";
              try { refl_opt = loptics->GetS("REFLECTIVITY_option"); } catch (...) {}
              
              std::vector<double> reflectivity = loptics->GetDArray("REFLECTIVITY_value2");
              std::vector<double> wavelengths_refl = loptics->GetDArray("REFLECTIVITY_value1");
              
              if (refl_opt == "wavelength" && wavelengths_refl.size() == reflectivity.size() && wavelengths_refl.size() > 0) {
                  std::vector<double> energies;
                  for (double w : wavelengths_refl) {
                      energies.push_back(1239.84193 * CLHEP::eV / w);
                  }
                  std::reverse(energies.begin(), energies.end());
                  std::reverse(reflectivity.begin(), reflectivity.end());
                  alMPT->AddProperty("REFLECTIVITY", energies.data(), reflectivity.data(), energies.size());
              } else {
                 try {
                     std::vector<double> energies = loptics->GetDArray("photon_energies");
                     std::vector<double> refl = loptics->GetDArray("REFLECTIVITY_value");
                     if (energies.size() == refl.size()) {
                         alMPT->AddProperty("REFLECTIVITY", energies.data(), refl.data(), energies.size());
                     }
                 } catch (...) {}
              }

              std::string eff_opt = "";
              try { eff_opt = loptics->GetS("EFFICIENCY_option"); } catch (...) {}
              std::vector<double> efficiency = loptics->GetDArray("EFFICIENCY_value2");
              std::vector<double> wavelengths_eff = loptics->GetDArray("EFFICIENCY_value1");
              
              if (eff_opt == "wavelength" && wavelengths_eff.size() == efficiency.size() && wavelengths_eff.size() > 0) {
                  std::vector<double> energies;
                  for (double w : wavelengths_eff) {
                      energies.push_back(1239.84193 * CLHEP::eV / w);
                  }
                  std::reverse(energies.begin(), energies.end());
                  std::reverse(efficiency.begin(), efficiency.end());
                  alMPT->AddProperty("EFFICIENCY", energies.data(), efficiency.data(), energies.size());
              } else {
                  try {
                     std::vector<double> energies = loptics->GetDArray("photon_energies");
                     std::vector<double> eff = loptics->GetDArray("EFFICIENCY_value");
                     if (energies.size() == eff.size()) {
                         alMPT->AddProperty("EFFICIENCY", energies.data(), eff.data(), energies.size());
                     }
                 } catch (...) {}
              }
               
              alSurface->SetMaterialPropertiesTable(alMPT);

          } else {
              throw DBNotFoundError("OPTICS", "AluminumSurface", ""); 
          }
      } catch (DBNotFoundError &e) {
          alSurface = new G4OpticalSurface("AluminumSurface");
          alSurface->SetType(dielectric_metal);
          alSurface->SetModel(unified);
          alSurface->SetFinish(polished);
          
          G4MaterialPropertiesTable* alMPT = new G4MaterialPropertiesTable();
          const G4int num = 2;
          G4double pp[num] = {2.038*CLHEP::eV, 4.144*CLHEP::eV}; 
          G4double reflectivity[num] = {0.9, 0.9};
          G4double efficiency[num] = {0.0, 0.0};
          
          alMPT->AddProperty("REFLECTIVITY", pp, reflectivity, num);
          alMPT->AddProperty("EFFICIENCY", pp, efficiency, num);
          alSurface->SetMaterialPropertiesTable(alMPT);
      } catch (...) {
      }
  }

  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  int alVolCount = 0;
  for (auto vol : *store) {
      if (vol->GetMaterial() && vol->GetMaterial()->GetName() == "G4_Al") {
          if (!G4LogicalSkinSurface::GetSurface(vol)) {
              new G4LogicalSkinSurface("Al_Skin_" + vol->GetName(), vol, alSurface);
              alVolCount++;
          }
      }
  }


  return phys;
}



}  // namespace RAT
