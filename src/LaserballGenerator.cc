#include <LaserballGenerator.hh>
#include <iostream>

#include <G4VSolid.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Ellipsoid.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4VisAttributes.hh>
#include <G4PVPlacement.hh>

#include <CLHEP/Units/SystemOfUnits.h>

namespace END {

LaserballGenerator::LaserballGenerator() {
    // Build a laserball geometry
    // Set configurables from laserball ratdb
    // Set sphere radius
    // Set sphere z-position
    laserballDB = RAT::DB::Get()->GetLink("LASERBALL", "end");
    // Load the led generator parameters
    std::string laserball_laser = laserballDB->GetS("led_table_index");
    RAT::Gen_LED::SetState(laserball_laser);

    // Override specifics
    this->led_z[0] = laserballDB->GetD("z_position");

    // Create the laserball geometry
    CreateGeometry();
}


void LaserballGenerator::SetState(G4String state) {

}

void LaserballGenerator::CreateGeometry() {
    const double maximum_extent = 1000.0;
    const double z_position = laserballDB->GetD("z_position");
    const std::string motherName = laserballDB->GetS("mother");
    const double sphere_radius = laserballDB->GetD("sphere_radius");
    const double deployment_rod_radius = laserballDB->GetD("deployment_rod_radius");
    const std::string deployment_rod_material = laserballDB->GetS("deployment_rod_material");
    const double holder_rod_radius = laserballDB->GetD("holder_rod_radius");
    const std::string holder_rod_material = laserballDB->GetS("holder_rod_material");
    const double holder_rod_length = laserballDB->GetD("holder_rod_length");
    const std::string sphere_material = laserballDB->GetS("sphere_material");

    G4VPhysicalVolume *motherPhys = RAT::GeoFactory::FindPhysMother(motherName);
    G4LogicalVolume *motherLog = RAT::GeoFactory::FindMother(motherName);
    G4RotationMatrix *rotation = motherPhys->GetObjectRotation();
    G4ThreeVector position = motherPhys->GetObjectTranslation();
    // Stick
    G4Material *g4_deployment_rod_material = G4Material::GetMaterial(deployment_rod_material);
    const double rod_length = maximum_extent - z_position - holder_rod_length - sphere_radius;
    G4Tubs *deployment_rod_body = new G4Tubs("deployment_rod_body", 
            0.0,
            deployment_rod_radius * CLHEP::mm,
            rod_length/2,
            0.0,
            CLHEP::twopi);
    G4LogicalVolume *deployment_rod_log = new G4LogicalVolume(
            deployment_rod_body, 
            g4_deployment_rod_material, 
            "deployment_rod_log");
    G4VPhysicalVolume *deployment_rod_phys = new G4PVPlacement(
            rotation,
            position + G4ThreeVector(0, 0, rod_length/2 + z_position + sphere_radius + holder_rod_length),
            "deployment_rod",
            deployment_rod_log,
            motherPhys,
            false,
            0);
    SetVis(deployment_rod_log, laserballDB->GetDArray("deployment_rod_color"));
    // Little stick
    G4Material *g4_holder_rod_material = G4Material::GetMaterial(holder_rod_material);
    G4Tubs *holder_rod_body = new G4Tubs("holder_rod_body", 
            0.0,
            holder_rod_radius * CLHEP::mm,
            holder_rod_length/2,
            0.0,
            CLHEP::twopi);
    G4LogicalVolume *holder_rod_log = new G4LogicalVolume(
            holder_rod_body, 
            g4_holder_rod_material, 
            "holder_rod_log");
    G4VPhysicalVolume *holder_rod_phys = new G4PVPlacement(
            rotation,
            position + G4ThreeVector(0, 0, holder_rod_length/2 + z_position + sphere_radius),
            "holder_rod",
            holder_rod_log,
            motherPhys,
            false,
            0);
    SetVis(holder_rod_log, laserballDB->GetDArray("holder_rod_color"));
    // Sphere
    G4Material *g4_sphere_material = G4Material::GetMaterial(sphere_material);
    G4Sphere *sphere_body = new G4Sphere("sphere_body", 
            //(sphere_radius - 0.1) * CLHEP::mm,
            0.0,
            sphere_radius * CLHEP::mm,
            0.0,
            CLHEP::twopi,
            0.0,
            CLHEP::pi);
    G4LogicalVolume *sphere_log = new G4LogicalVolume(
            sphere_body, 
            g4_sphere_material, 
            "sphere_log");
    G4VPhysicalVolume *sphere_phys = new G4PVPlacement(
            rotation,
            position + G4ThreeVector(0, 0, z_position),
            "sphere",
            sphere_log,
            motherPhys,
            false,
            0);
    SetVis(sphere_log, laserballDB->GetDArray("sphere_color"));
}

void LaserballGenerator::SetVis(G4LogicalVolume *volume, std::vector<double> color) {
  G4VisAttributes *att = GetColor(color);
  att->SetForceSolid(true);
  volume->SetVisAttributes(att);
}

G4VisAttributes* LaserballGenerator::GetColor(std::vector<double> color) {
  if (color.size() == 4) {
    return new G4VisAttributes(G4Color(color[0], color[1], color[2], color[3]));
  } else if (color.size() == 3) {
    return new G4VisAttributes(G4Color(color[0], color[1], color[2]));
  } else {
    return new G4VisAttributes(G4Color());
  }
}

LaserballGenerator::~LaserballGenerator() {

}

} // namespace END
