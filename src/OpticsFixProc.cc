#include <OpticsFixProc.hh>
#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4Material.hh>
#include <iostream>
#include <vector>
#include <algorithm>

namespace END {

OpticsFixProc::OpticsFixProc() : RAT::Processor("opticsfix") {
    fprintf(stderr, "OpticsFixProc: Instantiated.\n");
}

void OpticsFixProc::BeginOfRun(RAT::DS::Run *run) {
    fprintf(stderr, "OpticsFixProc: BeginOfRun called. Applying optics fix...\n");
    
    // Apply optics fix in BeginOfRun when geometry is available
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
                }
                  
                alSurface->SetMaterialPropertiesTable(alMPT);
                fprintf(stderr, "OpticsFixProc: Loaded AluminumSurface from RATDB.\n");
            } else {
                fprintf(stderr, "OpticsFixProc: AluminumSurface not found in OPTICS table.\n");
            }
        } catch (RAT::DBNotFoundError &e) {
            fprintf(stderr, "OpticsFixProc: DBNotFoundError - AluminumSurface not in DB. Using fallback.\n");
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
            fprintf(stderr, "OpticsFixProc: ERROR: Unknown exception.\n");
        }
    }

    if (alSurface) {
        G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
        int alVolCount = 0;
        for (auto vol : *store) {
            if (vol->GetMaterial() && vol->GetMaterial()->GetName() == "aluminum") {
                if (!G4LogicalSkinSurface::GetSurface(vol)) {
                    new G4LogicalSkinSurface("Al_Skin_" + vol->GetName(), vol, alSurface);
                    alVolCount++;
                }
            }
        }
        fprintf(stderr, "OpticsFixProc: Attached AluminumSurface to %d volumes.\n", alVolCount);
    }
    
    fprintf(stderr, "OpticsFixProc: BeginOfRun complete.\n");
}

RAT::Processor::Result OpticsFixProc::Event(RAT::DS::Root *ds, RAT::DS::EV *ev) {
    return RAT::Processor::OK;
}

}  // namespace END
