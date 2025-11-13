#include "TGeoManager.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>
#include <string>

void display_gdml(const std::string& gdml_file_path) {
    // Check if the file exists
    if (gSystem->AccessPathName(gdml_file_path.c_str())) {
        std::cerr << "Error: GDML file not found at " << gdml_file_path << std::endl;
        return;
    }

    // Import the GDML file
    TGeoManager::Import(gdml_file_path.c_str());

    if (gGeoManager) {
        // Create a canvas to draw on
        TCanvas *c1 = new TCanvas("c1", "GDML Viewer", 800, 600);
        c1->cd();

        // Draw the top volume
        gGeoManager->GetTopVolume()->Draw();

        // Update the canvas and run the ROOT application to keep it open
        c1->Update();
        gApplication->Run();
    } else {
        std::cerr << "Error: TGeoManager not initialized after import." << std::endl;
    }
}

int main(int argc, char** argv) {
    // Initialize ROOT application
    TApplication app("GDMLViewer", &argc, argv);

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_gdml_file>" << std::endl;
        return 1;
    }

    std::string gdml_file = argv[1];
    display_gdml(gdml_file);

    return 0;
}