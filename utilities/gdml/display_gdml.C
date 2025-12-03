#include <TGeoManager.h>
#include <TBrowser.h>

void display_gdml() {
    // Create a new TGeoManager. The name is not critical.
    new TGeoManager("world", "The detector geometry");

    // Import the GDML file. ROOT automatically detects the file type.
    // Make sure "detector.gdml" is in the same directory where you run this script.
    gGeoManager->Import("detector.gdml");

    // Close the geometry definition. This is important for performance.
    gGeoManager->CloseGeometry();

    // Get the top volume (the world volume).
    TGeoVolume *top = gGeoManager->GetTopVolume();

    // Draw the geometry. "ogl" specifies the OpenGL viewer.
    if (top) {
        top->Draw("ogl");
    } else {
        printf("Error: Could not get the top volume. Check the GDML file.\n");
    }

    // You can also open a TBrowser to inspect the geometry hierarchy.
    // new TBrowser();
}
