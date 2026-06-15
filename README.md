# ENDSim — END Neutrino Detector Simulation

Geant4 simulation of the END (Experiment for Neutrino Detection), an underwater neutrino telescope using KM3NeT-style multi-PMT digital optical modules (DOMs). Built on [ratpac-two](https://github.com/rat-pac/ratpac-two).

## Prerequisites

- ROOT 6.25+
- Geant4 11.0+
- ratpac-two (installed and sourced)

## Installation

```bash
source /path/to/ratpac-two/install/bin/ratpac.sh
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j$(nproc) && make install
```

## Quick Start

```bash
source install/bin/end.sh
end vis.mac              # Interactive visualization
end numu.mac             # Neutrino simulation (batch)
end muon.mac             # Muon simulation (batch)
```

## Detector Geometry

The detector consists of:

- **Water volume**: 50m x 50m x 500m box filled with lake water
- **Rock bed**: positioned below the detector
- **5 Detection Unit (DU) strings**: aluminum A-frame structures imported from GDML
- **90 DOMs**: each a VITROVEX 17-inch borosilicate glass sphere containing 31 PMTs

### DOM Structure

Each DOM is a nested Geant4 volume:

```
detector (lake_water)
  +-- dom_glass_N   (G4Sphere, R=216mm, borosilicate_glass)
  |     +-- dom_inner_N  (G4Sphere, R=203mm, lake_water)
  |           +-- 31 x R12199 PMTs (oriented radially outward)
  +-- gdml_du_1..5   (aluminum DU structures)
  +-- rock
```

**Glass properties (VITROVEX / Schott DURAN 3.3 borosilicate):**

| Property | Value |
|---|---|
| Outer radius | 216 mm (17 inch) |
| Wall thickness | 13 mm |
| Density | 2.23 g/cm3 |
| Refractive index | 1.473 at 587 nm |
| Transmission (400 nm, 14 mm) | ~88% |
| Composition | SiO2 80.6%, B2O3 13%, Na2O 2.1%, K2O 2%, Al2O3 2.3% |

### PMT Layout

- 90 DOMs total, arranged in a hex pattern across 4 vertical strings
- 31 Hamamatsu R12199 (3-inch) PMTs per DOM
- PMT positions defined in `ratdb/End/PMTINFO.ratdb`, table `DOMINFO_aframe_spacing5m_hex_v3`
- PMT IDs are sequential: DOM 0 has PMTs 0-30, DOM 1 has PMTs 31-61, etc.

## Running Simulations

### 1. Single Event with Visualization

```bash
source install/bin/end.sh
end vis.mac
```

Opens the Geant4 Qt viewer. Rotate, zoom, and inspect the geometry interactively.

### 2. Neutrino Events (Batch Mode)

Edit `numu.mac` to set the input vertex file and number of events:

```
/generator/add vertexfile nuInput/out.numu_100000.rootracker.0-1.randomVtx_0_100.root
/run/beamOn 100
```

Then run:

```bash
source install/bin/end.sh
end numu.mac
```

Output goes to `output.root`. To specify an output file:

```bash
end -o my_output.root numu.mac
```

### 3. Muon Events

For a single muon with configurable energy and direction, edit `muon.mac`:

```
/generator/vtx/set mu- 0.0 0.0 1.0 1000.0    # direction (x,y,z) and KE in MeV
/generator/pos/set 0.0 0.0 0.0                 # start position in mm
/run/beamOn 10
```

### 4. Cosmic Ray Muons

For realistic cosmic muon flux using the CRY generator:

```
/rat/db/set CRY latitude 37.8715
/rat/db/set CRY subboxLength 10
/generator/add combo cry:point
/generator/pos/set 0.0 0.0 2000.0
```

### 5. Parallel Batch Processing

Use `cosmic_gen.sh` or `generate_mac_files.py` to create many macro files, then run in parallel:

```bash
python3 generate_mac_files.py          # Creates macro files in a directory
bash cosmic_gen.sh                     # Runs all macros in parallel using GNU Parallel
```

Edit `cosmic_gen.sh` to set the input directory and number of cores.

## Output Files

The simulation produces ROOT files containing:

- **`T` tree** (event data):
  - `ds.mc` — Monte Carlo truth (particles, tracks, PMT hits)
  - `ds.ev` — Triggered event data (PMT times, charges)

- **`meta` tree** (run metadata):
  - `geo_file`, `experiment` — geometry configuration
  - PMT positions, source positions

### Reading Output

Using Python with `uproot`:

```python
import uproot
import awkward as ak

f = uproot.open("output.root")
tree = f["T"]

# Get MC PMT hit data
mc = tree["ds.mc"].array()
```

Using ROOT C++:

```cpp
TFile *f = new TFile("output.root");
TTree *T = (TTree*)f->Get("T");
RAT::DS::Root *ds = new RAT::DS::Root();
T->SetBranchAddress("ds", &ds);
```

### Event Display

Use `displayEvents.py` to visualize individual events:

```bash
python3 displayEvents.py
```

Edit the file paths at the top of the script to point to your output ROOT file.

## Macro Reference

All macros begin with the detector setup:

```
/rat/db/set DETECTOR experiment "End"
/rat/db/set DETECTOR geo_file "End/END.geo"
/run/initialize
```

### Processors

Add these after `/run/initialize` to control output:

| Processor | Description |
|---|---|
| `/rat/proc simpledaq` | Simple DAQ — sums charge, records earliest hit time per PMT |
| `/rat/proc count` | Print event counter (`/rat/procset update N` for every N events) |
| `/rat/proclast outroot` | Write ROOT output file |

### Generators

| Generator | Example | Description |
|---|---|---|
| `gun` | `/generator/vtx/set mu- 0 0 1 1000` | Single particle, fixed direction/energy |
| `cry` | `/generator/add combo cry:point` | Cosmic ray shower (needs CRY library) |
| `vertexfile` | `/generator/add vertexfile file.root` | Pre-generated vertices (GENIE output) |
| `pbomb` | `/generator/vtx/set 1000 400` | Photon bomb (N photons, wavelength in nm) |

## File Structure

```
ENDSim/
  ratdb/
    End/
      END.geo              # Main geometry (detector + DOMs + DU structures)
      END_dom.geo          # DOM glass sphere definitions (auto-generated)
      PMTINFO.ratdb        # PMT position tables
    MATERIALS_borosilicate.ratdb   # Glass material definition
    OPTICS_borosilicate.ratdb      # Glass optical properties
    PMT.ratdb              # PMT model definitions
  src/                     # C++ source code
  macros/                  # Example macros
  vis.mac                  # Visualization macro
  numu.mac                 # Neutrino simulation macro
  muon.mac                 # Muon simulation macro
  generate_dom_geometry.py # Regenerate DOM geometry from PMTINFO
  generate_mac_files.py    # Generate batch macro files
  cosmic_gen.sh            # Parallel batch runner
```
