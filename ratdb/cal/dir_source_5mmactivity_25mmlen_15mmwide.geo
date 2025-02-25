// Created:   01.08.2023
// Modified:  09.28.2023

// Downward facing directional electron source
// Exit located at origin
// Shaved conical top, fully encased in a sealant material

// 1-5mm diameter conical exit borehole
// Subgeometries from top (+z) to bottom

// Modifications:
// borehole shape changed to: 1-5mm cone, 25mm long
// back diameter decreased to: 15mm
// outer plastic shielding thickness increased to: 1mm
// fiber position moved to: 4mm from exit

{ // mother volume to be inserted in End, material will be set in macro to match "end_inner"
  name: "GEO",
  index: "source_mother",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "",
  type: "box",
  size: [50.0, 50.0, 50.0],
  position: [0.0, 0.0, 0.0], // Make the source exit at z = 0
  rotation: [0.0, 0.0, 0.0],
  material: "cryostat_vacuum",
}

{ //rear sealing, 1mm thick around shielding, z=45-25
  name: "GEO",
  index: "rear_seal",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "source_mother",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [7.5, 7.5],
  z: [25.0, 45.0],
  position: [0.0, 0.0, 0.0],
  rotation: [0.0, 0.0, 0.0],
  material: "polypropylene",
}

{ //rear outer shielding, reaches to z=44-25
  name: "GEO",
  index: "rear_outershield",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "rear_seal",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [6.5, 6.5],
  z: [25.0, 44.0]
  position: [0.0, 0.0, 0.0],
  rotation: [0.0, 0.0, 0.0],
  material: "tungsten",
}

{ //rear inner shielding, reaches to z=44-25
  name: "GEO",
  index: "rear_innershield",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "rear_outershield",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [5.0, 5.0],
  z: [25.0, 44.0],
  position: [0.0, 0.0, 0.0],
  rotation: [0.0, 0.0, 0.0],
  material: "polypropylene",
}

{ //radioactive source 0.1mm thick z=25.1-25.0
  name: "GEO",
  index: "dir_source",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "rear_innershield",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [2.5, 2.5],
  z: [25.0, 25.1]
  position: [0.0, 0.0, 0.0],
  rotation: [0.0, 0.0, 0.0],
  material: "polypropylene",
}

// "Front/Bottom portion"
{ //rear sealing, 1mm thick around shielding
  name: "GEO",
  index: "front_seal",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "source_mother",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [1.5, 7.5],
  z: [0.0, 25.0],
  position: [0.0, 0.0, 0.0], 
  material: "polypropylene",
}

{ 
  name: "GEO",
  index: "front_outershield",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "front_seal",
  type: "revolve",
  r_min: [0.5, 2.5],//1-5mm diameter cone for exit
  r_max: [0.5, 6.5],
  z: [0.1, 25.0],
  position: [0.0, 0.0, 0.0], 
  material: "tungsten",
}

{ 
  name: "GEO",
  index: "front_innershield",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "front_outershield",
  type: "revolve",
  r_min: [1.6967871485943775, 2.5],//1-5mm diameter cone for exit
  r_max: [2.5, 5.0],
  z: [15.0, 25.0],
  position: [0.0, 0.0, 0.0], 
  material: "polypropylene",
}

// Source Exit 1-5mm diameter
{
  name: "GEO",
  index: "borehole",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "front_seal",
  type: "revolve",
  r_min: [0.0, 0.0]
  r_max: [0.5, 2.5],
  z: [0.1, 25.0],
  position: [0.0, 0.0, 0.0],
  material: "cryostat_vacuum",
  color: [13.0, 255.0, 60.0, 1.0]
}

{ // scintillating fiber gate, .2mm thick
  name: "GEO",
  index: "fiber",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "borehole",
  type: "revolve",
  r_min: [0.0, 0.0]
  r_max: [0.8052208835341366, 0.821285140562249],
  z: [3.9, 4.1],
  position: [0.0, 0.0, 0.0],
  material: "polypropylene",
}
