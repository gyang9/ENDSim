// Downward facing directional electron source
// Shaved conical top, fully encased in a sealant material


/// REGION 1: sealed point, lead point, vacuum borehole (1/3)

{
  name: "GEO",
  index: "source_mother",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "end_inner",
  type: "box",
  size: [80.0, 80.0, 80.0],
  position: [0.0, 0.0, 0.0], // Make the source exit at z = 0
  rotation: [0.0, 0.0, 0.0],
  material: "water",
}


{
  name: "GEO",
  index: "sealed_point",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "source_mother",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [10.1, 20.1],
  z: [0.0, 6.1],
  position: [0.0, 0.0, 0.0], // Make the source exit at z = 0
  material: "polypropylene",
  color: [13.0, 153.0, 255.0, 0.2],
}

{
  name: "GEO",
  index: "lead_point",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "sealed_point",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [10.0, 20.0],
  z: [0.0, 6.0],
  position: [0.0, 0.0, 0.1], 
  material: "lead",
  color: [13.0, 153.0, 255.0, 0.2],
}

{
  name: "GEO",
  index: "bh1",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "lead_point",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [1.5, 1.5],
  z: [-3.0, 3.0],
  position: [0.0, 0.0, 3.0],
  material: "cryostat_vacuum", 
  color: [230.0, 0.0, 0.0, 0.3],
}

{
  name: "GEO",
  index: "fiber",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "bh1",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [1.5, 1.5],
  z: [-0.1, 0.1],
  position: [0.0, 0.0, -0.9],
  material: "polypropylene", 
  color: [230.0, 0.0, 0.0, 0.3],
}



// REGION 2: (sealed middle, lead middle, borehole (2/3)

{
  name: "GEO",
  index: "sealed_middle",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "source_mother",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [20.1, 20.1],
  z: [-0.625, 0.625],
  position: [0.0, 0.0, 6.725],
  material: "polypropylene", 
  color: [13.0, 153.0, 255.0, 0.2]
}

{
  name: "GEO",
  index: "lead_middle",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "sealed_middle",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [20.0, 20.0],
  z: [-0.625, 0.625],
  position: [0.0, 0.0, 0.0],
  material: "lead", 
  color: [13.0, 153.0, 255.0, 0.2]
}

{
  name: "GEO",
  index: "bh2",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "lead_middle",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [1.5, 1.5],
  z: [-0.625, 0.625],
  position: [0.0, 0.0, 0.0],
  material: "cryostat_vacuum", 
  color: [230.0, 0.0, 0.0, 0.3]
}


// REGION 3: (sealed bottom, lead bottom, plastic case, bh3, source)

{
  name: "GEO",
  index: "sealed_bottom",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "source_mother",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [20.1, 20.1],
  z: [-16.425, 16.425],
  position: [0.0, 0.0, 23.775],
  material: "polypropylene", 
  color: [13.0, 153.0, 255.0, 0.2],
}

{
  name: "GEO",
  index: "lead_bottom",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "sealed_bottom",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [20.0, 20.0],
  z: [-16.375, 16.375],
  position: [0.0, 0.0, -0.05],
  material: "lead", 
  color: [13.0, 153.0, 255.0, 0.2],
}

{
  name: "GEO",
  index: "plastic_case",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "lead_bottom",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [14.0, 14.0],
  z: [-12.75, 12.75],
  position: [0.0, 0.0, -3.625],
  material: "polypropylene", 
  color: [0.0, 150, 0.0, 0.3],
}

{
  name: "GEO",
  index: "bh3",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "plastic_case",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [1.5, 1.5],
  z: [-6.375, 6.375],
  position: [0.0, 0.0, -6.375],
  material: "cryostat_vacuum",
  color: [230.0, 0.0, 0.0, 0.3],
}

{
  name: "GEO",
  index: "dir_source",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "bh3",
  type: "revolve",
  r_min: [0.0, 0.0],
  r_max: [1.5, 1.5],
  z: [-0.25, 0.25],
  position: [0.0, 0.0, 6.125],
  material: "cryostat_vacuum", 
  color: [6.0, 9.0, 240.0, 0.3],
}
