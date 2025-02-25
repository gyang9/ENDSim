{
  name: "GEO",
  index: "source_mother",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "end_inner",
  type: "tube",
  r_max: 1.5,
  size_z: 3.0,
  position: [0.0, 0.0, 0.0],
  material: "stainless_steel",
}

{
  name: "GEO",
  index: "active_source_element",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "source_mother",
  type: "tube",
  r_max: 1.0,
  size_z: 2.5,
  position: [0.0, 0.0, 0.0],
  material: "cryostat_vacuum",
  color: [6.0, 9.0, 240.0, 0.3],
}
