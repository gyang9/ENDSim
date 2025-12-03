// EndDetectorFactory recomputes the following geometry components to match
// the values specified in the End_PARAMS table

{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [30000.0, 40000.0, 30000.0], // mm, half-length
material: "water",
invisible: 1,
}

//{
//name: "GEO",
//index: "tank",
//valid_begin: [0, 0],
//valid_end: [0, 0],
//mother: "world",
//type: "box",
//size: [10900, 35400, 9200],
//position: [0.0, 0.0, 0.0],
//material: "water",
//color: [1.0, 0.0, 0.0, 0.1],
//drawstyle: "solid"
//}

//{
//name: "GEO",
//index: "detector",
//valid_begin: [0, 0],
//valid_end: [0, 0],
//mother: "tank",
//type: "box",
//size: [10700,35200,9700],
//position: [0.0, 0.0, 0.0],
//material: "wbls_3pct",
//color: [0.4, 0.4, 0.6, 0.3],
//drawstyle: "solid"
//}

{
name: "GEO",
index: "inner_surface",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "box",
//size: [10100,34600,9100],
size: [10500,35100,9600],
//size: [10000,34500,9000],
//size: [10000,30000,10000],
position: [0.0, 0.0, 0.0],
material: "wbls_3pct",
color: [0.4, 0.4, 0.6, 0.3],
drawstyle: "solid"
}

{
name: "GEO",
index: "inner_pmts",
enable: 1,
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "inner_surface",
type: "pmtarray",
pmt_model: "r7081",
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 1.000,
pos_table: "PMTINFO_25kt", //generated on the fly
//start_idx: 0, //index of first inner pmt
//end_idx: 0, //index of last inner pmt
orientation: "manual",
}

{
name: "GEO",
index: "veto_pmts",
enable: 1,
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "inner_surface",
type: "pmtarray",
pmt_model: "r7081",
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/veto",
efficiency_correction: 1.000,
pos_table: "PMTINFO_25kt", //generated on the fly
//start_idx: 0, //index of first veto pmt
//end_idx: 0, //index of last veto pmt
orientation: "manual",
}

{
name: "GEO",
index: "fiducial_volume",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "inner_surface",
type: "box",
size: [8000,33000,7500],
position: [0.0, 0.0, 0.0],
material: "wbls_5pct",
color: [0.4, 0.4, 0.6, 0.3],
drawstyle: "solid"
}
