TensorFit
=========
Simple training routine to create models for tensorfit. Ideally models should
be stored in `models/fitters/tensorfit/<geometry>/<material>/` for each target
material and geometry. Run with `python train.py -h` for complete information.
This program expects the data to be in an ntuple such that uproot can read it.

```bash
python utilities/tensorfit/train.py trainingfile.ntuple.root -m direction -i 1000 -d models/fitters/tensorfit/end/water/direction
```
