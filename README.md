# End Simulations

End simulations using rat-pac.

## Installation
- Install Root 6.25 and Geant4 11.0
- Install ratpac-two
- Make sure you `source ratpac.sh` inside ratpac-two repository
- `mkdir build && cd build && make install`

## Running
```
export ENDDATA=$(pwd)/data
source end.sh
end <macro.mac> <rat options>
```
Use vis.mac as an example.

