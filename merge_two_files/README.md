# Merge_two_systems

## About this script
This python script merge two single frames in GRO format.

## Pre-requisites:
You need to install the following libraries:
* MDAnalysis

If you want to run the Jupyter notebook you also will need:
* nglview

## Inputs
In your working directory you need the following files:
* `system1.gro`
  * Single configuration in GRO format of the first system
* `system2.gro`
  * Single configuration in GRO format of the second system
* `merge_two_systems.py`
  * Python script
 
## Usage
1 . Decide across what axis you want merge the systems. Set `axis` variable to `0` for $x$, `1` for $y$ or `2` for $z$

2 . Set `gap` variable to the distance between two systems (in angstroms)

3 . Run the script:
```bash
python3 merge_two_systems.py
```

## Ouputs
* `out.gro`
  * File in GRO format with the atoms positions of the new system. 
