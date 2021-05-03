# conformers

## About this script
This python script calculates the conformer populations in a molecular dynamics simulations according to the following ranges of angles:
* t for 120°-240°
* g for 0°-120°
* g′ for 240°-360°.

## Pre-requisites:
You need to install the following libraries:
* MDAnalysis

## Inputs
In your working directory you need the following files:

* `file.trr`
  * Trajectory
* `file.gro`
  * Topology
* `myconfig.py`
  * Parameters to carry out the calculation
* `conformers.py`
  * Python script
  
## Usage
1 . Set the variables in `myconfig.py`:
 * `inicio`: first frame to read from trajectory
 * `fin`: last frame to read from trajectory
 * `num_diedros`: number of dihedros 
 * Now you have to create lists with the  indexes of the atoms of dihedrals of interest. You will need to substract  1 to the indexes of the topology. The names of the lists have to be diedro$n$ where *n* is the number of diedral.
 
Example:
I need calculate the populations of conformers from frame 40 to frame 100. I am interested in dihedrals whose indexes (in topology) are 1,2,3,4 and 2,3,4,5. My myconfig.py should look like:
```bash
inicio = 40          
fin =   100        
num_diedros = 2      
diedro1 = [0,1,2,3]
diedro2 = [1,2,3,4]
```

3 . Run the script:
```bash
python3 conformers.py
```

## Ouputs
 * The script prints the populations of conformers 
