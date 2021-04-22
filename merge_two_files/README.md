# Merge_two_systems

## About this script
This python script merge two single frames in GRO format in a new one. 

## Prerequisites:
You need to install the folowing libraries:
* MDAnalysis

If you want to run the Jupyter notebook you also will need:
* nglview

## Useage
To run this script are necessary two files in GRO format with a single configuration each one. The names of this files have to be  `system1.gro` and `system2.gro`.

```bash
python3 merge_two_systems.py
```

## Ouputs
* `out.gro`
  * File in GRO format with the positions of the new system. 
