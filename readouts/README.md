# Detector performance data files

## Folders
There are 4 folders for each readout technology performance curve as a function of recoil energy (keVnr)

* Energy resolution [number from 0 to 1, to be multiplied by E_r]
* Efficiency [number from 0 to 1]
* Angular resolution [in degrees]
* Head tail efficiency [number from 0 to 1]

## Format
The format of each file are three columns with:

`[Recoil energy (keV),		Value for Fluorine,	Value for Helium]`

All have 1000 rows for energies between 0 and 200 keV.

## Readouts 
Are (in the order the code defines them)
* 'Ideal' = Idealised detector perfect in everything
* 'Pixel' = pixel grid e.g. CCDs
* 'Strip' = electronic strip based readout e.g. MIMAC
* 'Optical' = e.g. fibres
* 'Wire' = e.g. in DRIFT
* 'Pad' = don't remember what this is
* 'Planar' = nor this..
* 'Nondirectional' = our toy non-directional detector (currently it's idealised)
