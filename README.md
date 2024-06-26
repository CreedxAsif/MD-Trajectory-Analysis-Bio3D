# MD-Trajectory-Analysis-Bio3D

After successful MD simulation, 
Create a new folder named analysis within the working folder

> mkdir analysis
> 
> cd analysis

### Trajectory correction 
> gmx trjconv -s ../md.tpr -f ../md.xtc -o md_center.xtc -center -pbc -ur mol -ur compact

Select: protein

Select: system

collect the initial structure
> gmx trjconv -s ../md.tpr -f md_center.xtc -o start.gro -dump 0
> 
> gmx trjconv -s ../md.tpr -f md_center.xtc -o start.pdb -dump 0 
                     (or) 
> gmx editconf -f start.gro -o start.pdb

[Note: you can either convert start.gro into start.pdb using gmx editconf, or directly download the start.pdb from the trajectory using -dump 0 as above mentioned] 

## Trajectory conversion
> mdconvert md_center.xtc -t md.gro -o md.dcd

## Required files prior to do analysis using Bio3D
1) md.dcd
2) Start.pdb

## Bio3D installation

> install.packages("bio3d", dependencies=TRUE)

> library(bio3d)

## Trajectory Preparation

> pdb<-read.pdb("start.pdb")

> dcd<-read.dcd("md.dcd", cell=FALSE)

Select CA atoms of the protein

> ca.inds<-atom.select(pdb,"protein",elety="CA")

Select all atoms of the ligand (assuming the ligand residue name is 'LIG', if not change it accordingly)

> ligand.inds<-atom.select(pdb, "ligand", rename="LIG")

Combine protein CA and ligand atom indices

> combined.inds<-c(ca.inds$xyz, ligand.inds$xyz)

Fit the trajectory to the combined selection of protein and ligand atoms

> trj<-fit.xyz(pdb$xyz, dcd, fixed.inds=combined.inds, mobile.inds=combined.inds)

Trim the PDB to include only the selected protein CA and ligand atoms

> protlig_pdb<-trim.pdb(pdb, atom.select(pdb, "combined", elety="CA", resname="LIG"))

> protlig_dcd<-trim(trj, col.inds=combined.inds)

Print the trimmed PDB coordinates

> print(protlig_pdb$xyz)

Print the trimmed trajectory coordinates

> print(protlig_dcd)

## Trajectory Analysis

### RMSD

> rd<-rmsd(protlig_pdb, protlig_dcd, fit=TRUE, ncore=1)

> plot.rmsd(rd)

### RMSF

> rf<-rmsf(protlig_dcd)

> plot.rmsf(rf)

### component Analysis (PCA)

> pc<-pca.xyz(protlig_dcd)

> plot.pca(pc)

### Dynamic Cross Correlation Matrix (DCCM)

> dc<-dccm.xyz(protlig_dcd)

> plot.dccm(dc)
