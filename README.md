# MD-Trajectory-Analysis-Bio3D

## Introduction

This repository contains documentation and example scripts for performing molecular dynamics (MD) trajectory analysis using the Bio3D package in R. Bio3D is an R package that provides tools for the analysis of biomolecular structure, sequence, and molecular dynamics trajectories.

After a successful MD simulation, 
Create a new folder named analysis within the working folder

> mkdir analysis

> cd analysis

### Trajectory correction 
> gmx trjconv -s ../md.tpr -f ../md.xtc -o md_center.xtc -center -pbc mol -ur compact

Select: protein

Select: system

collect the initial structure

> gmx trjconv -s ../md.tpr -f ../md.xtc -o md.gro -pbc mol -ur compact -center -dump 0

select: protein

select: system

> gmx trjconv -s ../md.tpr -f md_center.xtc -o md.pdb -pbc mol -ur compact -center -dump 0

select: protein

select: system
                     (or) 
> gmx editconf -f md.gro -o md.pdb

[Note: you can either convert md.gro into md.pdb using gmx editconf, or directly download the md.pdb from the trajectory using -dump 0 as above mentioned] 

## Trajectory conversion
> mdconvert md_center.xtc -t md.gro -o md.dcd

## Required files prior to do analysis using Bio3D
1) md.dcd
2) md.pdb

## Bio3D installation

> install.packages("bio3d", dependencies=TRUE)

> library(bio3d)

## Trajectory Preparation

> pdb<-read.pdb("md.pdb")

> dcd<-read.dcd("md.dcd", cell=FALSE)

Select CA atoms of the protein

> ca.inds<-atom.select(pdb,"protein",elety="CA")

Fit the trajectory to the selection of protein atoms

> trj<-fit.xyz(pdb$xyz, dcd, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)

Trim the PDB to include only the selected protein CA atoms

> protpdb<-trim.pdb(pdb,ca.inds)

> protdcd<-trim(trj, col.inds=ca.inds$xyz)

Print the trimmed PDB coordinates

> print(protpdb$xyz)

Print the trimmed trajectory coordinates

> print(protdcd)

## Trajectory Analysis

### RMSD

> rd<-rmsd(protpdb, protdcd, fit=TRUE, ncore=1)

> plot(rd, typ='l',ylab="RMSD (Å)",xlab="Frame No.")

> points(lowess(rd), typ='l', col="red", lty=1, lwd=3)

> dev.copy(jpeg,filename="rmsd.jpg")

> dev.off() 

### RMSF

> rf<-rmsf(protdcd)

> plot(rf, ylab="RMSF (Å)",xlab="Residue")

>dev.copy(jpeg,filename="rmsf.jpg") 

> dev.off()

### Principal Component Analysis (PCA)

> pc<-pca.xyz(protdcd)

> plot(pc,col=bwr.colors(nrow(dcd)))

> dev.copy(jpeg, filename="pca,jpg")

> dev.off()

![pca](https://github.com/CreedxAsif/MD-Trajectory-Analysis-Bio3D/assets/122298899/66b4f598-d38b-4af5-875e-8c6a1df8776c)


### Dynamic Cross Correlation Matrix (DCCM)

> dc<-dccm.xyz(protdcd)

> plot.dccm(dc)

<img width="850" alt="DCCM" src="https://github.com/CreedxAsif/MD-Trajectory-Analysis-Bio3D/assets/122298899/3f6b4526-05f5-4f7a-be0a-951b4ab1c627">

## Acknowledgements

The MD-Trajectory-Analysis-Bio3D was documented by Mohamed Asif. 

Feel free to contribute, report issues, or suggest enhancements by creating a [GitHub Issue](https://github.com/CreedxAsif/MD-Trajectory-Analysis-Bio3D/issues) or submitting a pull request.

Enjoy analyzing molecular dynamics trajectories with Bio3D! 😉

## Citation

1) Grant, B. J., Skjærven, L., & Yao, X. (2020). The Bio3D packages for structural bioinformatics. Protein Science, 30(1), 20–30. https://doi.org/10.1002/pro.3923
