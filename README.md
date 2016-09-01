
![alt text](https://github.com/dgallichan/retroMoCoBox/blob/master/images/retroMocoBox_logo_small.png?raw=true "Retro MocoBox Logo") 

## Retro-MoCo-Box


Matlab toolbox for retrospective motion-correction of 3D MRI k-space data - as used for my work using 3D FatNavs to obtain the motion information.

---

### Installation

You can choose to either clone the folder to allow easy updates via git:

```
$ git clone https://github.com/dgallichan/retroMoCoBox.git
```

Or you can download the package and update manually as required.

---

### Demo

`run_retroMocoDemo.m`

The simple demonstration code  uses a single volume of real data from our scanner to demonstrate how the measured motion parameters can be used to correct the 3D k-space. The translations correpond to a simple phase ramp in k-space, but the rotations move the sampling positions away from a simple Cartesian grid, so some form of gridding is necessary. This was implemented here using [Prof. Jeff Fessler's toolbox for Matlab](http://web.eecs.umich.edu/~fessler/code/index.html) - which is a prerequisite for using the RetroMoCoBox. 

The data for the demo can be [downloaded here]() (32 Mb).

---

### Full script for MP2RAGE with 3D FatNavs

`reconstructSiemensMP2RAGEwithFatNavs.m`

This is the code we currently use for the complete reconstruction pipeline, starting from the raw data. The full script will probably only be useful if you already have the MP2RAGE with FatNavs pulse sequence. If you have a 7T Siemens system and are interested in obtaining a C2P version of the sequence - [please let me know](daniel.gallichan@epfl.ch). If you are interested in using it for 3T, we don't yet have it ported beyond VB17, but this can be done and we hope to have it available soon. If you don't have a 3D FatNavs sequence, perhaps some parts of the pipeline are also useful for other research. 

The only part of my code that is not included here is the modified version of Philipp Ehses' `mapVBVD.m` code for reading the raw Siemens data into Matlab - which is freely available from the Siemens MR-IDEA online forum, but I have not included here as it appears to be 'Siemens sensitive'. If you would like this code too, please [email me](daniel.gallichan@epfl.ch).

I have included various open-source tools inside the toolbox - but  you will need to separately download and install both [Prof. Jeff Fessler's toolbox for Matlab](http://web.eecs.umich.edu/~fessler/code/index.html) (used for the NUFFT) and [SPM 12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (used for the co-registration of the individual FatNav images to estimate the motion parameters).

The full script takes rather a long time to run (one to several hours, depending on the number of CPUs you have and the size of the dataset) as it has to do the full reconstruction of the raw data (which typically also requires a GRAPPA reconstruction of the 3D dataset) and then perform the motion-correction step for each RF channel and each inversion time of the MP2RAGE scan. 

The script `run_SiemensMP2RAGErecon.m` gives an example of how to call the reconstruction code on your data.

