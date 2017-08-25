
![alt text](https://github.com/dgallichan/retroMoCoBox/blob/master/images/retroMocoBox_logo_small.png?raw=true "Retro MocoBox Logo") 

## Retro-MoCo-Box


Matlab toolbox for retrospective motion-correction of 3D MRI k-space data - as used for my work using 3D FatNavs to obtain the motion information. You can read more about FatNavs on [my research website](http://www.cardiff.ac.uk/people/view/507850-gallichan-daniel).


Used in:

* [Retrospective correction of involuntary microscopic head movement using highly accelerated fat image navigators (3D FatNavs) at 7T](http://doi.wiley.com/10.1002/mrm.25670), _Gallichan, Marques and Gruetter_, MRM 2015
* [Optimizing the acceleration and resolution of three-dimensional fat image navigators for high-resolution motion correction at 7T](http://doi.wiley.com/10.1002/mrm.26127), _Gallichan and Marques,_ MRM 2016
* [Motion-Correction Enabled Ultra-High Resolution In-Vivo 7T-MRI of the Brain](http://dx.plos.org/10.1371/journal.pone.0154974), _Federau and Gallichan,_ Plos One 2016

All of the work cited above is based on 3D FatNavs data from a Siemens 7T scanner - but in collaboration with the [Spinoza Center for Neuroimaging](https://www.spinozacentre.nl/) we have also demonstrated that the same code can be used to correct data collected with 3D EPI FatNavs from a Philips 7T scanner (_Fast and Flexible 3D-EPI Fat Navigators for High-Resolution Brain Imaging at 7 Tesla_ Buur et al, Proc ISMRM 2016). If you are also interested in using the code adapted for Philips scanners (not yet included in this Github), [please contact me](mailto:gallichand@cardiff.ac.uk).


---

### Installation

You can choose to either download the last 'release' (which hopefully won't contain any accidental commits which break it ;) ) or you can clone the folder to allow easy updates via git:

```
$ git clone https://github.com/dgallichan/retroMoCoBox.git
```

---

### Demo

`run_retroMocoDemo.m`

The simple demonstration code  uses a single volume of real data from our scanner to demonstrate how the measured motion parameters can be used to correct the 3D k-space. The translations correpond to a simple phase ramp in k-space, but the rotations move the sampling positions away from a simple Cartesian grid, so some form of gridding is necessary. This was implemented here using the [Michigan Image Reconstruction Toolbox (MIRT) for Matlab](http://web.eecs.umich.edu/~fessler/code/index.html) - which is a prerequisite for using the RetroMoCoBox. 

The data for the demo can be downloaded [here for the 1 mm dataset with large motion](http://goo.gl/ERULZA) (32 Mb), [here for the 600 um dataset (medium motion)](http://goo.gl/wto1MK) (86 Mb) and [here for the 1mm dataset with small motion](https://goo.gl/oEnLgQ) (32 Mb).

---

### Full script for MP2RAGE with 3D FatNavs

`reconstructSiemensMP2RAGEwithFatNavs.m`

This is the code we currently use for the complete reconstruction pipeline, starting from the raw data. The full script will probably only be useful if you already have the MP2RAGE with FatNavs pulse sequence. If you have a 7T Siemens system and are interested in obtaining a C2P version of the sequence - [please let me know](mailto:gallichand@cardiff.ac.uk). If you are interested in using it for 3T, it is also available for VB17 and VE11B/C - and hopefully some other baselines soon. If you don't have a 3D FatNavs sequence, perhaps some parts of the pipeline are also useful for other research. 

The only part of my code that is not included here is the FatNav-modification version of Philipp Ehses' `mapVBVD.m` code for reading the raw Siemens data into Matlab. Philipp's original code is freely available from the Siemens MR-IDEA online forum, but I have not included the FatNav-modified code here as it appears to be 'Siemens sensitive'. If you would like this code too, please [email me](mailto:gallichand@cardiff.ac.uk).

I have included various open-source tools inside the toolbox - but you will need to separately download and install both the [Michigan Image Reconstruction Toolbox (MIRT) for Matlab](http://web.eecs.umich.edu/~fessler/code/index.html) (used for the NUFFT) and [SPM 12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (used for the co-registration of the individual FatNav images to estimate the motion parameters).

The full script takes rather a long time to run (from around 10 mins through to several hours, depending on the number of CPUs you have and the size of the dataset, as well as the amount of RAM you have available) as it has to do the full reconstruction of the raw data (which typically also requires a GRAPPA reconstruction of the 3D dataset) and then perform the motion-correction step for each RF channel and each inversion time of the MP2RAGE scan. 

The script `run_SiemensMP2RAGErecon.m` gives an example of how to call the reconstruction code on your data. Also, feel free to get in touch if you would like some full example datasets to test - the raw data for 32 RF channels at 1mm resolution is 4.3 Gb and for 600 um resolution is 7.9 Gb so I haven't put them online by default, but anonymized raw data can also be shared if you are interested.

---

### Additional things to note

Be aware that because the FatNavs are never acquired coincidentally with the host data, the quality of the correction they provide will depend on the kind of motion that took place. Slow smooth motion can be corrected well, as can a few instances where there was sudden motion. However, if the subject moves continuously and in a 'random' way, then the FatNav motion estimates will not correspond to the head position when the host data were acquired. I have not noticed that this causes a problem for healthy volunteers, but it could be important for scanning certain patient populations. It also means that when testing that the FatNavs are working, it is best to perform a small slow movement throughout the scan as the 'deliberate motion' to be corrected.
