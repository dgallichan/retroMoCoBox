
![alt text](https://github.com/dgallichan/retroMoCoBox/blob/master/images/retroMocoBox_logo_small.png?raw=true "Retro MocoBox Logo") 

## Retro-MoCo-Box


Matlab toolbox for retrospective motion-correction of 3D MRI k-space data - as used for my work using 3D FatNavs to obtain the motion information. You can read more about FatNavs on [my research website](http://www.cardiff.ac.uk/people/view/507850-gallichan-daniel).


Used in:

* [Retrospective correction of involuntary microscopic head movement using highly accelerated fat image navigators (3D FatNavs) at 7T](http://doi.wiley.com/10.1002/mrm.25670), _Gallichan, Marques and Gruetter_, MRM 2015
* [Optimizing the acceleration and resolution of three-dimensional fat image navigators for high-resolution motion correction at 7T](http://doi.wiley.com/10.1002/mrm.26127), _Gallichan and Marques,_ MRM 2016
* [Motion-Correction Enabled Ultra-High Resolution In-Vivo 7T-MRI of the Brain](http://dx.plos.org/10.1371/journal.pone.0154974), _Federau and Gallichan,_ Plos One 2016
* [Evaluation of 3D fat-navigator based retrospective motion correction in the clinical setting of patients with brain tumors](https://doi.org/10.1007/s00234-019-02160-w), _Glessgen, Gallichan, Moor, Hainc and Federau_, Neuroradiology 2019
* [Fat navigators and Moiré phase tracking comparison for motion estimation and retrospective correction](https://doi.org/10.1002/mrm.27908), _Gretsch, Mattern, Gallichan, Speck_, MRM 2020
* [Sharpness in motion corrected quantitative imaging at 7T](https://doi.org/10.1016/j.neuroimage.2020.117227), _Bazin, Nijsse, van der Zwaag, Gallichan, Alkemade, Vos, Forstmann, Caan_, NeuroImage 2020  


Most of the work cited above (all except the Bazin 2020 paper!) is based on data collected on Siemens scanners. If you are also interested in using the code adapted for Philips scanners (not yet included in this Github), [please contact me](mailto:gallichand@cardiff.ac.uk).


---

### Installation

You can choose to either download the last 'release' or you can download the latest version packaged as a zip file (which hopefully won't contain any accidental commits which break it ;) ). You can also clone the folder to allow easy updates via git:

```
$ git clone https://github.com/dgallichan/retroMoCoBox.git
```

---

### Requirements

The only part of my code that is not included here is the FatNav-modification version of Philipp Ehses' `mapVBVD.m` code for reading the raw Siemens data into Matlab. Philipp's original code is freely available from the Siemens MR-IDEA online forum, but I have not included the FatNav-modified code here as it appears to be 'Siemens sensitive'. If you would like this code too, please [email me](mailto:gallichand@cardiff.ac.uk).

I have included various open-source tools inside the toolbox - but you will need to separately download and install both the [Michigan Image Reconstruction Toolbox (MIRT) for Matlab](http://web.eecs.umich.edu/~fessler/code/index.html) (used for the NUFFT) and [SPM 12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (used for the co-registration of the individual FatNav images to estimate the motion parameters).

---

### Demos

#### `run_retroMocoDemo_NUFFT.m`

The simple demonstration code  uses a single volume of real data from our scanner to demonstrate how the measured motion parameters can be used to correct the 3D k-space. The translations correpond to a simple phase ramp in k-space, but the rotations move the sampling positions away from a simple Cartesian grid, so some form of gridding is necessary. This was implemented here using the [Michigan Image Reconstruction Toolbox (MIRT) for Matlab](http://web.eecs.umich.edu/~fessler/code/index.html) - the necessary part of which is now included here in the 'mirt_nufft' folder.

The data for the demo can be downloaded [here for the 1 mm dataset with large motion](http://goo.gl/ERULZA) (32 MB), [here for the 600 um dataset (medium motion)](http://goo.gl/wto1MK) (86 MB) and [here for the 1mm dataset with small motion](https://goo.gl/oEnLgQ) (32 MB).

#### `run_retroMocoDemo_simulateMotion.m`

This script starts with any example volume (I found the [original Colin27 brain](http://www.bic.mni.mcgill.ca/ServicesAtlases/Colin27) to be handy for this) and then simulates the motion artifacts caused by different motion profiles. By then attempting to apply motion-correction to this simulated data it is possible to explore how well retrospective correction can be expected to work in different motion regimes. In the case of very 'rough' motion the retrospective correction still has noticeable artifacts - which can be almost completely suppressed by using an iterative NUFFT operation instead of the single-step approach used by default. 

#### `run_retroMocoDemo_FatNavRecon.m`

This script loads an example of real 3D FatNavs data (15 volumes at 2mm resolution with 4x4 GRAPPA acceleration) and does the GRAPPA reconstruction. The example data [can be downloaded from here](https://goo.gl/1qYjsc) (186 MB).


---

### Full script for MP2RAGE with 3D FatNavs

`reconstructSiemensMP2RAGEwithFatNavs.m`

This is the code we currently use for the complete reconstruction pipeline, starting from the raw data. The full script will probably only be useful if you already have the MP2RAGE with FatNavs pulse sequence. If you have a 7T Siemens system and are interested in obtaining a C2P version of the sequence - [please let me know](mailto:gallichand@cardiff.ac.uk). If you are interested in using it for 3T, it is also available for VB17 and VE11B/C - and hopefully some other baselines soon. If you don't have a 3D FatNavs sequence, perhaps some parts of the pipeline are also useful for other research. 

The full script takes rather a long time to run (from around 10 mins through to several hours, depending on the number of CPUs you have and the size of the dataset, as well as the amount of RAM you have available) as it has to do the full reconstruction of the raw data (which typically also requires a GRAPPA reconstruction of the 3D dataset) and then perform the motion-correction step for each RF channel and each inversion time of the MP2RAGE scan. 

The script `run_SiemensMP2RAGErecon.m` gives an example of how to call the reconstruction code on your data. Also, feel free to get in touch if you would like some full example datasets to test - the raw data for 32 RF channels at 1mm resolution is 4.3 GB and for 600 um resolution is 9.3 GB so I haven't put them online by default, but anonymized raw data can also be shared if you are interested.

---

### Additional things to note

## "Correctable" motion 
Be aware that because the FatNavs are never acquired coincidentally with the host data, the quality of the correction they provide will depend on the kind of motion that took place. Slow smooth motion can be corrected well, as can a few instances where there was sudden motion. However, if the subject moves continuously and in a 'random' way, then the FatNav motion estimates will not correspond to the head position when the host data were acquired. I have not noticed that this causes a problem for healthy volunteers, but it could be important for scanning certain patient populations. It also means that when testing that the FatNavs are working, it is best to perform a small slow movement throughout the scan as the 'deliberate motion' to be corrected.

## Phase or slice oversampling
Currently my code for the Siemens reconstruction will not properly handle the data if phase oversampling or slice oversampling is on (it runs but gets the voxel dimension wrong and so the correction won't be quite right either). This is not actually particularly complicated to fix - but would require a fair amount of time to implement and debug. I therefore recommend increasing the FOV instead of oversampling. Please let me know if anyone would *require* oversampling and we can look at changing the code to handle it properly.
