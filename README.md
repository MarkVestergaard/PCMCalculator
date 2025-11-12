# PCMCalculator
GUI for analyzing 2D phase contrast mapping (PCM) MRI data to measure trough-plane flow. 

Written for python 3.8. 
Tested on PCM data from 3T Philips dSTREAM Achieva MRI and 3T Siemens Biograph mMR hybrid PET/MR system. 


Input for Philips scanner data: **.PAR file <n>  
Input for Siemens scanner data: nifti-file converted from dicom using dcmniix (https://github.com/rordenlab/dcm2niix) 

Region of interest (ROI) is manual delineated or automatically drawn based using region-growin algoritm.
Data saved as csv file. 
  
  
For in-depth description of analysis see: <n> 
Vestergaard et al. Cerebral Cortex, Volume 32, Issue 6, 15 March 2022, 1295–1306, doi:https://doi.org/10.1093/cercor/bhab294 <n> or <n>
Vestergaard et al.  Journal of Cerebral Blood Flow & Metabolism 2019, Vol. 39(5) 834–848, doi:https://doi.org/10.1177/0271678X17737909


 Mark B. Vestergaard <n>  
 Functional Imaging Unit, <n>  
 Department of Clinical Physiology and Nuclear Medicine <n>  
 Rigshospitalet <n> 
 Copenhagen, Denmark <n>  
 mark.bitsch.vestergaard@regionh.dk
 

<img width="1554" height="985" alt="PCMCalculator" src="https://github.com/user-attachments/assets/3b7aac6a-69c5-423c-af3f-320fd987d1de" />
