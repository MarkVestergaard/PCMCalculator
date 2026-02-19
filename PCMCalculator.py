#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCM Flow Analysis GUI
=====================

A graphical user interface for analyzing 2D phase contrast mapping (PCM) MRI data 
to measure through-plane blood flow in cerebral vessels.

Description
-----------
This application provides tools for:
    - Loading and displaying PCM MRI data from Philips PAR/REC or NIfTI formats
    - Manual and automatic ROI delineation for vessel segmentation
    - Blood flow quantification (ml/min) from velocity-encoded MRI
    - Pulsatility analysis including PI (Pulsatility Index) and ΔV (Volume of 
      Arterial Pulsatility)
    - Data export to CSV with optional ROI masks

Supported Input Formats
-----------------------
    - Philips PAR/REC files
    - NIfTI files (.nii): Converted from DICOM using dcm2niix 
      (https://github.com/rordenlab/dcm2niix). Requires accompanying JSON sidecar.

System Requirements
-------------------
    - Python 3.13
    - Dependencies: tkinter, matplotlib, nibabel, numpy, pandas, scipy, 
      scikit-image, PIL

Tested Scanners
---------------
    - 3T Philips dSTREAM Achieva MRI
    - 3T Siemens Biograph mMR hybrid PET/MR system

Usage
-----
Command line:
    $ python PCM_GUI_ver10.py
    
With pre-specified files:
    $ python PCM_GUI_ver10.py --img <path_to_PAR_file>
    $ python PCM_GUI_ver10.py --img_nii_vel <vel.nii> --img_nii_mod <mod.nii> --img_nii_mag <mag.nii>

Keyboard Shortcuts
------------------
    Ctrl+1-5    : Switch image display (velocity, modulus, magnitude, ROI, mean velocity)
    Ctrl+Q/W/E  : Change colormap (jet, gray, viridis)
    Ctrl+N      : Load new PAR/REC file
    Ctrl+M      : Load new NIfTI file
    Ctrl+R      : Load saved ROI file
    Ctrl+S      : Save flow data
    Ctrl+P      : Calculate pulsatility
    Up/Down     : Navigate through frames

Pulsatility Measures
--------------------
    PI (Pulsatility Index):
        PI = (Qmax - Qmin) / Qmean
        - Dimensionless measure of flow waveform shape
        - Reflects downstream vascular resistance
        
    ΔV (Volume of Arterial Pulsatility):
        ΔV = max(∫[Q-Q̄]dt) - min(∫[Q-Q̄]dt)
        - Volume (mL) of cyclic arterial distension per heartbeat
        - Represents downstream arterial buffering capacity
        - Requires heart rate input for correct temporal scaling
        - When combined with information on pressure the arterial compliance can be calculated

Output Files
------------
    - CSV file: Flow (ml/min), Velocity (cm/s), Cross-sectional area (mm²) per frame
    - NPZ file (optional): ROI polygon coordinates and binary masks
    - NIfTI file (optional): ROI mask in original image space
    - GIF file (optional): Animation of ROI across cardiac cycle


Author
------
    Mark B. Vestergaard
    Functional Imaging Unit
    Department of Clinical Physiology, Nuclear Medicine and PET
    Rigshospitalet Copenhagen, Denmark
    mark.bitsch.vestergaard@regionh.dk

"""

# =============================================================================
# IMPORTS
# =============================================================================
import tkinter as tk
from tkinter import filedialog, messagebox, END, Menu, LabelFrame, Label, Entry, Button, Scale, Frame, Grid, Checkbutton
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import nibabel as nib
import argparse
import numpy as np
import numpy.matlib
import pandas as pd
from matplotlib.widgets import PolygonSelector
from matplotlib.patches import Polygon 
from matplotlib.path import Path
from skimage import measure
import os 
import scipy
from functools import partial
from PIL import Image 
import glob
from skimage.morphology import binary_dilation, disk

# =============================================================================
# COMMAND LINE ARGUMENT PARSING
# =============================================================================
# Command line arguments for batch processing or scripted usage
parser = argparse.ArgumentParser(
    description='Calculate blood flow in vessels from Phase Contrast MRI data',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog='''
Examples:
  %(prog)s                                    # Interactive file selection
  %(prog)s --img scan.PAR                     # Load Philips PAR/REC file
    '''
)
parser.add_argument('--img', dest='img_input', 
                    help='Input PCM image (PAR/REC format)')
parser.add_argument('--img_nii_vel', dest='img_input_nii_vel', 
                    help='Input velocity NIfTI image (phase data)')
parser.add_argument('--img_nii_mod', dest='img_input_nii_mod', 
                    help='Input modulus NIfTI image')
parser.add_argument('--img_nii_mag', dest='img_input_nii_mag', 
                    help='Input magnitude NIfTI image')

input_files = parser.parse_args()

# =============================================================================
# FILE LOADING
# =============================================================================
# Flag to track if working with NIFTI files (affects processing pipeline)
nifti_files = None

# Load PCM data - either from command line argument or via file dialog
# Load PCM data - either from command line argument or via file dialog
if input_files.img_input_nii_vel is not None:
    # NIfTI files provided via command line
    nifti_files = True
    import json
    raw_img_filename = input_files.img_input_nii_vel
    raw_img_filename_mod = input_files.img_input_nii_mod
    raw_img_filename_mag = input_files.img_input_nii_mag

elif input_files.img_input is not None:
    # PAR/REC file provided via command line
    raw_img_filename = input_files.img_input

else:
    # No command line arguments - open file dialog
    print('Select PCM file')
    raw_img_filename_tmp = filedialog.askopenfilename(
        initialdir="/", title="Select file",
        filetypes=(("Philips PAR", "*.PAR"), ("NIFTI", "*.nii"), ("all files", "*.*")),
        multiple=True
    )
    raw_img_filename = raw_img_filename_tmp[0]
    file_type = os.path.splitext(raw_img_filename)[1]
    if file_type == '.nii':
        nifti_files = True
        import json


# =============================================================================
# GUI WINDOW INITIALIZATION
# =============================================================================
# Create the main tkinter GUI window with responsive sizing
window = tk.Tk()
window.title('PCMCalculator')
GUI_width = window.winfo_screenwidth()*0.9
GUI_height = window.winfo_screenheight()*0.9
window.geometry( str(int(GUI_width)) +'x'+ str(int(GUI_height)) )
window.resizable(True,True)

# DPI scaling factor for different screen resolutions
# Scales UI elements relative to a 1080p baseline for consistent appearance
dpi_scale = max(GUI_width / 1920, GUI_height / 1080)
fig_dpi = 100
fig_size_main = max(3.5, 4.0 * dpi_scale)    # Main image figure size
fig_size_flow = max(1.8, 2.0 * dpi_scale)    # Flow plot figure size
slider_length = int(GUI_height * 0.4)         # Slider length proportional to screen height
button_width = max(10, int(10 * dpi_scale))   # Button width scaled
button_height = max(2, int(2 * dpi_scale))    # Button height scaled


# =============================================================================
# GLOBAL STATE VARIABLES
# =============================================================================
# Module-level state variables for image display
Disp_image_str = 'vel'      # Current image type: 'vel', 'mod', 'mag', 'ROI', 'Mean_vel'
colormap_str = 'jet'        # Current colormap: 'jet', 'gray', 'viridis'
imgFrame = 1                # Current frame index (1-based)


# =============================================================================
# IMAGE DATA CLASSES
# =============================================================================

class Img_data_parrec():
    """
    Class to handle Philips PAR/REC format PCM image data.
    
    This class loads and processes phase contrast MRI data from Philips scanners,
    extracting velocity, modulus, and magnitude images from the multi-volume PAR/REC
    format.
    
    Attributes
    ----------
    raw_img : nibabel.parrec.PARRECImage
        Raw loaded image object from nibabel
    raw_img_filename : str
        Path to the source PAR file
    img : nibabel image
        Canonical orientation image
    nifti_file : nibabel.Nifti1Image
        NIfTI representation for compatibility
    Venc : float
        Velocity encoding value (cm/s) - extracted from header or set to default
    Venc_from_hdr : bool
        Whether VENC was successfully read from header
    Vel_image : ndarray
        3D velocity image array (x, y, cardiac_phase) in cm/s
    Mod_image : ndarray
        3D modulus image array (anatomical reference)
    Thres_image : ndarray
        3D phase magnitude image (for vessel visualization)
    Mean_vel_image : ndarray
        Mean velocity across specified frames (for background subtraction)
    Inv_image_indx : int
        Multiplier for velocity direction (+1 or -1)
    
    Methods
    -------
    set_new_data()
        Load a new PAR/REC file interactively
    """
    def __init__(self, raw_img_filename):
        self.raw_img = nib.load(raw_img_filename)# Loads PCM par file nibabel load
        self.raw_img_filename=raw_img_filename
        self.img = nib.as_closest_canonical(self.raw_img) # Loads PCM par file
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        file_type = os.path.splitext(raw_img_filename)[1]
        self.Inv_image_indx=1
        
        if file_type=='.PAR':
            self.Venc=max(self.img.header.general_info.get('phase_enc_velocity')) # Find venc in header
            self.Venc_from_hdr=True
        else:
            self.Venc=100
            self.Venc_from_hdr=False

        self.thres_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==4)
        self.mod_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==2)
        self.phase_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==3,self.img.header.image_defs['scanning sequence']==4)

        self.img.header.image_defs['recon resolution']
        self.img.dataobj_sq=np.squeeze(self.img.dataobj)

        no_dyn=self.img.header.image_defs['index in REC file'][-1]+1
        self.img_data_indx=np.where(self.img.dataobj_sq.shape==no_dyn)
        self.img_data_indx=np.where(self.img.dataobj_sq.shape==no_dyn)

        if self.img_data_indx[0][0]==0:
            self.Thres_image=(np.rot90(self.img.dataobj_sq[self.thres_image_idx,:,:])/self.img.dataobj_sq[self.thres_image_idx,:,:].max()*self.Venc*0.7)
            self.Mod_image=(np.rot90(self.img.dataobj_sq[self.mod_image_idx,:,:])/self.img.dataobj_sq[self.mod_image_idx,:,:].max()*self.Venc*0.9)
            self.Vel_image=np.rot90(self.img.dataobj_sq[self.phase_image_idx,:,:])
        
        if self.img_data_indx[0][0]==1:
            self.Thres_image=(np.rot90(self.img.dataobj_sq[:,self.thres_image_idx,:])/self.img.dataobj_sq[:,self.thres_image_idx,:].max()*self.Venc*0.7)
            self.Mod_image=(np.rot90(self.img.dataobj_sq[:,self.mod_image_idx,:])/self.img.dataobj_sq[:,self.mod_image_idx,:].max()*self.Venc*0.9)
            self.Vel_image=np.rot90(self.img.dataobj_sq[:,self.phase_image_idx,:])

        if self.img_data_indx[0][0]==2:
            self.Thres_image=(np.rot90(self.img.dataobj_sq[:,:,self.thres_image_idx])/self.img.dataobj_sq[:,:,self.thres_image_idx].max()*self.Venc*0.7)
            self.Mod_image=(np.rot90(self.img.dataobj_sq[:,:,self.mod_image_idx])/self.img.dataobj_sq[:,:,self.mod_image_idx].max()*self.Venc*0.9)
            self.Vel_image=np.rot90(self.img.dataobj_sq[:,:,self.phase_image_idx])


    def set_new_data(self):
        """
        Load a new PAR/REC file via file dialog.
        
        Opens a file selection dialog for the user to choose a new PAR file.
        Reloads all image data and resets ROI masks if the image dimensions change.
        Updates the GUI display and resets the frame slider.
        """
         
        self.raw_img_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("PAR files","*.PAR"),("all files","*.*")))        
        self.raw_img = nib.load(self.raw_img_filename)
        self.img = nib.as_closest_canonical(self.raw_img) # Loads PCM par file
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        self.Venc=self.img.header.general_info.get('phase_enc_velocity')[2]
        
        self.thres_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==4)
        self.mod_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==2)
        self.phase_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==3,self.img.header.image_defs['scanning sequence']==4)

        if self.img.header.general_info.get('prep_direction')=='Right-Left':
            self.Thres_image=(np.rot90(self.img.dataobj[:,:,0,self.thres_image_idx])/self.img.dataobj[:,:,0,self.phase_image_idx].max()*self.Venc*0.2)
            self.Mod_image=(np.rot90(self.img.dataobj[:,:,0,self.mod_image_idx])/self.img.dataobj[:,:,0,self.mod_image_idx].max()*self.Venc*0.5)
            self.Vel_image=np.rot90(self.img.dataobj[:,:,0,self.phase_image_idx])
        elif self.img.header.general_info.get('prep_direction')=='Anterior-Posterior':
            self.Thres_image=(np.rot90(np.transpose(self.img.dataobj[0,:,:,self.thres_image_idx],(1,2,0)))/self.img.dataobj[0,:,:,self.phase_image_idx].max()*self.Venc*0.2)
            self.Mod_image=(np.rot90(np.transpose(self.img.dataobj[0,:,:,self.mod_image_idx],(1,2,0))) /self.img.dataobj[0,:,:,self.mod_image_idx].max()*self.Venc*0.5)
            self.Vel_image=np.rot90(np.transpose(self.img.dataobj[0,:,:,self.phase_image_idx],(1,2,0)))

        global Disp_image_str
        global colormap_str
        global imgFrame
        imgFrame=1
        slider_scale.set(1)
        slider_scale['to'] = 1
        slider_scale['from'] = self.Vel_image.shape[2]
        #slider_scale.ax.set_xlim(slider_scale.valmin,slider_scale.valmax)
        save_str=[Img_data.raw_img_filename[0:-4]+'_Flow_data.csv']

        entry_save_filename.delete(0, 'end')
        entry_save_filename.insert(END, save_str[0])
        change_image(Disp_image_str,colormap_str)
        global nifti_files
        nifti_files=None
        

        if ROI_art.BWMask.shape!=self.Vel_image.shape:
            ROI_art.polygon=[None]*Img_data.Vel_image.shape[2]
            ROI_art.BWMask=Img_data.Vel_image*False
            ROI_art.flag=[0]*(Img_data.Vel_image.shape[2])


        

# =============================================================================
# NIFTI IMAGE DATA CLASS
# =============================================================================

class Img_data_nii():
    """
    Class to handle NIfTI format PCM image data (e.g., from Siemens scanners).
    
    This class loads and processes phase contrast MRI data from NIfTI files,
    typically converted from DICOM using dcm2niix. Requires three separate NIfTI
    files: phase (velocity), modulus, and magnitude images, each with accompanying
    JSON sidecar files containing metadata.
    
    Attributes
    ----------
    raw_img : nibabel.Nifti1Image
        Primary loaded NIfTI image
    raw_img_filename : str
        Path to the primary NIfTI file
    img : nibabel image
        Canonical orientation image
    nifti_file : nibabel.Nifti1Image
        NIfTI image object
    img_b, img_c : nibabel image
        Second and third NIfTI images (different contrasts)
    json_header : dict
        Metadata from JSON sidecar file
    ImageType : list
        List of image types ['P', 'M', 'MAG'] indicating phase, modulus, magnitude
    Venc : float
        Velocity encoding value (cm/s)
    phase_scale : float
        Scaling factor for phase data normalization
    Vel_image : ndarray
        3D velocity image array (x, y, cardiac_phase) in cm/s
    Mod_image : ndarray
        3D modulus image array
    Thres_image : ndarray
        3D phase magnitude image
    
    Methods
    -------
    set_new_data()
        Load new NIfTI files interactively
        
    Notes
    -----
    The phase_scale is automatically determined based on the data range:
        - 12-bit data (max ~4096): Siemens typical
        - 16-bit signed (max ~32768): Alternative encoding
        - Normalized float (max ~1.0): Pre-scaled data
    """
    def __init__(self, raw_img_filename):
        self.raw_img = nib.load(raw_img_filename[0])
        self.raw_img_filename = raw_img_filename[0]
        self.img = nib.as_closest_canonical(self.raw_img)
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        with open(os.path.splitext(raw_img_filename[0])[0]+'.json') as f:
            self.img.json_header = json.load(f)
        ImageType_a = self.img.json_header['ImageType'][2]
        
        self.raw_img_b = nib.load(raw_img_filename[1])
        self.raw_img_filename_b = raw_img_filename[1]
        self.img_b = nib.as_closest_canonical(self.raw_img_b)
        self.nifti_file_b = nib.Nifti1Image(self.img_b.dataobj, self.img_b.affine, header=self.img_b.header)
        with open(os.path.splitext(raw_img_filename[1])[0]+'.json') as f:
            self.img_b.json_header = json.load(f)
        ImageType_b = self.img_b.json_header['ImageType'][2]
        
        self.raw_img_c = nib.load(raw_img_filename[2])
        self.raw_img_filename_c = raw_img_filename[2]
        self.img_c = nib.as_closest_canonical(self.raw_img_c)
        self.nifti_file_c = nib.Nifti1Image(self.img_c.dataobj, self.img_c.affine, header=self.img_c.header)
        with open(os.path.splitext(raw_img_filename[2])[0]+'.json') as f:
            self.img_c.json_header = json.load(f)
        ImageType_c = self.img_c.json_header['ImageType'][2]
        
        self.ImageType=[ImageType_a, ImageType_b, ImageType_c]
        self.Inv_image_indx=1
        file_type = os.path.splitext(raw_img_filename[0])[1]
        if file_type=='.PAR':
            self.Venc=self.img.header.general_info.get('phase_enc_velocity')[2] # Find venc in header
            self.Venc_from_hdr=True
        else:
            self.Venc=100
            self.Venc_from_hdr=False
        
        # Determine phase scale factor from the phase image data
        n_tmp=self.ImageType.index('P')
        if n_tmp==0:
            phase_data = np.squeeze(self.nifti_file.dataobj)
            phase_nifti = self.nifti_file
        elif n_tmp==1:
            phase_data = np.squeeze(self.nifti_file_b.dataobj)
            phase_nifti = self.nifti_file_b
        elif n_tmp==2:
            phase_data = np.squeeze(self.nifti_file_c.dataobj)
            phase_nifti = self.nifti_file_c
        
        # Calculate phase_scale from data characteristics
        # Method 1: Use max absolute value (works for symmetric data around 0)
        # Method 2: Use data type bit depth
        # Method 3: Check JSON for scaling info
        phase_max_abs = max(abs(phase_data.max()), abs(phase_data.min()))
        
        # Check if data appears to be 12-bit (common for Siemens)
        if phase_max_abs > 2048 and phase_max_abs <= 4096:
            self.phase_scale = 4096
        # Check if data appears to be 16-bit signed
        elif phase_max_abs > 4096 and phase_max_abs <= 32768:
            self.phase_scale = 32768
        # Check if data is already normalized (float between -1 and 1)
        elif phase_max_abs <= 1.0:
            self.phase_scale = 1.0
        # Otherwise use the actual max value for normalization
        else:
            self.phase_scale = phase_max_abs
        
        print(f"Phase scale factor determined: {self.phase_scale} (data max: {phase_max_abs})")
            
        
        n_tmp=self.ImageType.index('M')
        if n_tmp==0:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file.dataobj/self.nifti_file.dataobj.max())*self.Venc*0.7, k=-1) 
        elif n_tmp==1:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file_b.dataobj/self.nifti_file_b.dataobj.max())*self.Venc*0.7, k=-1) 
        elif n_tmp==2:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file_c.dataobj/self.nifti_file_c.dataobj.max())*self.Venc*0.7,k=-1)   
            
        n_tmp=self.ImageType.index('MAG')
        if n_tmp==0:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file.dataobj/self.nifti_file.dataobj.max())*self.Venc*0.9 ,k=-1) 
        elif n_tmp==1:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file_b.dataobj/self.nifti_file_b.dataobj.max())*self.Venc*0.9 ,k=-1)
        elif n_tmp==2:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file_c.dataobj/self.nifti_file_c.dataobj.max())*self.Venc*0.5 ,k=-1)
            
        n_tmp=self.ImageType.index('P')
        if n_tmp==0:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file.dataobj), k=-1 ))/self.phase_scale)*self.Venc
        elif n_tmp==1:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_b.dataobj), k=-1 ))/self.phase_scale)*self.Venc
        elif n_tmp==2:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_c.dataobj), k=-1 ))/self.phase_scale)*self.Venc


    
    
    def set_new_data(self):
        """Load new NIFTI data files."""
        raw_img_filename = filedialog.askopenfilename(
            initialdir="/", title="Select file",
            filetypes=(("NIFTI", "*.nii"), ("all files", "*.*")), multiple=True
        )
        self.raw_img = nib.load(raw_img_filename[0])
        self.raw_img_filename = raw_img_filename[0]
        self.img = nib.as_closest_canonical(self.raw_img)
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        with open(os.path.splitext(raw_img_filename[0])[0]+'.json') as f:
            self.img.json_header = json.load(f)
        ImageType_a = self.img.json_header['ImageType'][2]

        self.raw_img_b = nib.load(raw_img_filename[1])
        self.raw_img_filename_b = raw_img_filename[1]
        self.img_b = nib.as_closest_canonical(self.raw_img_b)
        self.nifti_file_b = nib.Nifti1Image(self.img_b.dataobj, self.img_b.affine, header=self.img_b.header)
        with open(os.path.splitext(raw_img_filename[1])[0]+'.json') as f:
            self.img_b.json_header = json.load(f)
        ImageType_b = self.img_b.json_header['ImageType'][2]
    
        self.raw_img_c = nib.load(raw_img_filename[2])
        self.raw_img_filename_c = raw_img_filename[2]
        self.img_c = nib.as_closest_canonical(self.raw_img_c)
        self.nifti_file_c = nib.Nifti1Image(self.img_c.dataobj, self.img_c.affine, header=self.img_c.header)
        with open(os.path.splitext(raw_img_filename[2])[0]+'.json') as f:
            self.img_c.json_header = json.load(f)
        ImageType_c = self.img_c.json_header['ImageType'][2]
        self.ImageType = [ImageType_a, ImageType_b, ImageType_c]
    
        self.Inv_image_indx = 1

        file_type = os.path.splitext(raw_img_filename[0])[1]
        if file_type == '.PAR':
            self.Venc = self.img.header.general_info.get('phase_enc_velocity')[2]
            self.Venc_from_hdr = True
        else:
            self.Venc = 100
            self.Venc_from_hdr = False
    
        
        n_tmp = self.ImageType.index('P')
        if n_tmp == 0:
            phase_data = np.squeeze(self.nifti_file.dataobj)
        elif n_tmp==1:
            phase_data = np.squeeze(self.nifti_file_b.dataobj)
        elif n_tmp==2:
            phase_data = np.squeeze(self.nifti_file_c.dataobj)
        
        # Determine phase scale factor from the phase image data
        phase_max_abs = max(abs(phase_data.max()), abs(phase_data.min()))
        
        if phase_max_abs > 2048 and phase_max_abs <= 4096:
            self.phase_scale = 4096
        elif phase_max_abs > 4096 and phase_max_abs <= 32768:
            self.phase_scale = 32768
        elif phase_max_abs <= 1.0:
            self.phase_scale = 1.0
        else:
            self.phase_scale = phase_max_abs
        
        print(f"Phase scale factor determined: {self.phase_scale} (data max: {phase_max_abs})")
    
        n_tmp=self.ImageType.index('M')
        if n_tmp==0:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file.dataobj/self.nifti_file.dataobj.max())*self.Venc*0.7, k=-1) 
        elif n_tmp==1:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file_b.dataobj/self.nifti_file_b.dataobj.max())*self.Venc*0.7, k=-1) 
        elif n_tmp==2:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file_c.dataobj/self.nifti_file_c.dataobj.max())*self.Venc*0.7,k=-1)   
        
        n_tmp=self.ImageType.index('MAG')
        if n_tmp==0:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file.dataobj/self.nifti_file.dataobj.max())*self.Venc*0.9 ,k=-1) 
        elif n_tmp==1:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file_b.dataobj/self.nifti_file_b.dataobj.max())*self.Venc*0.9 ,k=-1)
        elif n_tmp==2:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file_c.dataobj/self.nifti_file_c.dataobj.max())*self.Venc*0.9 ,k=-1)
        
        n_tmp=self.ImageType.index('P')
        if n_tmp==0:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file.dataobj), k=-1 ))/self.phase_scale)*self.Venc
        elif n_tmp==1:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_b.dataobj), k=-1 ))/self.phase_scale)*self.Venc
        elif n_tmp==2:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_c.dataobj), k=-1 ))/self.phase_scale)*self.Venc
    
        global nifti_files
        nifti_files=True
    
        global Disp_image_str
        global colormap_str
        global imgFrame
        imgFrame=1
        slider_scale.set(1)
        slider_scale['to'] = 1
        slider_scale['from'] = self.Vel_image.shape[2]
    
        save_str=[self.raw_img_filename[0:-4]+'_Flow_data.csv']
        entry_save_filename.delete(0, 'end')
        entry_save_filename.insert(END, save_str[0])
    
    # Reset ROIs if shape changes
        if ROI_art.BWMask.shape != self.Vel_image.shape:
            ROI_art.polygon=[None]*self.Vel_image.shape[2]
            ROI_art.BWMask=self.Vel_image*False
            ROI_art.flag=[0]*(self.Vel_image.shape[2])
    
    # Calculate mean velocity image
        self.Mean_vel_image_tmp=np.mean(self.Vel_image[:,:,0:3],2)
        self.Mean_vel_image=np.repeat(self.Mean_vel_image_tmp[:,:,np.newaxis], self.Vel_image.shape[2], axis=2)
    
        change_image(Disp_image_str,colormap_str)


if not nifti_files:
    Img_data=Img_data_parrec(raw_img_filename)  
    Img_data.Mean_vel_image_tmp=np.mean(Img_data.Vel_image[:,:,0:3],2)
    Img_data.Mean_vel_image=np.repeat( Img_data.Mean_vel_image_tmp[:,:,np.newaxis],Img_data.Vel_image.shape[2], axis=2)    

if nifti_files:
    Img_data = Img_data_nii([raw_img_filename, raw_img_filename_mod, raw_img_filename_mag])
    Img_data.Mean_vel_image_tmp = np.mean(Img_data.Vel_image[:,:,0:3], 2)
    Img_data.Mean_vel_image = np.repeat(Img_data.Mean_vel_image_tmp[:,:,np.newaxis], Img_data.Vel_image.shape[2], axis=2)
    


# =============================================================================
# IMAGE DISPLAY FUNCTIONS
# =============================================================================

def displayed_image(disp_name):
    """
    Return the appropriate image array based on display type selection.
    
    Parameters
    ----------
    disp_name : str
        Image type identifier:
        - 'vel': Velocity image (phase data converted to cm/s)
        - 'mod': Modulus image (anatomical reference)
        - 'thres': Phase magnitude image (for vessel identification)
        - 'ROI': Binary ROI mask (scaled for visibility)
        - 'Mean_vel': Mean velocity map across selected frames
    
    Returns
    -------
    ndarray
        3D image array (x, y, frames) for the selected display type
    """
    if disp_name == 'vel':
        return Img_data.Vel_image
    elif disp_name == 'mod':
        return Img_data.Mod_image
    elif disp_name == 'thres':
        return Img_data.Thres_image
    elif disp_name == 'ROI':
        return ROI_art.BWMask*100
    elif disp_name == 'Mean_vel':
        return Img_data.Mean_vel_image
        
    
Disp_image=displayed_image(Disp_image_str)


# -----------------------------------------------------------------------------
# Image Type Change Functions
# -----------------------------------------------------------------------------
# These functions handle switching between different image display modes
# via menu selection or keyboard shortcuts

def change_image_type_str_vel(self=''):
    """Switch display to velocity image (Ctrl+1)."""
    global Disp_image_str
    Disp_image_str='vel'
    change_image(Disp_image_str,colormap_str)

def change_image_type_str_mod(self=''):
    """Switch display to modulus image (Ctrl+2)."""
    global Disp_image_str
    Disp_image_str='mod'
    change_image(Disp_image_str,colormap_str)

def change_image_type_str_thres(self=''):
    """Switch display to threshold image (Ctrl+3)."""
    global Disp_image_str
    Disp_image_str='thres'
    change_image(Disp_image_str,colormap_str)
    
def change_image_type_str_ROI(self=''):
    """Switch display to ROI mask overlay (Ctrl+4)."""
    global Disp_image_str
    Disp_image_str='ROI'
    change_image(Disp_image_str,colormap_str)

def change_image_type_str_Mean_vel_map(self=''):
    """Switch display to mean velocity map (Ctrl+5)."""
    global Disp_image_str
    Disp_image_str='Mean_vel'
    change_image(Disp_image_str,colormap_str)


# -----------------------------------------------------------------------------
# Colormap Change Functions
# -----------------------------------------------------------------------------
# These functions handle switching between different colormaps

def change_cmap_jet(self='jet'):
    """Switch to Jet colormap (Ctrl+Q)."""
    global colormap_str
    colormap_str='jet'
    change_image(Disp_image_str,colormap_str)

def change_cmap_gray(self='gray'):
    """Switch to Grayscale colormap (Ctrl+W)."""
    global colormap_str
    colormap_str='gray'
    change_image(Disp_image_str,colormap_str)

def change_cmap_viridis(self='viridis'):
    """Switch to Viridis colormap (Ctrl+E)."""
    global colormap_str
    colormap_str='viridis'
    change_image(Disp_image_str,colormap_str)
    
def popup_change_colorbar():
    """
    Open a popup dialog to adjust colorbar limits.
    
    Creates a small window with entry fields for minimum and maximum
    colorbar values. Changes are applied when the 'Update Colorbar' 
    button is clicked.
    """
    popup_change_colorbar = tk.Toplevel(window)
    popup_change_colorbar.title('Colorbar limits')
    popup_cb_width = max(150, int(window.winfo_screenwidth()*0.10))
    popup_cb_height = max(180, int(window.winfo_screenheight()*0.22))
    popup_change_colorbar.geometry( str(int(popup_cb_width)) +'x'+ str(int(popup_cb_height)) )
    popup_change_colorbar.resizable(True,True)
    
    entry_width = max(7, int(7 * dpi_scale))
    Min_entry_text = Label(popup_change_colorbar, text="Min:",width = entry_width,padx=0)
    Min_entry_text.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='n',pady=10,padx=0)
    Min_entry = Entry(popup_change_colorbar,width=entry_width)
    Min_entry.grid(row=1,rowspan=1,column=1,columnspan=1,sticky='nw',pady=10,padx=0)
    Min_entry.insert(END, '-15')
    
    Max_entry_text = Label(popup_change_colorbar, text="Max:",width = entry_width,padx=0)
    Max_entry_text.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='n',pady=0,padx=0)
    Max_entry = Entry(popup_change_colorbar,width=entry_width)
    Max_entry.grid(row=0,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
    Max_entry.insert(END, '15')
    
    def update_colorbar():
        climits.lims=[float(Min_entry.get()),float(Max_entry.get())]
        change_image(Disp_image_str,colormap_str)
        
    
    UpdateCB_button = Button(master = popup_change_colorbar,
                      height = button_height,
                      width = max(9, int(9 * dpi_scale)),
                     text = "Update Colorbar", command=update_colorbar)
    UpdateCB_button.grid(row=2,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)


# =============================================================================
# MENU BAR SETUP
# =============================================================================
## Create topmenu 
menubar = Menu(window)
window.config(menu=menubar)


# New file menu
analysismenu = Menu(menubar, tearoff=1)
analysis_submenu = Menu(analysismenu, tearoff=1)
menubar.add_cascade(label="File", menu=analysismenu)
analysismenu.add_command(label="Load new PAR/REC file", command=partial(Img_data_parrec.set_new_data, Img_data),accelerator="Control+n")
window.bind_all('<Control-Key-n>', func=partial(Img_data_parrec.set_new_data, Img_data))

analysismenu.add_command(label="Load new NIFTI file", command=partial(Img_data_nii.set_new_data, Img_data),accelerator="Control+m")
window.bind_all('<Control-Key-m>', func=partial(Img_data_nii.set_new_data, Img_data) )


# Image type topmenu
imagemenu = Menu(menubar, tearoff=1)
submenu = Menu(imagemenu, tearoff=1)
menubar.add_cascade(label="Image", menu=imagemenu)
imagemenu.add_cascade(label="Change image type", menu=submenu)

submenu.add_radiobutton(label="Velocity",accelerator="Control+1",command = change_image_type_str_vel)
submenu.add_radiobutton(label="Modulus",accelerator="Control+2",command = change_image_type_str_mod)
submenu.add_radiobutton(label="Magnitude" ,command  = change_image_type_str_thres,accelerator="Control+3")
submenu.add_radiobutton(label="ROI mask" ,command  = change_image_type_str_ROI,accelerator="Control+4")
submenu.add_radiobutton(label="Mean velocity" ,command  = change_image_type_str_Mean_vel_map,accelerator="Control+5")

window.bind_all('<Control-Key-1>', func=change_image_type_str_vel ) # Binds keyboard shortcut to topbar 
window.bind_all('<Control-Key-2>', func=change_image_type_str_mod )
window.bind_all('<Control-Key-3>', func=change_image_type_str_thres )
window.bind_all('<Control-Key-4>', func=change_image_type_str_ROI )
window.bind_all('<Control-Key-5>', func=change_image_type_str_Mean_vel_map )

# Colormap topmenu
subsubmenu = Menu(submenu, tearoff=1)
imagemenu.add_cascade(label="Change colorbar", menu=subsubmenu)
subsubmenu.add_command(label="Jet",accelerator="Control+q", command = change_cmap_jet)
subsubmenu.add_command(label="Gray",accelerator="Control+w", command = change_cmap_gray)
subsubmenu.add_command(label="Viridis",accelerator="Control+e" ,command = change_cmap_viridis)
subsubmenu.add_command(label="Change colorbar limits",accelerator="Control+l" ,command = popup_change_colorbar)

window.bind_all('<Control-Key-q>', func=change_cmap_jet ) # Binds keyboard shortcut to topbar 
window.bind_all('<Control-Key-w>', func=change_cmap_gray )
window.bind_all('<Control-Key-e>', func=change_cmap_viridis )

# Create grid for GUI with proper weights for resizing
Grid.rowconfigure(window, 0, weight=0)  # Header row - fixed
Grid.rowconfigure(window, 1, weight=1)  # Main content row - expandable
Grid.rowconfigure(window, 2, weight=1)  # Content row - expandable
Grid.rowconfigure(window, 3, weight=1)  # Content row - expandable
Grid.rowconfigure(window, 4, weight=1)  # Content row - expandable
Grid.rowconfigure(window, 5, weight=0)  # Toolbar row - fixed
Grid.columnconfigure(window, 0, weight=0)  # Button column - fixed
Grid.columnconfigure(window, 1, weight=1)  # Image column - expandable
Grid.columnconfigure(window, 4, weight=0)  # Slider column - fixed
Grid.columnconfigure(window, 5, weight=1)  # Flow plot column - expandable

    
# Create figure frame for displaying MRI image
fig = plt.figure(figsize=(fig_size_main, fig_size_main), dpi=fig_dpi) # Scaled to screen size
ax = fig.add_subplot(111)

class climits():
    lims=[-15,40]

global pcm_plot
pcm_plot=ax.imshow(Disp_image[:,:,imgFrame-1],cmap=plt.get_cmap(colormap_str),vmin=climits.lims[0],vmax=climits.lims[1], interpolation='none')
ax.set_xticks([]) # Removes x axis ticks
ax.set_yticks([]) # Removes x axis ticks
fig.tight_layout() # Tight layout

canvas = FigureCanvasTkAgg(fig, master=window)  # Draws the figure in the tkinter GUI window
canvas.draw()
canvas.get_tk_widget().grid(row=1,rowspan=4,column=1,columnspan=1,padx=0,pady=0,sticky='nsew') # place in grid with sticky for resizing

# Create line data for plotting ROI data
ROI_line_plot, = ax.plot([], [], '.w-')
fig_flow = plt.figure(figsize=(fig_size_flow, fig_size_flow), dpi=110)  # Scaled to screen size
ax_flow = fig_flow.add_subplot(111)


# Create line plot for flow data
ax_flow.set_position([0.2, 0.15, 0.6, 0.6])
ax_flow.tick_params(labelsize=4)
ax_flow.set_ylabel('Flow [ml/min]', fontsize = 5.0) # Y label
ax_flow.set_xlabel('Index', fontsize = 5.0) # Y label
Flow_line_plot, = ax_flow.plot([], [], '.r-')

canvas_flow = FigureCanvasTkAgg(fig_flow, master=window)  # Draws the figure in the tkinter GUI window
canvas_flow.draw()
canvas_flow.get_tk_widget().grid(row=1,rowspan=4,column=5,columnspan=1,padx=20,pady=50,sticky='nsew') # place in grid with sticky


# Create frame for the navigation toolbar
ToolbarFrame = Frame(window)
ToolbarFrame.grid(row=5, column=1,rowspan=1,columnspan=1,padx=0,pady=0,sticky='sw')
toobar = NavigationToolbar2Tk(canvas, ToolbarFrame)

# =============================================================================
# IMAGE UPDATE FUNCTIONS
# =============================================================================

def change_image(image_str, cmp_str):
    """
    Update the displayed image with new type and/or colormap.
    
    Parameters
    ----------
    image_str : str
        Image type identifier ('vel', 'mod', 'thres', 'ROI', 'Mean_vel')
    cmp_str : str
        Matplotlib colormap name ('jet', 'gray', 'viridis')
    
    Returns
    -------
    str
        Current display image string
    
    Notes
    -----
    Handles image dimension changes when loading new data by recreating
    the imshow object if necessary.
    """
    Disp_image=displayed_image(image_str)
    global pcm_plot

    if pcm_plot.get_array().shape!=Disp_image.shape[0:2]: 
        pcm_plot=ax.imshow(Disp_image[:,:,imgFrame-1],cmap=plt.get_cmap(colormap_str),vmin=climits.lims[0],vmax=climits.lims[1], interpolation='none')
        fig.tight_layout() # Tight layout
        ax.set_xticks([]) # Removes x axis ticks
        ax.set_yticks([]) # Removes x axis ticks
        print('Changed image')
    pcm_plot.set_data(Disp_image[:,:,int(imgFrame-1)])
#   ax.set_ybound(lower=1, upper=Img_data.Vel_image.shape[0])
#    ax.set_xbound(lower=1, upper=Img_data.Vel_image.shape[1])

    
    pcm_plot.set_cmap(cmp_str)
    pcm_plot.set_clim(climits.lims[0],climits.lims[1])
    canvas.draw()
    return Disp_image_str


def update_image(self):
    """
    Update display when cardiac frame changes via slider.
    
    Parameters
    ----------
    self : str or int
        Frame number from slider (1-based indexing)
    
    Returns
    -------
    int
        Current frame index
    
    Notes
    -----
    Also updates the ROI polygon overlay if one exists for the current frame.
    """
    global imgFrame
    imgFrame=int(self)
    Disp_image=displayed_image(Disp_image_str)
    pcm_plot.set_data(Disp_image[:,:,int(imgFrame-1)])
    
    ROI_as_array = np.array(ROI_art.polygon[int(imgFrame-1)])
    if np.equal(ROI_as_array, None).all():
        ROI_line_plot.set_ydata([])
        ROI_line_plot.set_xdata([])
    else:
        ROI_as_array_tmp=np.append(ROI_as_array[:,:], ROI_as_array[0,:]).reshape(ROI_as_array.shape[0]+1,ROI_as_array.shape[1])
        ROI_line_plot.set_ydata(ROI_as_array_tmp[:,1])
        ROI_line_plot.set_xdata(ROI_as_array_tmp[:,0])
        
    canvas.draw()    
    return imgFrame


# Create slider for changing frame.
slider_scale = Scale(window, to=1, from_=Img_data.Vel_image.shape[2], width=20, length=slider_length, command=update_image)
slider_scale.grid(row=1,rowspan=4,column=4,columnspan=1,padx=0,pady=0,sticky='ns')


def key_arrow_up(self=''):
    global imgFrame
    if imgFrame>(Img_data.Vel_image.shape[2]-1):
        imgFrame=imgFrame
    else:
        imgFrame=imgFrame+1
    update_image(imgFrame)
    slider_scale.set(imgFrame)

window.bind_all('<Up>', func=key_arrow_up )       
    
def key_arrow_down(self=''):
    global imgFrame
    if imgFrame<2:
        imgFrame=imgFrame
    else:
        imgFrame=imgFrame-1
    update_image(imgFrame)
    slider_scale.set(imgFrame)

window.bind_all('<Down>', func=key_arrow_down )    


# Add headline text
headline_text = Label(window, text="Draw ROI to calculate flow")
headline_text.grid(row=0,rowspan=1,column=0,columnspan=2,sticky='nw')

# =============================================================================
# ROI (REGION OF INTEREST) CLASSES
# =============================================================================

class ROI_art(object):
    """
    Class to store and manage ROI data for vessel segmentation.
    
    This class holds the polygon vertices and binary masks for ROIs
    across all cardiac frames. ROIs can be drawn manually or automatically
    using the region growing algorithm.
    
    Class Attributes
    ----------------
    polygon : list
        List of polygon vertex arrays, one per cardiac frame.
        Each element is either None or an Nx2 array of (x,y) coordinates.
    BWMask : ndarray
        3D binary mask array (x, y, frames) where True indicates pixels
        inside the ROI.
    flag : list
        List of integers indicating ROI status per frame:
        0 = no ROI defined, 1 = ROI defined
    
    Class Methods
    -------------
    set_polygon(verts)
        Set ROI polygon for current frame from vertex coordinates
    loaded_roi_from_file(Loaded_ROI_data)
        Load ROI data from a saved NPZ file
    """
    polygon=[None] * Img_data.Vel_image.shape[2]
    BWMask=Img_data.Vel_image*False
    flag=[0]*(Img_data.Vel_image.shape[2])

    def set_polygon(verts):
        """
        Set ROI polygon for the current frame.
        
        Parameters
        ----------
        verts : array-like
            Nx2 array of polygon vertex coordinates (x, y)
        
        Notes
        -----
        Creates a binary mask by testing which pixels fall inside the polygon
        using matplotlib.path.Path.contains_points().
        """
        print('Updating ROI')
        ROI_art.polygon[int(imgFrame-1)]=verts # ROI as polygon as input
        ROI_art.flag[int(imgFrame-1)]=1
        path = Path(verts)
        x, y = np.meshgrid(np.arange(Img_data.Vel_image.shape[1]), np.arange(Img_data.Vel_image.shape[0]))
        x, y = x.flatten(), y.flatten()
        points = np.vstack((x,y)).T
        grid = path.contains_points(points) # Find voxels inside polygon
        grid = grid.reshape((Img_data.Vel_image.shape[0],Img_data.Vel_image.shape[1]))        
        ROI_art.BWMask[:,:,int(imgFrame-1)]=grid
    def loaded_roi_from_file(Loaded_ROI_data):
        """
        Load ROI data from a saved NPZ file.
        
        Parameters
        ----------
        Loaded_ROI_data : dict-like
            NPZ file contents with keys:
            - 'PCMROI_poly': polygon vertices per frame
            - 'PCMROI_BWMask': binary mask array
            - 'PCMROI_flag': ROI status flags
        """
        ROI_art.polygon=Loaded_ROI_data['PCMROI_poly'].tolist()
        ROI_art.BWMask=Loaded_ROI_data['PCMROI_BWMask']
        ROI_art.flag=Loaded_ROI_data['PCMROI_flag'].tolist()
        global imgFrame
        update_image(imgFrame)
        

class ROIPolygon(object):
    """
    Interactive polygon drawing tool for manual ROI delineation.
    
    Wraps matplotlib's PolygonSelector to enable interactive drawing
    of vessel ROIs directly on the image display.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to attach the polygon selector to
    row, col : int
        Image dimensions (for reference)
    
    Attributes
    ----------
    canvas : matplotlib.backends.backend_tkagg.FigureCanvasTkAgg
        The figure canvas for drawing updates
    poly : matplotlib.widgets.PolygonSelector
        The polygon selector widget
    
    Methods
    -------
    onselect(verts)
        Callback when polygon selection is completed
    """
    def __init__(self, ax, row, col):
        self.canvas = ax.figure.canvas
        global PS
        PS = PolygonSelector(ax,
                                    self.onselect,
                                    props = dict(color = 'm', alpha = 1),
                                    handle_props = dict(mec = 'm', mfc = 'm', alpha = 1,markersize=2),grab_range=10)
        ROIPolygon.poly=PS 
    def onselect(self, verts):
        ROI_art.set_polygon(verts)
        self.canvas.draw_idle()
        
        ROI_as_array = np.array(PS.verts)
        ROI_as_array_tmp=np.append(ROI_as_array[:,:], ROI_as_array[0,:]).reshape(ROI_as_array.shape[0]+1,ROI_as_array.shape[1])
        ROI_line_plot.set_ydata(ROI_as_array_tmp[:,1])
        ROI_line_plot.set_xdata(ROI_as_array_tmp[:,0])
        global Disp_image_str
        global colormap_str
        change_image(Disp_image_str,colormap_str)

# Global polygon selector reference
global PS
PS=None

# -----------------------------------------------------------------------------
# ROI Control Functions
# -----------------------------------------------------------------------------
    
def add_roi_func():
    """
    Initialize a new ROI polygon for the current frame.
    
    Creates a new PolygonSelector widget and activates it for drawing.
    Any existing polygon selector is deactivated first.
    """
    global PS
    if PS is not None:
        PS.set_visible(False)
        PS.set_active(False)
    roi_poly=ROIPolygon(pcm_plot.axes, Img_data.Vel_image.shape[0], Img_data.Vel_image.shape[2])
    UpdateROI_button.configure(text='Stop ROI edit')
    
def copy_roi_func():
    """
    Copy the current frame's ROI to all other frames.
    
    If the current frame has an ROI, copies it to all frames.
    Otherwise, finds the first frame with an ROI and copies that.
    Useful for vessels that don't move significantly across the cardiac cycle.
    """
    first_indx=ROI_art.flag.index(1)
    if ROI_art.flag[int(imgFrame-1)]==1:
        for x in range(len(ROI_art.polygon)):
            ROI_art.polygon[x]=ROI_art.polygon[int(imgFrame-1)]
            ROI_art.BWMask[:,:,x]=ROI_art.BWMask[:,:,int(imgFrame-1)]
            ROI_art.flag[x]=1
    else:
        for x in range(len(ROI_art.polygon)): 
            ROI_art.polygon[x]=ROI_art.polygon[first_indx]
            ROI_art.BWMask[:,:,x]=ROI_art.BWMask[:,:,first_indx]
            ROI_art.flag[x]=1

def update_roi_func():
    """
    Toggle ROI editing mode on/off.
    
    When activated, allows interactive modification of the current ROI polygon.
    Button text changes to indicate current state.
    """
    global PS
    if PS.get_active():
        PS.set_visible(False)
        PS.set_active(False)
        canvas.draw() 
        UpdateROI_button.configure(text='Edit ROI')
    else:
        print('Updating ROI')

        PS.set_visible(True)
        PS.set_active(True)
        canvas.draw()
        UpdateROI_button.configure(text='Stop ROI edit')


# =============================================================================
# REGION GROWING ALGORITHM
# =============================================================================

class RegGrow():
    """
    Automatic ROI delineation using region growing algorithm.
    
    This class implements a seeded region growing algorithm for automatic
    vessel segmentation. The user clicks on a seed point, and the algorithm
    grows outward to include neighboring pixels that meet velocity threshold
    and distance criteria.
    
    Class Attributes
    ----------------
    nRow, nCol, nSlice : int
        Image dimensions
    qu : int or ndarray
        Queue of pixels to process
    ginput_input : tuple
        Mouse click coordinates
    btm_press_event : int
        Event connection ID for mouse callback
    
    Class Methods
    -------------
    create_mask(event)
        Process mouse click and run region growing
    Calc_Auto_ROI()
        Initialize automatic ROI mode (changes cursor to crosshair)
    
    Algorithm Parameters (from GUI entries)
    ---------------------------------------
    Max. dist (e1) : int
        Maximum distance from seed point (pixels)
    Threshold (e2) : float
        Minimum velocity threshold (cm/s) for inclusion
    
    Notes
    -----
    The algorithm:
    1. User clicks seed point on vessel
    2. Region grows to include neighbors with velocity >= threshold
    3. Growth limited by maximum distance from seed
    4. Process repeats for each frame, tracking the vessel
    5. Contours are extracted and converted to polygon ROIs
    """

    nRow, nCol, nSlice=Img_data.Vel_image.shape
    qu=1
    ginput_input=[]
    btm_press_event=[]
    
    def create_mask(event):
        RegGrow.nRow=Img_data.Vel_image.shape[0]
        RegGrow.nCol=Img_data.Vel_image.shape[1]
        RegGrow.nSlice==Img_data.Vel_image.shape[2]
        window.config(cursor='')
        ginput_input= event.xdata, event.ydata
        maxDist=int(e1.get())
        ThresVal=int(e2.get()) 
        p=0
        fig.canvas.callbacks.disconnect(RegGrow.btm_press_event)
        seed_tmp=np.round(ginput_input)
        seed=np.flip(seed_tmp)
        Reg_mask=np.zeros(Img_data.Vel_image.shape)
        for nn in range(Img_data.Vel_image.shape[2]):
            queue = seed
            Imax=int(seed[0])
            Jmax=int(seed[1])
            print(nn)
            while queue.any():
                if queue.ndim==1:
                    xv = int(queue[0])
                    yv = int(queue[1])
                else:
                    xv = int(queue[0][0])
                    yv = int(queue[0][1])
                    
                for n in [-1,0,1]:
                    # print(n)
                    for m in [-1,0,1]:
                        #print(m)
                        if xv+n > 0  and  xv+n <= RegGrow.nRow and yv+m > 0 and  yv+m <= RegGrow.nCol and any([n, m]) and Reg_mask[xv+n,yv+m,nn]==0 and np.sqrt( (xv+n-Imax)**2 + (yv+m-Jmax)**2 ) < maxDist and Img_data.Vel_image[xv+n,yv+m,nn] >= ThresVal:
                            Reg_mask[(xv+n, yv+m,nn)]=1
                            queue=np.vstack((queue, np.array([xv+n, yv+m])))
                        #print(queue)   
      
                queue = numpy.delete(queue, (0), axis=0)
                RegGrow.qu=queue  
                masked_inp_tmp=Img_data.Vel_image[:,:,nn]*Reg_mask[:,:,nn]
                New_seed = np.where(masked_inp_tmp == np.amax(masked_inp_tmp))
                seed=np.array( [New_seed[0][0], New_seed[1][0]])
                
        for mm in range(Img_data.Vel_image.shape[2]):
            found_contours=measure.find_contours(Reg_mask[:,:,mm], 0.5)
            cnt_length=np.zeros(len(found_contours))
            for bb in range(len(found_contours)):
                    cnt_length[bb]=found_contours[bb].shape[0]
            max_count_idx=np.where(cnt_length == np.amax(cnt_length))        
                    
            ROI_art.polygon[mm]=np.fliplr(found_contours[max_count_idx[0][0]])
            path = Path(np.fliplr(found_contours[max_count_idx[0][0]]))
            x, y = np.meshgrid(np.arange(Img_data.Vel_image.shape[1]), np.arange(Img_data.Vel_image.shape[0]))
            x, y = x.flatten(), y.flatten()
            points = np.vstack((x,y)).T
            grid = path.contains_points(points)
            grid = grid.reshape((Img_data.Vel_image.shape[0],Img_data.Vel_image.shape[1]))        
            ROI_art.BWMask[:,:,mm]=grid
            ROI_art.flag[mm]=1
            

        global PS 
        ROIPolygon(pcm_plot.axes, Img_data.Vel_image.shape[0], Img_data.Vel_image.shape[2])
        PS.set_visible(False)
        PS.set_active(False) 
        ROI_as_array = np.array(ROI_art.polygon[int(imgFrame-1)])
        ROI_line_plot.set_ydata(ROI_as_array[:,1])
        ROI_line_plot.set_xdata(ROI_as_array[:,0])
        fig.canvas.draw()
        

    def Calc_Auto_ROI() :
        window.config(cursor='plus green white')
        RegGrow.btm_press_event=fig.canvas.callbacks.connect('button_press_event', RegGrow.create_mask)
            

    global PS 
    ROIPolygon(pcm_plot.axes, Img_data.Vel_image.shape[0], Img_data.Vel_image.shape[2])
    PS.set_visible(False)
    PS.set_active(False) 
    ROI_as_array = np.array(ROI_art.polygon[int(imgFrame-1)])
    #ROI_line_plot.set_ydata(ROI_as_array[:,1])
    #ROI_line_plot.set_xdata(ROI_as_array[:,0])
    fig.canvas.draw()
        

    def Calc_Auto_ROI() :
        window.config(cursor='plus green white')
        RegGrow.btm_press_event=fig.canvas.callbacks.connect('button_press_event', RegGrow.create_mask)
    
    

# Create buttons for ROI controls 
ROI_button_group = LabelFrame(window, text='ROI analysis', borderwidth=2,relief='solid')
ROI_button_group.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='nw',padx=10,pady=0)

# Add new ROI
AddROI_button = Button(master = ROI_button_group,
                     height = button_height,
                     width = button_width,
                    text = "Add ROI", command=add_roi_func)
AddROI_button.grid(row=0,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)

# Update ROI
UpdateROI_button = Button(master = ROI_button_group,
                     height = button_height,
                     width = button_width,
                    text = "Edit ROI", command=update_roi_func)
UpdateROI_button.grid(row=1,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)

# Copy ROI to remaning frames 
CopyROI_button = Button(master = ROI_button_group,
                     height = button_height,
                     width = button_width,
                    text = "Copy ROI to \n all frames", command=copy_roi_func)
CopyROI_button.grid(row=3,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)

# Call region growing algorithm
AutoROI_button = Button(master = ROI_button_group,
                     height = button_height + 1,
                     width = button_width,
                    text = "Automatic \n delineation", command=RegGrow.Calc_Auto_ROI)
AutoROI_button.grid(row=4,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)
Max_dist_text = Label(ROI_button_group, text="Max. dist:",padx=0)
Max_dist_text.grid(row=5,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=0)

# Entries for region growing algorithm 
e1 = Entry(ROI_button_group,width=4)
e1.grid(row=5,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
e1.insert(END, '8')

Threshold_text = Label(ROI_button_group, text="Threshold:")
Threshold_text.grid(row=6,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=10)
e2 = Entry(ROI_button_group,width=4)
e2.grid(row=6,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
e2.insert(END, '10')

MF_text = Label(ROI_button_group, text="Avg. Vel. Frames:")
MF_text.grid(row=7,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=10)
e3 = Entry(ROI_button_group,width=4)
e3.grid(row=7,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
e3.insert(END, '0,1,2')


def clear_all_rois():
    """
    Clear all ROIs and reset the analysis.
    
    Prompts for confirmation, then:
    - Resets all ROI polygons and masks to empty
    - Clears the flow plot
    - Resets flow output data structures
    - Deactivates any active polygon selector
    - Refreshes the display
    """
    result = messagebox.askyesno("Clear ROI", 
                                "Are you sure you want to clear ROI?")
    
    if result:
        # Reset ROI data structures
        ROI_art.polygon = [None] * Img_data.Vel_image.shape[2]
        ROI_art.BWMask = Img_data.Vel_image * False
        ROI_art.flag = [0] * (Img_data.Vel_image.shape[2])
        
        # Clear ROI line plot
        ROI_line_plot.set_ydata([])
        ROI_line_plot.set_xdata([])
        
        # Clear flow plot
        ax_flow.clear()
        ax_flow.set_ylabel('Flow [ml/min]', fontsize=5.0)
        ax_flow.set_xlabel('Index', fontsize=5.0)
        ax_flow.set_xlim([1, Img_data.Vel_image.shape[2]+1])
        
        # Reset flow output data
        Flow_output.Flows = np.empty((1, Img_data.Vel_image.shape[2],))
        Flow_output.Flows[:] = np.nan
        Flow_output.Velocity = np.empty((1, Img_data.Vel_image.shape[2],))
        Flow_output.Velocity[:] = np.nan
        Flow_output.CS_area = np.empty((1, Img_data.Vel_image.shape[2],))
        Flow_output.CS_area[:] = np.nan
        Flow_output.flow_str = ''
        
        # Clear flow text display
        Flow_text_str.configure(text='')
        
        
        # Deactivate polygon selector if active
        global PS
        if PS is not None:
            PS.set_visible(False)
            PS.set_active(False)
            UpdateROI_button.configure(text='Edit ROI')
        
        # Refresh display
        canvas.draw()
        canvas_flow.draw()
        
        # Update the current frame display
        update_image(imgFrame)
        
        print("ROI cleared successfully")
        
        # Update status message
        Data_saved_txt.configure(text="ROI cleared - ready for new analysis")

       


# Clear all ROIs button
ClearAllROI_button = Button(master=ROI_button_group,
                           height=button_height,
                           width=button_width,
                           text="Clear ROI",
                           command=clear_all_rois,
                           highlightbackground='#ffcccc')  # Light red background to indicate destructive action
ClearAllROI_button.grid(row=8, rowspan=1, column=0, columnspan=2, sticky='nw', pady=2, padx=10)



def inverse_vel_image():
    """
    Invert the velocity image polarity.
    
    Multiplies velocity values by -1 to flip flow direction.
    Useful when vessel flow appears negative due to acquisition orientation.
    Also inverts the mean velocity map for consistency.
    """
    Img_data.Inv_image_indx=Img_data.Inv_image_indx*-1
    if nifti_files:
        Img_data.Venc=float(Venc_ent.get())
    
        n_tmp=Img_data.ImageType.index('P')
        if n_tmp==0:
            Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file.dataobj), k=-1 ))/Img_data.phase_scale)*Img_data.Venc*Img_data.Inv_image_indx
        elif n_tmp==1:
            Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_b.dataobj), k=-1 ))/Img_data.phase_scale)*Img_data.Venc*Img_data.Inv_image_indx
        elif n_tmp==2:
            Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_c.dataobj), k=-1 ))/Img_data.phase_scale)*Img_data.Venc*Img_data.Inv_image_indx
    else: 
        if Img_data.img_data_indx[0][0]==0:
               Img_data.Vel_image=np.rot90(Img_data.img.dataobj_sq[Img_data.phase_image_idx,:,:])*Img_data.Inv_image_indx
        
        if Img_data.img_data_indx[0][0]==1:
            Img_data.Vel_image=np.rot90(Img_data.img.dataobj_sq[:,Img_data.phase_image_idx,:])*Img_data.Inv_image_indx

        if Img_data.img_data_indx[0][0]==2:
              Img_data.Vel_image=np.rot90(Img_data.img.dataobj_sq[:,:,Img_data.phase_image_idx])*Img_data.Inv_image_indx
    
    
    Img_data.Mean_vel_image=Img_data.Mean_vel_image*-1
    
    update_image(imgFrame)


Inv_image_button_group = LabelFrame(window, borderwidth=2,relief='solid')
Inv_image_button_group.grid(row=4,rowspan=1,column=5,columnspan=2,sticky='sw',padx=0,pady=0)
Inv_image_button = Button(master = Inv_image_button_group,
                      height = button_height,
                      width = button_width,
                      text = "Inv. Vel. image", command=inverse_vel_image)
Inv_image_button.grid(row=0,rowspan=1,column=0, columnspan=1,sticky='nw',pady=0,padx=0)  



# =============================================================================
# FLOW CALCULATION
# =============================================================================

class Flow_output:
    """
    Class to store and calculate blood flow measurements.
    
    Computes flow from velocity data within the ROI mask, accounting for
    pixel size to convert velocity (cm/s) to volumetric flow (mL/min).
    
    Class Attributes
    ----------------
    Vel_image_tmp : ndarray
        Velocity image masked by ROI
    Flows : ndarray
        Flow values (mL/min) per cardiac frame
    Velocity : ndarray
        Mean velocity (cm/s) per cardiac frame
    CS_area : ndarray
        Cross-sectional area (mm²) per cardiac frame
    flow_str : str
        Formatted mean flow string for display
    
    Class Methods
    -------------
    Calc_flow()
        Calculate flow from current ROI and velocity data
    
    Flow Calculation
    ----------------
    Flow [mL/min] = mean_velocity [cm/s] × area [mm²] × 60 / 100
    
    Where:
        - mean_velocity: Average velocity within ROI
        - area: Number of ROI pixels × pixel area (from header)
        - Factor 60: converts seconds to minutes
        - Factor 100: converts cm³ to mL
    """
    Vel_image_tmp=ROI_art.BWMask*Img_data.Vel_image
    Flows= np.empty((1,Vel_image_tmp.shape[2],))
    Flows[:] = np.nan
    Velocity= np.empty((1,Vel_image_tmp.shape[2],))
    Velocity[:] = np.nan
    flow_str=''
    CS_area= np.empty((1,Vel_image_tmp.shape[2],))
    CS_area[:] = np.nan
    
    
    def Calc_flow():
        """
        Calculate blood flow from velocity data within the ROI.
        
        For each cardiac frame:
        1. Extracts velocity values within the ROI mask
        2. Computes mean velocity
        3. Calculates cross-sectional area from pixel count
        4. Converts to volumetric flow (mL/min)
        
        Results are stored in class attributes and displayed in the flow plot.
        
        Notes
        -----
        - For NIfTI data, VENC is read from the entry field
        - Pixel size is obtained from image header metadata
        """
    
        if nifti_files:
            Img_data.Venc=float(Venc_ent.get())
            n_tmp=Img_data.ImageType.index('P')
            
            if n_tmp==0:
                Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file.dataobj), k=-1 ))/Img_data.phase_scale)*Img_data.Venc*Img_data.Inv_image_indx
            elif n_tmp==1:
                Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_b.dataobj), k=-1 ))/Img_data.phase_scale)*Img_data.Venc*Img_data.Inv_image_indx
            elif n_tmp==2:
                Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_c.dataobj), k=-1 ))/Img_data.phase_scale)*Img_data.Venc*Img_data.Inv_image_indx
            
            
            MF_indx=[int(x) for x in e3.get().split(',')]
            Img_data.Mean_vel_image_tmp=np.mean(Img_data.Vel_image[:,:,MF_indx],2)
            Img_data.Mean_vel_image=np.repeat( Img_data.Mean_vel_image_tmp[:,:,np.newaxis],Img_data.Vel_image.shape[2], axis=2)
            #Img_data.Vel_image =Img_data.Vel_image

            change_image(Disp_image_str,colormap_str)

        if np.mean(ROI_art.BWMask)==0:
            print('Please draw ROI for calculation flow')

        else:
            Flow_output.Vel_image_tmp=ROI_art.BWMask[:,:,range(0,Img_data.Vel_image.shape[2])]*Img_data.Vel_image
            Flow_output.Flows=np.empty((1,Flow_output.Vel_image_tmp.shape[2],))
            Flow_output.Velocity=np.empty((1,Flow_output.Vel_image_tmp.shape[2],))
            Flow_output.CS_area=np.empty((1,Flow_output.Vel_image_tmp.shape[2],))

            if nifti_files:
                sort_indx=np.argsort(Img_data.img.header.get_data_shape())
                vox_size=Img_data.img.header.get_zooms()[sort_indx[-1]]*Img_data.img.header.get_zooms()[sort_indx[-2]]
            else: 
                vox_size=np.prod( Img_data.img.header.image_defs['pixel spacing'][0:2,0] )
               
                
            for x in range(0,Img_data.Vel_image.shape[2]):
                #ROI_indx, ROI_indy= np.nonzero(Flow_output.Vel_image_tmp[:,:,x])
                ROI_indx, ROI_indy= np.nonzero(ROI_art.BWMask[:,:,x])

                ROI_indx_BW, ROI_indy_BW= np.nonzero(ROI_art.BWMask[:,:,x])
                Velocity_values=Flow_output.Vel_image_tmp[ROI_indx,ROI_indy,x]
                Flow_output.Flows[0][x]=(Velocity_values.mean())*( (ROI_indx.shape[0]*vox_size) /1e2)*60  
                Flow_output.Velocity[0][x]=Velocity_values.mean()
                Flow_output.CS_area[0][x]=ROI_indx.shape[0]*vox_size

               
            print('Flow:')
            print(Flow_output.Flows)
            ax_flow.clear()

            ax_flow.plot(range(1,Img_data.Vel_image.shape[2]+1), Flow_output.Flows[0],'-ro',linewidth=1,markersize=2)
            ax_flow.set_xlim([1,Img_data.Vel_image.shape[2]+1])
            ax_flow.set_ylim([np.min(Flow_output.Flows)*0.9,np.max(Flow_output.Flows)*1.1])

            ax_flow.set_ylabel('Flow [ml/min]', fontsize = 5.0) # Y label
            ax_flow.set_xlabel('Index', fontsize = 5.0) # Y label

            
            canvas_flow.draw()
            Flow_output.flow_str="%5.4f" % Flow_output.Flows.mean()
            Flow_text_str.configure(text=Flow_output.flow_str)


# =============================================================================
# DATA EXPORT
# =============================================================================

class save_ouput_data:
    """
    Class to handle saving analysis results to files.
    
    Exports flow data to CSV and optionally saves:
    - ROI masks as NPZ (numpy archive)
    - ROI masks as NIfTI (for NIfTI input files)
    - Animated GIF of ROI across cardiac cycle
    
    Class Methods
    -------------
    save_data()
        Save all selected output files
    
    Output CSV Columns
    ------------------
    - Flow: Volumetric flow (mL/min) per frame
    - Velocity: Mean velocity (cm/s) per frame  
    - CSarea: Cross-sectional area (mm²) per frame
    """
    
    def save_data(self=''):
        """
        Save flow data and optional ROI files.
        
        Reads filename from entry field and saves:
        - CSV with flow, velocity, and area data
        - NPZ with ROI polygons and masks (if checkbox selected)
        - NIfTI ROI mask (if NIfTI input and checkbox selected)
        - GIF animation (if checkbox selected)
        """
        output_file=entry_save_filename.get()

    
    # Create the dataframe with proper handling
        d = {
            'Flow': Flow_output.Flows[0].tolist(), 
            'Velocity': Flow_output.Velocity[0].tolist(), 
            'CSarea': Flow_output.CS_area[0].tolist()
        }
    
        df = pd.DataFrame(data=d)
        df.to_csv(output_file)
        Data_saved_txt.configure(text='Data saved:'+output_file)

    
        print('Data saved:'+output_file)
    
        if npz_roi.status:
            ROI_filename = os.path.splitext(output_file)[0]+'_ROIs'
            np.savez(ROI_filename, PCMROI_poly=np.array(ROI_art.polygon,dtype='object'), 
                 PCMROI_BWMask=ROI_art.BWMask, PCMROI_flag=ROI_art.flag)
    
        if gif_roi.status:
            if os.path.isdir(os.path.splitext(output_file)[0]+'_ROIgif')==0: 
                os.mkdir(os.path.splitext(output_file)[0]+'_ROIgif')
            fig_flow.savefig(os.path.splitext(output_file)[0]+'_ROIgif/'+os.path.splitext((output_file.replace('/',' ').split(' ')[-1]))[0]+'_Flow.png') 
            for i in range(1,Img_data.Vel_image.shape[2]+1): 
                update_image(i) 
                fig.savefig(os.path.splitext(output_file)[0]+'_ROIgif/'+os.path.splitext((output_file.replace('/',' ').split(' ')[-1]))[0]+'_frame'+str(i)+'.png') 
            create_gif(os.path.splitext(output_file)[0]+'_ROIgif/'+os.path.splitext((output_file.replace('/',' ').split(' ')[-1]))[0]) 
            update_image(imgFrame)
            slider_scale.set(imgFrame) 

        if nifti_files:
            if nii_roi.status:
                raw_img=nib.load(Img_data.raw_img_filename)
                ROI_nii = nib.Nifti1Image(numpy.expand_dims(np.flipud(np.rot90(ROI_art.BWMask)),2), affine=raw_img.affine) 
                nib.save(ROI_nii, os.path.splitext(output_file)[0]+'_ROIs')

def create_gif(path):
    """
    Create an animated GIF from saved frame images.
    
    Parameters
    ----------
    path : str
        Base path for frame images (without _frameN.png suffix)
    
    Notes
    -----
    Expects PNG files named {path}_frame1.png, {path}_frame2.png, etc.
    """
    frames = []
    imgs = glob.glob(path+"_frame*.png")
    for i in imgs:
        new_frame = Image.open(i)
        frames.append(new_frame)

    frames[0].save(path+'.gif', format='GIF',
                   append_images=frames[1:],
                   save_all=True,
                   duration=100, loop=0)
    return


# Button groups for saving data
Save_button_group = LabelFrame(window, text='Save data', borderwidth=2,relief='solid')
Save_button_group.grid(row=7,rowspan=1,column=0,columnspan=5,sticky='nw',padx=10,pady=2)
Save_button = Button(master = Save_button_group,
                     height = button_height,
                     width = max(8, int(8 * dpi_scale)),
                    text = "Save data", command=save_ouput_data.save_data)
Save_button.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='sw',pady=2,padx=10)



class roi_save:
    """
    Helper class to track checkbox state for save options.
    
    Used for the NPZ, GIF, and NIfTI save checkboxes.
    
    Attributes
    ----------
    status : bool
        Current checkbox state (True = checked)
    
    Methods
    -------
    change_status()
        Toggle the status between True and False
    """
    def __init__(self):
        self.status=False
    def change_status(self):
        if self.status==True:
            self.status=False
            # print('Status changes to false')
        else:
            self.status=True
            # print('Status changes to True')
    
npz_roi=roi_save()
npz_roi.change_status()
gif_roi=roi_save()

if nifti_files:
    nii_roi=roi_save()
    nii_roi.change_status()

Save_button.Save_check_roi_npz = Checkbutton(Save_button_group, text='.npz',variable=1, onvalue=1,offvalue=0, command=npz_roi.change_status)
Save_button.Save_check_roi_npz.grid(row=0,rowspan=1,column=2,columnspan=1,sticky='sw',pady=2,padx=0)
Save_button.Save_check_roi_npz.select()
      
Save_check_roi_gif = Checkbutton(Save_button_group, text='.gif',variable=3, onvalue=1,offvalue=0, command=gif_roi.change_status)
Save_check_roi_gif.grid(row=0,rowspan=1,column=4,columnspan=1,sticky='sw',pady=2,padx=0)

if nifti_files:
    Save_button.Save_check_roi_nii = Checkbutton(Save_button_group, text='.nii',variable=2, onvalue=1,offvalue=0, command=nii_roi.change_status)
    Save_button.Save_check_roi_nii.grid(row=0,rowspan=1,column=3,columnspan=1,sticky='sw',pady=2,padx=0)
    Save_button.Save_check_roi_nii.select()


entry_save_filename=Entry(Save_button_group,width=max(80, int(80 * dpi_scale)))
entry_save_filename.grid(row=3,rowspan=1,column=0,columnspan=12,sticky='nw',pady=0,padx=0)
save_str=[Img_data.raw_img_filename[0:-4]+'_Flow_data.csv']
entry_save_filename.insert(END, save_str[0])

Data_saved_txt = Label(Save_button_group)
Data_saved_txt.grid(row=2,rowspan=1,column=0,columnspan=12,sticky='se',pady=0,padx=0)

# Button group for calculating flow
calc_button_group = LabelFrame(window, text='Calculate flow', borderwidth=2,relief='solid')
calc_button_group.grid(row=4,rowspan=1,column=0,columnspan=1,sticky='sw',padx=10,pady=10)

Venc_str = Label(calc_button_group, text="Venc:")
Venc_str.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=0)
if nifti_files:
    Venc_ent = Entry(calc_button_group,width=4)
    Venc_ent.grid(row=0,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
    Venc_ent.insert(END, '100')
else:
    Venc_ent = Label(calc_button_group, text=Img_data.Venc)
    Venc_ent.grid(row=0,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)

CalcFlow_button = Button(master = calc_button_group,
                      height = button_height,
                      width = button_width,
                      text = "Calculate flow", command=Flow_output.Calc_flow)

CalcFlow_button.grid(row=1,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)
Flow_text = Label(calc_button_group, text="Mean flow:")
Flow_text.grid(row=2,rowspan=1,column=0,columnspan=1,sticky='nw')

Flow_text_str = Label(calc_button_group)
Flow_text_str.grid(row=2,rowspan=1,column=1,columnspan=1,sticky='ne')


def load_ROI_file(self=''):
    """
    Load previously saved ROI data from NPZ file.
    
    Opens file dialog to select an NPZ file containing:
    - PCMROI_poly: Polygon vertices per frame
    - PCMROI_BWMask: Binary mask array
    - PCMROI_flag: ROI status flags
    
    Keyboard shortcut: Ctrl+R
    """
    ROI_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("ROI files","*.npz"),("all files","*.*")))
    Loaded_ROI = np.load(ROI_filename,allow_pickle=True)
    ROI_art.loaded_roi_from_file(Loaded_ROI)


# Load ROI shortcut
analysismenu.add_command(label="Load ROI file", command=load_ROI_file, accelerator="Control+r")
window.bind_all('<Control-Key-r>', func=load_ROI_file )

# Save data shortcut
analysismenu.add_command(label="Save flow", command=save_ouput_data.save_data, accelerator="Control+s")
window.bind_all('<Control-Key-s>', func=save_ouput_data.save_data)


# =============================================================================
# PULSATILITY ANALYSIS
# =============================================================================

def CalcPulsatility():
    """
    Open pulsatility analysis window and calculate pulsatility measures.
    
    Calculates two complementary measures of arterial pulsatility:
    
    PI (Pulsatility Index)
    -------------------------------
    PI = (Qmax - Qmin) / Qmean
    
    - Dimensionless ratio describing flow waveform shape
    - Higher values indicate more pulsatile flow
    - Reflects downstream vascular resistance and arterial stiffness
    - Normal cerebral values: 0.6-1.2
    - Does NOT require heart rate input
    
    ΔV (Volume of Arterial Pulsatility)
    ------------------------------------
    ΔV = max(∫[Q-Q̄]dt) - min(∫[Q-Q̄]dt)
    
    - Volume (mL) of cyclic arterial distension per heartbeat
    - Represents the buffering capacity of downstream arterial tree
    - Quantifies how much arteries expand/contract to smooth pulsatile flow
    - Requires heart rate input for correct temporal scaling
    
    
    GUI Elements
    ------------
    - Heart rate entry: Required for ΔV calculation (default: 60 bpm)
    - Flow waveform plot: Shows Q(t) with Qmax, Qmin, Qmean lines
    - Cumulative volume plot: Shows ∫(Q-Q̄)dt with ΔV annotation
    - Results panel: Displays PI, ΔV, and component values
    

    """
    PulsatilityWindow = tk.Toplevel(window)
    puls_width = window.winfo_screenwidth()*0.65
    puls_height = window.winfo_screenheight()*0.65
    PulsatilityWindow.geometry( str(int(puls_width)) +'x'+ str(int(puls_height)) )
    PulsatilityWindow.resizable(True,True)
    PulsatilityWindow.wm_title("Pulsatility analysis")
    
    # Scale figure size based on window dimensions
    puls_fig_size = max(3.0, 3.5 * dpi_scale)
    
    fig_puls = plt.figure(figsize=(puls_fig_size, puls_fig_size), dpi=100)
    ax_puls = fig_puls.add_subplot(211)
    ax_cums = fig_puls.add_subplot(212)
    canvas_puls = FigureCanvasTkAgg(fig_puls, master=PulsatilityWindow) 
    canvas_puls.draw()
    canvas_puls.get_tk_widget().grid(row=0, rowspan=2, column=0, columnspan=1, padx=0, pady=0, sticky='nsew')

    # Get flow data
    Flows_data = Flow_output.Flows
    if Flows_data.ndim > 1:
        Flows_arr = Flows_data[0]
    else:
        Flows_arr = Flows_data
    
    n_frames = len(Flows_arr)
    
    # Calculate basic flow parameters
    Qmean = np.mean(Flows_arr)
    Qmax = np.max(Flows_arr)
    Qmin = np.min(Flows_arr)
    
    # Calculate Pulsatility Index (PI)
    # PI = (Qmax - Qmin) / Qmean
    PI = (Qmax - Qmin) / Qmean if Qmean != 0 else 0
    
    # Create results panel (before plots so we can access HR entry)
    Puls_group = LabelFrame(PulsatilityWindow, text='Pulsatility Measures', borderwidth=2, relief='solid')
    Puls_group.grid(row=0, rowspan=1, column=1, columnspan=1, sticky='nw', padx=10, pady=20)
    
    # Heart rate input
    HR_text = Label(Puls_group, text="Heart rate [bpm]:", width=18, anchor='w')
    HR_text.grid(row=0, column=0, sticky='w', pady=2, padx=5)
    etr_HR = Entry(Puls_group, width=8)
    etr_HR.insert(END, '60')
    etr_HR.grid(row=0, column=1, sticky='e', pady=2, padx=5)
    
    # Separator
    separator1 = Label(Puls_group, text="─" * 25)
    separator1.grid(row=1, column=0, columnspan=2, pady=5)
    
    # Display Pulsatility Index (PI) - does not depend on HR
    PI_text = Label(Puls_group, text="PI = ", width=18, anchor='w')
    PI_text.grid(row=2, column=0, sticky='w', pady=2, padx=5)
    PI_str = "%.3f" % PI
    PI_val = Label(Puls_group, text=PI_str, width=10, anchor='e')
    PI_val.grid(row=2, column=1, sticky='e', pady=2, padx=5)
    
    # Display ΔV (Volume of Arterial Pulsatility) - placeholder, updated by button
    DeltaV_text = Label(Puls_group, text="\u0394V [ml] = ", width=18, anchor='w')
    DeltaV_text.grid(row=3, column=0, sticky='w', pady=2, padx=5)
    DeltaV_val = Label(Puls_group, text="---", width=10, anchor='e')
    DeltaV_val.grid(row=3, column=1, sticky='e', pady=2, padx=5)
    
    # Separator
    separator2 = Label(Puls_group, text="─" * 25)
    separator2.grid(row=4, column=0, columnspan=2, pady=5)
    
    # Display component values
    Qmax_text = Label(Puls_group, text="Qmax [ml/min] = ", width=18, anchor='w')
    Qmax_text.grid(row=5, column=0, sticky='w', pady=2, padx=5)
    Qmax_str = "%.3f" % Qmax
    Qmax_val = Label(Puls_group, text=Qmax_str, width=10, anchor='e')
    Qmax_val.grid(row=5, column=1, sticky='e', pady=2, padx=5)
    
    Qmin_text = Label(Puls_group, text="Qmin [ml/min] = ", width=18, anchor='w')
    Qmin_text.grid(row=6, column=0, sticky='w', pady=2, padx=5)
    Qmin_str = "%.3f" % Qmin
    Qmin_val = Label(Puls_group, text=Qmin_str, width=10, anchor='e')
    Qmin_val.grid(row=6, column=1, sticky='e', pady=2, padx=5)
    
    Qmean_text = Label(Puls_group, text="Qmean [ml/min] = ", width=18, anchor='w')
    Qmean_text.grid(row=7, column=0, sticky='w', pady=2, padx=5)
    Qmean_str = "%.3f" % Qmean
    Qmean_val = Label(Puls_group, text=Qmean_str, width=10, anchor='e')
    Qmean_val.grid(row=7, column=1, sticky='e', pady=2, padx=5)
    
    # Display dt (time per frame) - updated by button
    dt_text = Label(Puls_group, text="dt [s/frame] = ", width=18, anchor='w')
    dt_text.grid(row=8, column=0, sticky='w', pady=2, padx=5)
    dt_val = Label(Puls_group, text="---", width=10, anchor='e')
    dt_val.grid(row=8, column=1, sticky='e', pady=2, padx=5)
    
    # Add formula descriptions
    formula_group = LabelFrame(PulsatilityWindow, text='Formulas', borderwidth=2, relief='solid')
    formula_group.grid(row=1, column=1, sticky='nw', padx=10, pady=10)
    
    PI_formula = Label(formula_group, text="PI = (Qmax - Qmin) / Qmean", font=('Courier', 12), anchor='w')
    PI_formula.grid(row=0, column=0, sticky='w', pady=2, padx=5)
    
    DeltaV_formula = Label(formula_group, text="ΔV = max(∫[Q-Qmean]dt) - min(∫[Q-Qmean]dt)", font=('Courier', 12), anchor='w')
    DeltaV_formula.grid(row=1, column=0, sticky='w', pady=2, padx=5)
    
    
    # Function to calculate and update ΔV based on heart rate
    def update_DeltaV():
        try:
            HR = float(etr_HR.get())
            if HR <= 0:
                raise ValueError("HR must be positive")
        except ValueError:
            DeltaV_val['text'] = "Invalid HR"
            dt_val['text'] = "---"
            return
        
        # Calculate time step: dt = cardiac_cycle_duration / n_frames
        # cardiac_cycle_duration = 60 / HR (in seconds)
        cardiac_cycle_duration = 60.0 / HR  # seconds per beat

        
        dt = cardiac_cycle_duration / n_frames  # seconds per frame
        
        # Calculate cumulative volume with correct time step
        # Flow is in ml/s, dt is in seconds, so integral gives ml
        time_seconds = np.arange(n_frames) * dt
        Vmean = scipy.integrate.cumulative_trapezoid((Flows_arr - Qmean)/60, dx=dt, axis=-1, initial=0)
        
        # ΔV = max(V) - min(V) : Volume of arterial pulsatility in ml
        DeltaV = Vmean.max() - Vmean.min()
        
        # Update display
        DeltaV_val['text'] = "%.3f" % DeltaV
        dt_val['text'] = "%.4f" % dt
        
        # Update plots
        ax_puls.clear()
        ax_cums.clear()
        
        # Plot flow waveform with time in seconds
        ax_puls.tick_params(labelsize=4)
        ax_puls.set_ylabel('Flow [ml/min]', fontsize=5.0)
        ax_puls.set_xlabel('Time [s]', fontsize=5.0)
        ax_puls.plot(time_seconds, Flows_arr, '.r-', linewidth=0.8, markersize=3)
        
        # Add horizontal lines for Qmax, Qmin, Qmean
        ax_puls.axhline(y=Qmax, color='b', linestyle='--', linewidth=0.5, label=f'Qmax={Qmax:.2f}')
        ax_puls.axhline(y=Qmin, color='g', linestyle='--', linewidth=0.5, label=f'Qmin={Qmin:.2f}')
        ax_puls.axhline(y=Qmean, color='k', linestyle='-', linewidth=0.5, label=f'Qmean={Qmean:.2f}')
        ax_puls.legend(fontsize=4, loc='upper right')
        ax_puls.set_title(f'Flow waveform (HR={HR:.0f} bpm)', fontsize=6)
        
        # Plot cumulative volume
        ax_cums.tick_params(labelsize=4)
        ax_cums.set_ylabel('Cumulative volume [ml]', fontsize=5.0)
        ax_cums.set_xlabel('Time [s]', fontsize=5.0)
        ax_cums.plot(time_seconds, Vmean, '.r-', linewidth=0.8, markersize=3)
        
        # Add horizontal lines for max and min volume
        ax_cums.axhline(y=Vmean.min(), color='k', linestyle=':', linewidth=0.5)
        ax_cums.axhline(y=Vmean.max(), color='k', linestyle=':', linewidth=0.5)
        
        # Add annotation for DeltaV
        mid_idx = len(time_seconds) // 2
        mid_time = time_seconds[mid_idx]
        ax_cums.annotate('', xy=(mid_time, Vmean.max()), xytext=(mid_time, Vmean.min()),
                         arrowprops=dict(arrowstyle='<->', color='blue', lw=0.8))
        ax_cums.text(mid_time + dt, (Vmean.max() + Vmean.min())/2, f'ΔV={DeltaV:.2f} ml', fontsize=5, color='blue')
        ax_cums.set_title(f'Cumulative volume (ΔV={DeltaV:.3f} ml)', fontsize=6)
        
        fig_puls.tight_layout()
        canvas_puls.draw()
    
    # Calculate button
    CalcPuls_button = Button(master=Puls_group,
                             height=button_height,
                             width=button_width,
                             text="Calculate ΔV",
                             command=update_DeltaV)
    CalcPuls_button.grid(row=9, column=0, columnspan=2, sticky='ew', pady=10, padx=5)
    
    # Initial calculation with default HR
    update_DeltaV()


calcmenu = Menu(menubar, tearoff=1)
menubar.add_cascade(label="Analysis", menu=calcmenu)
calcmenu.add_command(label="Pulsatility", command=CalcPulsatility,accelerator="Control+p")
window.bind_all('<Control-Key-p>', func=Img_data.set_new_data )

# =============================================================================
# HELP MENU
# =============================================================================

def popup_help():
    """
    Display the About dialog with author and usage information.
    
    Keyboard shortcut: Ctrl+A
    """
    popup = tk.Toplevel(window)
    GUI_width = window.winfo_screenwidth()*0.35
    GUI_height = window.winfo_screenheight()*0.25
    popup.geometry( str(int(GUI_width)) +'x'+ str(int(GUI_height)) )
    popup.resizable(True,True)
    popup.wm_title("About me")
    help_str='GUI for calculating flow in blood vessel from PCM-images. \n'
    name_str='Mark B. Vestergaard \n mark.bitsch.vestergaard@regionh.dk, \n Functional Imaging Unit \n Department of Clinical Physiology and Nuclear Medicine \n Rigshospitalet, Glostrup, Denmark \n 2026.' 
    text_title = Label(popup,text=help_str,anchor="w", background='white')
    text_title.pack(side="top", fill="x", pady=10)
    text_name = Label(popup,text=name_str,justify="left",anchor="w")
    text_name.pack(side="top", fill="x", pady=10)
    
    B1 = Button(popup, text="Close", command = popup.destroy)
    B1.pack(side="top", pady=10)
    popup.mainloop()

analysismenu.add_command(label="About", command=popup_help, accelerator="Control+a")
window.bind_all('<Control-Key-a>', func=popup_help)


# =============================================================================
# MAIN EVENT LOOP
# =============================================================================
window.config(menu=menubar)  # Attach menu bar to window
window.mainloop()            # Start the GUI event loop
