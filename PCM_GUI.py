#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 21:33:05 2021

GUI for drawing ROI on PCM - images and calculate Flow.

@author: Mark B. Vestegaard, June 2021, mark.bitsch.vestergaard@regionh.dk 
"""

from tkinter import *
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import nibabel as nib
import argparse
import numpy as np
import numpy.matlib
import pandas as pd
#from functools import partial    
from matplotlib.widgets import PolygonSelector
from matplotlib.patches import Polygon 
from matplotlib.path import Path
from skimage import measure
#import csv
import os 
import scipy
from functools import partial


# parser for loading PAR/REC file
parser = argparse.ArgumentParser(description='Calculate flow in vessel')
parser.add_argument('--img', dest='img_input', help='Input PCM image')
parser.add_argument('--img_nii_vel', dest='img_input_nii_vel', help='Input velocity nifti image')
parser.add_argument('--img_nii_mod', dest='img_input_nii_mod', help='Input modolus nifti image')
parser.add_argument('--img_nii_thres', dest='img_input_nii_thres', help='Input thresholded phase nifti image')

input_files = parser.parse_args()

global nifti_files
nifti_files=None

# Load PCM data
if input_files.img_input==None: # Open dialog if no argument
    print('Select PCM file')
    raw_img_filename_tmp = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("Philips PAR","*.PAR"),("NIFTI","*.nii"),("all files","*.*")),multiple=True)        
    raw_img_filename=raw_img_filename_tmp[0]
    file_type = os.path.splitext(raw_img_filename)[1]
    if file_type=='.nii':
        nifti_files=True
        import json
else: # Set argument as filename
    raw_img_filename = input_files.img_input # 


# The main tkinter GUI window
window = Tk()
window.title('Calculate PCM flow')
GUI_width = window.winfo_screenwidth()*0.9
GUI_height = window.winfo_screenheight()*1
window.geometry( str(int(GUI_width)) +'x'+ str(int(GUI_height)) )
window.resizable(True,True)


# Function for specifying image data type
global Disp_image_str
Disp_image_str='vel'
global colormap_str
colormap_str='jet'
global imgFrame
imgFrame=1



# Class of PCM image data
class Img_data_parrec():
    def __init__(self, raw_img_filename):
        self.raw_img = nib.load(raw_img_filename)# Loads PCM par file nibabel load
        self.raw_img_filename=raw_img_filename
        self.img = nib.as_closest_canonical(self.raw_img) # Loads PCM par file
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        file_type = os.path.splitext(raw_img_filename)[1]
        self.Inv_image_indx=1
        
        if file_type=='.PAR':
            self.Venc=self.img.header.general_info.get('phase_enc_velocity')[2] # Find venc in header
            self.Venc_from_hdr=True
        else:
            self.Venc=100
            self.Venc_from_hdr=False

        self.thres_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==4)
        self.mod_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==2)
        self.phase_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==3,self.img.header.image_defs['scanning sequence']==4)
        
        self.Thres_image=(np.rot90(self.img.dataobj[:,:,0,self.thres_image_idx])/self.img.dataobj[:,:,0,self.phase_image_idx].max()*self.Venc*0.7)
        self.Mod_image=(np.rot90(self.img.dataobj[:,:,0,self.mod_image_idx])/self.img.dataobj[:,:,0,self.mod_image_idx].max()*self.Venc*0.9)
        self.Vel_image=np.rot90(self.img.dataobj[:,:,0,self.phase_image_idx])

    def set_new_data(self): #(self=Img_data): # Set new data
         
        self.raw_img_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("PAR files","*.PAR"),("all files","*.*")))        
        self.raw_img = nib.load(self.raw_img_filename)
        self.img = nib.as_closest_canonical(self.raw_img) # Loads PCM par file
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        self.Venc=self.img.header.general_info.get('phase_enc_velocity')[2]
        
        self.thres_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==4)
        self.mod_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==0,self.img.header.image_defs['scanning sequence']==2)
        self.phase_image_idx=np.logical_and(self.img.header.image_defs['image_type_mr']==3,self.img.header.image_defs['scanning sequence']==4)


        self.Thres_image=(np.rot90(self.img.dataobj[:,:,0,self.thres_image_idx])/self.img.dataobj[:,:,0,self.phase_image_idx].max()*self.Venc*0.2)
        self.Mod_image=(np.rot90(self.img.dataobj[:,:,0,self.mod_image_idx])/self.img.dataobj[:,:,0,self.mod_image_idx].max()*self.Venc*0.5)
        self.Vel_image=np.rot90(self.img.dataobj[:,:,0,self.phase_image_idx])
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


        

# NIFTI img class 

class Img_data_nii():
    def __init__(self, raw_img_filename):
        self.raw_img = nib.load(raw_img_filename[0])# Loads PCM par file nibabel load
        self.raw_img_filename=raw_img_filename[0]
        self.img = nib.as_closest_canonical(self.raw_img) # Loads PCM par file
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        fname_json_tmp = open(os.path.splitext(raw_img_filename[0])[0]+'.json')
        self.img.json_header = json.load(fname_json_tmp)     
        ImageType_a=self.img.json_header['ImageType'][2]
        
        self.raw_img_b = nib.load(raw_img_filename[1])# Loads PCM par file nibabel load
        self.raw_img_filename_b=raw_img_filename[1]
        self.img_b = nib.as_closest_canonical(self.raw_img_b) # Loads PCM par file
        self.nifti_file_b = nib.Nifti1Image(self.img_b.dataobj, self.img_b.affine, header=self.img_b.header)
        fname_json_tmp = open(os.path.splitext(raw_img_filename[1])[0]+'.json')
        self.img_b.json_header = json.load(fname_json_tmp)     
        ImageType_b=self.img_b.json_header['ImageType'][2]
        
        self.raw_img_c = nib.load(raw_img_filename[2])# Loads PCM par file nibabel load
        self.raw_img_filename_c=raw_img_filename[2]
        self.img_c = nib.as_closest_canonical(self.raw_img_c) # Loads PCM par file
        self.nifti_file_c = nib.Nifti1Image(self.img_c.dataobj, self.img_c.affine, header=self.img_c.header)
        fname_json_tmp = open(os.path.splitext(raw_img_filename[2])[0]+'.json')
        self.img_c.json_header = json.load(fname_json_tmp)     
        ImageType_c=self.img_c.json_header['ImageType'][2]
        
        
        self.ImageType=[ImageType_a, ImageType_b, ImageType_c]
        self.Inv_image_indx=1
        file_type = os.path.splitext(raw_img_filename[0])[1]
        if file_type=='.PAR':
            self.Venc=self.img.header.general_info.get('phase_enc_velocity')[2] # Find venc in header
            self.Venc_from_hdr=True
        else:
            self.Venc=100
            self.Venc_from_hdr=False
            
        
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
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file.dataobj), k=-1 ))/4096)*self.Venc
        elif n_tmp==1:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_b.dataobj), k=-1 ))/4096)*self.Venc
        elif n_tmp==2:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_c.dataobj), k=-1 ))/4096)*self.Venc

        
    
    
    def set_new_data(self): #(self=Img_data): # Set new data
        raw_img_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("NIFTI","*.nii"),("all files","*.*")),multiple=True)       
        self.raw_img = nib.load(raw_img_filename[0])# Loads PCM par file nibabel load
        self.raw_img_filename=raw_img_filename[0]
        self.img = nib.as_closest_canonical(self.raw_img) # Loads PCM par file
        self.nifti_file = nib.Nifti1Image(self.img.dataobj, self.img.affine, header=self.img.header)
        fname_json_tmp = open(os.path.splitext(raw_img_filename[0])[0]+'.json')
        self.img.json_header = json.load(fname_json_tmp)     
        ImageType_a=self.img.json_header['ImageType'][2]

        
        self.raw_img_b = nib.load(raw_img_filename[1])# Loads PCM par file nibabel load
        self.raw_img_filename_b=raw_img_filename[1]
        self.img_b = nib.as_closest_canonical(self.raw_img_b) # Loads PCM par file
        self.nifti_file_b = nib.Nifti1Image(self.img_b.dataobj, self.img_b.affine, header=self.img_b.header)
        fname_json_tmp = open(os.path.splitext(raw_img_filename[1])[0]+'.json')
        self.img_b.json_header = json.load(fname_json_tmp)     
        ImageType_b=self.img_b.json_header['ImageType'][2]
        
        self.raw_img_c = nib.load(raw_img_filename[2])# Loads PCM par file nibabel load
        self.raw_img_filename_c=raw_img_filename[2]
        self.img_c = nib.as_closest_canonical(self.raw_img_c) # Loads PCM par file
        self.nifti_file_c = nib.Nifti1Image(self.img_c.dataobj, self.img_c.affine, header=self.img_c.header)
        fname_json_tmp = open(os.path.splitext(raw_img_filename[2])[0]+'.json')
        self.img_c.json_header = json.load(fname_json_tmp)     
        ImageType_c=self.img_c.json_header['ImageType'][2]
        ImageType=[ImageType_a, ImageType_b, ImageType_c]
        
        file_type = os.path.splitext(raw_img_filename[0])[1]
        if file_type=='.PAR':
            self.Venc=self.img.header.general_info.get('phase_enc_velocity')[2] # Find venc in header
            self.Venc_from_hdr=True
        else:
            self.Venc=100
            self.Venc_from_hdr=False
            
        
           
        n_tmp=ImageType.index('M')
        if n_tmp==0:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file.dataobj/self.nifti_file.dataobj.max())*self.Venc*0.7, k=-1) 
        elif n_tmp==1:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file_b.dataobj/self.nifti_file_b.dataobj.max())*self.Venc*0.7, k=-1) 
        elif n_tmp==2:
            self.Mod_image=np.rot90(np.squeeze(self.nifti_file_c.dataobj/self.nifti_file_c.dataobj.max())*self.Venc*0.7,k=-1)   
            
        n_tmp=ImageType.index('MAG')
        if n_tmp==0:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file.dataobj/self.nifti_file.dataobj.max())*self.Venc*0.9 ,k=-1) 
        elif n_tmp==1:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file_b.dataobj/self.nifti_file_b.dataobj.max())*self.Venc*0.9 ,k=-1)
        elif n_tmp==2:
            self.Thres_image=np.rot90(np.squeeze(self.nifti_file_c.dataobj/self.nifti_file_c.dataobj.max())*self.Venc*0.9 ,k=-1)
            
        n_tmp=ImageType.index('P')
        if n_tmp==0:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file.dataobj), k=-1 ))/4096)*self.Venc
        elif n_tmp==1:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_b.dataobj), k=-1 ))/4096)*self.Venc
        elif n_tmp==2:
            self.Vel_image=((np.rot90(np.squeeze(self.nifti_file_c.dataobj), k=-1 ))/4096)*self.Venc

        global nifti_files
        nifti_files=True
        
        global Disp_image_str
        global colormap_str
        global imgFrame
        imgFrame=1
        slider_scale.set(1)
        slider_scale['to'] = 1
        slider_scale['from'] = Img_data.Vel_image.shape[2]
        #slider_scale.ax.set_xlim(slider_scale.valmin,slider_scale.valmax)
        change_image(Disp_image_str,colormap_str)
        
        save_str=[Img_data.raw_img_filename[0:-4]+'_Flow_data.csv']
        entry_save_filename.delete(0, 'end')
        entry_save_filename.insert(END, save_str[0])
        change_image(Disp_image_str,colormap_str)


if not nifti_files:
    Img_data=Img_data_parrec(raw_img_filename)        

if nifti_files:
    Img_data=Img_data_nii(raw_img_filename_tmp)  





def displayed_image(disp_name): # Returns Disp_image variable with shown data
    if disp_name == 'vel':
        return Img_data.Vel_image
    elif disp_name == 'mod':
        return Img_data.Mod_image
    elif disp_name == 'thres':
        return Img_data.Thres_image
    elif disp_name == 'ROI':
        return ROI_art.BWMask*100
    
Disp_image=displayed_image(Disp_image_str)


# Function for capturing change of image from topmenu or keyboard (could be changed to a class)
def change_image_type_str_vel(self=''):
    global Disp_image_str
    Disp_image_str='vel'
    change_image(Disp_image_str,colormap_str)

def change_image_type_str_mod(self=''):
    global Disp_image_str
    Disp_image_str='mod'
    change_image(Disp_image_str,colormap_str)

def change_image_type_str_thres(self=''):
    global Disp_image_str
    Disp_image_str='thres'
    change_image(Disp_image_str,colormap_str)
    
def change_image_type_str_ROI(self=''):
    global Disp_image_str
    Disp_image_str='ROI'
    change_image(Disp_image_str,colormap_str)




# Function for capturing change of colorbar from topmenu or keyboard (could be changed to a class)
def change_cmap_jet(self='jet'):
    global colormap_str
    colormap_str='jet'
    change_image(Disp_image_str,colormap_str)

def change_cmap_gray(self='gray'):
    global colormap_str
    colormap_str='gray'
    change_image(Disp_image_str,colormap_str)

def change_cmap_viridis(self='viridis'):
    global colormap_str
    colormap_str='viridis'
    change_image(Disp_image_str,colormap_str)
    
def popup_change_colorbar():
    popup_change_colorbar = Tk()
    popup_change_colorbar.title('Colorbar limits')
    GUI_width = window.winfo_screenwidth()*0.10
    GUI_height = window.winfo_screenheight()*0.22
    popup_change_colorbar.geometry( str(int(GUI_width)) +'x'+ str(int(GUI_height)) )
    popup_change_colorbar.resizable(True,True)
    
    #SinusROI_button_group = LabelFrame(window, text='Sinus ROI', borderwidth=2,relief='solid')
    #SinusROI_button_group.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='nw',padx=10,pady=20)

    Min_entry_text = Label(popup_change_colorbar, text="Min:",width = 7,padx=0)
    Min_entry_text.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='n',pady=10,padx=0)
    Min_entry = Entry(popup_change_colorbar,width=7)
    Min_entry.grid(row=1,rowspan=1,column=1,columnspan=1,sticky='nw',pady=10,padx=0)
    Min_entry.insert(END, '-5')
    
    Max_entry_text = Label(popup_change_colorbar, text="Max:",width = 7,padx=0)
    Max_entry_text.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='n',pady=0,padx=0)
    Max_entry = Entry(popup_change_colorbar,width=7)
    Max_entry.grid(row=0,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
    Max_entry.insert(END, '100')
    
    def update_colorbar():
        climits.lims=[float(Min_entry.get()),float(Max_entry.get())]
        change_image(Disp_image_str,colormap_str)
        
    
    UpdateCB_button = Button(master = popup_change_colorbar,
                      height = 2,
                      width = 9,
                     text = "Update Colorbar", command=update_colorbar)
    UpdateCB_button.grid(row=2,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)



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
submenu.add_radiobutton(label="Threshold" ,command  = change_image_type_str_thres,accelerator="Control+3")
submenu.add_radiobutton(label="ROI mask" ,command  = change_image_type_str_ROI,accelerator="Control+4")

window.bind_all('<Control-Key-1>', func=change_image_type_str_vel ) # Binds keyboard shortcut to topbar 
window.bind_all('<Control-Key-2>', func=change_image_type_str_mod )
window.bind_all('<Control-Key-3>', func=change_image_type_str_thres )
window.bind_all('<Control-Key-4>', func=change_image_type_str_ROI )

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




# Create grid for GUI
Grid.rowconfigure(window, 0, weight=0)
Grid.columnconfigure(window, 0, weight=0)

    
# Create figure frame for displaying MRI image
fig = plt.figure(figsize=(5.5, 5.5), dpi=100)
ax = fig.add_subplot(111)

class climits():
    lims=[-5,100]

global pcm_plot
pcm_plot=ax.imshow(Disp_image[:,:,imgFrame-1],cmap=plt.get_cmap(colormap_str),vmin=climits.lims[0],vmax=climits.lims[1], interpolation='none')
ax.set_xticks([]) # Removes x axis ticks
ax.set_yticks([]) # Removes x axis ticks
fig.tight_layout() # Tight layout

canvas = FigureCanvasTkAgg(fig, master=window)  # Draws the figure in the tkinter GUI window
canvas.draw()
canvas.get_tk_widget().grid(row=1,rowspan=1,column=1,columnspan=1,padx=0,pady=0) # place in grid

# Create line data for plotting ROI data
ROI_line_plot, = ax.plot([], [], '.w-')
fig_flow = plt.figure(figsize=(3,3), dpi=110)
ax_flow = fig_flow.add_subplot(111)


# Create line plot for flow data
ax_flow.set_position([0.2, 0.15, 0.7, 0.7])
ax_flow.tick_params(labelsize=7)
ax_flow.set_ylabel('Flow [ml/min]', fontsize = 8.0) # Y label
ax_flow.set_xlabel('Index', fontsize = 8.0) # Y label
Flow_line_plot, = ax_flow.plot([], [], '.r-')

canvas_flow = FigureCanvasTkAgg(fig_flow, master=window)  # Draws the figure in the tkinter GUI window
canvas_flow.draw()
canvas_flow.get_tk_widget().grid(row=1,rowspan=1,column=5,columnspan=1,padx=20,pady=50) # place in grid


# Create frame for the navigation toolbar
ToolbarFrame = Frame(window)
ToolbarFrame.grid(row=2, column=1,rowspan=1,columnspan=1,padx=0,pady=0,sticky='nw')
toobar = NavigationToolbar2Tk(canvas, ToolbarFrame)

# Function for changing image
def change_image(image_str,cmp_str): # Changed displayed image
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
    #breakpoint()
    return Disp_image_str


# Function for changing the image data and ROI data
def update_image(self):
    global imgFrame
    imgFrame=int(self)
    Disp_image=displayed_image(Disp_image_str)
    pcm_plot.set_data(Disp_image[:,:,int(imgFrame-1)])
    
    ROI_as_array = np.array(ROI_art.polygon[int(imgFrame-1)])
    if ROI_as_array.any() == None:
        ROI_line_plot.set_ydata([])
        ROI_line_plot.set_xdata([])
    else:
        ROI_as_array_tmp=np.append(ROI_as_array[:,:], ROI_as_array[0,:]).reshape(ROI_as_array.shape[0]+1,ROI_as_array.shape[1])
        ROI_line_plot.set_ydata(ROI_as_array_tmp[:,1])
        ROI_line_plot.set_xdata(ROI_as_array_tmp[:,0])
        
    canvas.draw()    
    return imgFrame


# Create slider for changing frame.
slider_scale = Scale(window, to=1, from_=Img_data.Vel_image.shape[2], width=20,length=400 ,command=update_image)
slider_scale.grid(row=1,rowspan=1,column=4,columnspan=1,padx=0,pady=0,sticky='w')


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
headline_text = Label(window, text="Draw ROI to calculate blood flow in the blood vessel")
headline_text.grid(row=0,rowspan=1,column=0,columnspan=2,sticky='nw')


# Class of ROI data
class ROI_art(object):
    polygon=[None] * Img_data.Vel_image.shape[2]
    BWMask=Img_data.Vel_image*False
    flag=[0]*(Img_data.Vel_image.shape[2])

    def set_polygon(verts):
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
        ROI_art.polygon=Loaded_ROI_data['PCMROI_poly'].tolist()
        ROI_art.BWMask=Loaded_ROI_data['PCMROI_BWMask']
        ROI_art.flag=Loaded_ROI_data['PCMROI_flag']
        global imgFrame
        update_image(imgFrame)
        

# Class for drawing ROI as a polygon
global PS
PS=None
class ROIPolygon(object):
    def __init__(self, ax, row, col):
        self.canvas = ax.figure.canvas
        global PS
        PS = PolygonSelector(ax,
                                    self.onselect,
                                    lineprops = dict(color = 'm', alpha = 1),
                                    markerprops = dict(mec = 'm', mfc = 'm', alpha = 1),vertex_select_radius=10)
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

    
# Function for initiating ROI polygon 
def add_roi_func():
    global PS
    if PS!=None:
        PS.set_visible(False)
        PS.set_active(False)
    roi_poly=ROIPolygon(pcm_plot.axes, Img_data.Vel_image.shape[0], Img_data.Vel_image.shape[2])
    UpdateROI_button.configure(text='Stop ROI edit')
    
# Copies ROIs to all frames
def copy_roi_func():
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

# Activates polygon drawing for updating ROI polygon 
def update_roi_func():
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


# Function for automatic delineation of ROI by region growing algorithm
class RegGrow():

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


# Create buttons for ROI controls 
ROI_button_group = LabelFrame(window, text='ROI analysis', borderwidth=2,relief='solid')
ROI_button_group.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='nw',padx=10,pady=20)

# Add new ROI
AddROI_button = Button(master = ROI_button_group,
                     height = 2,
                     width = 10,
                    text = "Add ROI", command=add_roi_func)
AddROI_button.grid(row=0,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)

# Update ROI
UpdateROI_button = Button(master = ROI_button_group,
                     height = 2,
                     width = 10,
                    text = "Edit ROI", command=update_roi_func)
UpdateROI_button.grid(row=1,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)

# Copy ROI to remaning frames 
CopyROI_button = Button(master = ROI_button_group,
                     height = 2,
                     width = 10,
                    text = "Copy ROI to \n all frames", command=copy_roi_func)
CopyROI_button.grid(row=3,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)

# Call region growing algorithm
AutoROI_button = Button(master = ROI_button_group,
                     height = 3,
                     width = 10,
                    text = "Automatic \n delineation", command=RegGrow.Calc_Auto_ROI)
AutoROI_button.grid(row=4,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)
Max_dist_text = Label(ROI_button_group, text="Max. dist:",width = 7,padx=0)
Max_dist_text.grid(row=5,rowspan=1,column=0,columnspan=1,sticky='n',pady=0,padx=0)

# Entries for region growing algorithm 
e1 = Entry(ROI_button_group,width=4)
e1.grid(row=5,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
e1.insert(END, '8')

Threshold_text = Label(ROI_button_group, text="Threshold:")
Threshold_text.grid(row=6,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=10)
e2 = Entry(ROI_button_group,width=4)
e2.grid(row=6,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
e2.insert(END, '10')

# Function for calculating flow
if nifti_files:
    vox_size=np.prod(Img_data.img.header.get_zooms()[0:2])
else: 
    vox_size=np.prod(Img_data.img.header.get_zooms()[0:2])

def inverse_vel_image():    
    Img_data.Inv_image_indx=Img_data.Inv_image_indx*-1
    if nifti_files:
        Img_data.Venc=float(Venc_ent.get())
    
        n_tmp=Img_data.ImageType.index('P')
        if n_tmp==0:
            Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file.dataobj), k=-1 ))/4096)*Img_data.Venc*Img_data.Inv_image_indx
        elif n_tmp==1:
            Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_b.dataobj), k=-1 ))/4096)*Img_data.Venc*Img_data.Inv_image_indx
        elif n_tmp==2:
            Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_c.dataobj), k=-1 ))/4096)*Img_data.Venc*Img_data.Inv_image_indx
    else: 
        Img_data.Vel_image=np.rot90(Img_data.img.dataobj[:,:,0,Img_data.phase_image_idx])*Img_data.Inv_image_indx

    change_image(Disp_image_str,colormap_str)


Inv_image_button_group = LabelFrame(window, borderwidth=2,relief='solid')
Inv_image_button_group.grid(row=1,rowspan=1,column=5,columnspan=2,sticky='sw',padx=0,pady=0)
Inv_image_button = Button(master = Inv_image_button_group,
                      height = 2,
                      width = 10,
                      text = "Inv. Vel. image", command=inverse_vel_image)
Inv_image_button.grid(row=0,rowspan=1,column=0, columnspan=1,sticky='nw',pady=0,padx=0)  



class Flow_output:
    Vel_image_tmp=ROI_art.BWMask*Img_data.Vel_image
    Flows= np.empty((1,Vel_image_tmp.shape[2],))
    Flows[:] = np.nan
    Velocity= np.empty((1,Vel_image_tmp.shape[2],))
    Velocity[:] = np.nan
    flow_str=''
    CS_area= np.empty((1,Vel_image_tmp.shape[2],))
    CS_area[:] = np.nan
    def Calc_flow():
        
        if nifti_files:
            Img_data.Venc=float(Venc_ent.get())
            n_tmp=Img_data.ImageType.index('P')
            if n_tmp==0:
                Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file.dataobj), k=-1 ))/4096)*Img_data.Venc*Img_data.Inv_image_indx
            elif n_tmp==1:
                Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_b.dataobj), k=-1 ))/4096)*Img_data.Venc*Img_data.Inv_image_indx
            elif n_tmp==2:
                Img_data.Vel_image=((np.rot90(np.squeeze(Img_data.nifti_file_c.dataobj), k=-1 ))/4096)*Img_data.Venc*Img_data.Inv_image_indx
            change_image(Disp_image_str,colormap_str)

        Flow_output.Vel_image_tmp=ROI_art.BWMask[:,:,range(0,Img_data.Vel_image.shape[2])]*Img_data.Vel_image
        Flow_output.Flows=np.empty((1,Flow_output.Vel_image_tmp.shape[2],))
        Flow_output.Velocity=np.empty((1,Flow_output.Vel_image_tmp.shape[2],))
        Flow_output.CS_area=np.empty((1,Flow_output.Vel_image_tmp.shape[2],))
        for x in range(0,Img_data.Vel_image.shape[2]):
            ROI_indx, ROI_indy= np.nonzero(Flow_output.Vel_image_tmp[:,:,x])
            Velocity_values=Flow_output.Vel_image_tmp[ROI_indx,ROI_indy,x]
            Flow_output.Flows[0][x]=Velocity_values.mean()*( (ROI_indx.shape[0]*vox_size) /1e2)*60  
            Flow_output.Velocity[0][x]=Velocity_values.mean()
            Flow_output.CS_area[0][x]=ROI_indx.shape[0]*vox_size
        print('Flow:')
        print(Flow_output.Flows)
        Flow_line_plot.set_ydata(Flow_output.Flows)
        Flow_line_plot.set_xdata(range(1,Img_data.Vel_image.shape[2]+1))
        ax_flow.set_xlim([1,Img_data.Vel_image.shape[2]+1])
        ax_flow.set_ylim([np.min(Flow_output.Flows)*0.9,np.max(Flow_output.Flows)*1.1])
        canvas_flow.draw()    
        Flow_output.flow_str="%5.2f" % Flow_output.Flows.mean()
        Flow_text_str.configure(text=Flow_output.flow_str)

# Function for saving output
class save_ouput_data:
    #output_file=args.img_input[0:-4]+'_Flow_data.csv'
    
    def save_data(self=''):
        output_file=entry_save_filename.get()

        d={'Flow': Flow_output.Flows[0].tolist(), 'Velocity': Flow_output.Velocity[0].tolist(), 'CSarea':Flow_output.CS_area[0].tolist()}
        df = pd.DataFrame(data=d)
        df.to_csv(output_file)
        Data_saved_txt.configure(text='Data saved:'+output_file)
        
        print('Data saved:'+output_file)
        ROI_filename=os.path.splitext(output_file)[0]+'_ROIs' # Also save ROI as npz data
        np.savez(ROI_filename, PCMROI_poly=np.array(ROI_art.polygon,dtype='object'), PCMROI_BWMask=ROI_art.BWMask, PCMROI_flag=ROI_art.flag)



# Button groups for saving data
Save_button_group = LabelFrame(window, text='Save data', borderwidth=2,relief='solid')
Save_button_group.grid(row=3,rowspan=1,column=0,columnspan=5,sticky='nw',padx=10,pady=2)
Save_button = Button(master = Save_button_group,
                     height = 2,
                     width = 8,
                    text = "Save data", command=save_ouput_data.save_data)
Save_button.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='sw',pady=2,padx=10)

# Entry for output file ROI_button_group = LabelFrame(window, text='ROI analysis', borderwidth=2,relief='solid')



#save_str=StringVar()
entry_save_filename=Entry(Save_button_group,width=80)
entry_save_filename.grid(row=2,rowspan=1,column=0,columnspan=4,sticky='nw',pady=0,padx=0)
save_str=[Img_data.raw_img_filename[0:-4]+'_Flow_data.csv']
entry_save_filename.insert(END, save_str[0])

#save_str.set(Img_data.raw_img_filename[0:-4]+'_Flow_data.csv')

Data_saved_txt = Label(Save_button_group)
Data_saved_txt.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='se',pady=0,padx=0)


    
# Button group for calculating flow
calc_button_group = LabelFrame(window, text='Calculate flow test', borderwidth=2,relief='solid')
calc_button_group.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='sw',padx=10,pady=10)

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
                      height = 2,
                      width = 10,
                      text = "Calculate flow", command=Flow_output.Calc_flow)

CalcFlow_button.grid(row=1,rowspan=1,column=0,columnspan=2,sticky='nw',pady=2,padx=10)
Flow_text = Label(calc_button_group, text="Mean flow:")
Flow_text.grid(row=2,rowspan=1,column=0,columnspan=1,sticky='nw')

Flow_text_str = Label(calc_button_group)
Flow_text_str.grid(row=2,rowspan=1,column=1,columnspan=1,sticky='ne')


#VENC_text = Label(calc_button_group, text="VENC:")
#VENC_text.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=10)
#e_venc = Entry(calc_button_group,width=4)
#e_venc.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=120)
#e_venc.insert(END, str(Img_data.Venc))



def load_ROI_file(self=''): # Load ROI as npz data
    ROI_filename = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("ROI files","*.npz"),("all files","*.*")))
    Loaded_ROI = np.load(ROI_filename,allow_pickle=True)
    ROI_art.loaded_roi_from_file(Loaded_ROI)


# Load ROI shortcut
analysismenu.add_command(label="Load ROI file", command=load_ROI_file, accelerator="Control+r")
window.bind_all('<Control-Key-r>', func=load_ROI_file )

# Save data shortcut
analysismenu.add_command(label="Save flow", command=save_ouput_data.save_data, accelerator="Control+s")
window.bind_all('<Control-Key-s>', func=save_ouput_data.save_data)


# Help menu 

def popup_help():
    popup = Tk()
    GUI_width = window.winfo_screenwidth()*0.35
    GUI_height = window.winfo_screenheight()*0.25
    popup.geometry( str(int(GUI_width)) +'x'+ str(int(GUI_height)) )
    popup.resizable(True,True)
    popup.wm_title("About me")
    help_str='GUI for calculating flow in blood vessel from PCM-images. \n Useable for PAR/REC philips file.'
    name_str='Mark B. Vestergaard \n mark.bitsch.vestergaard@regionh.dk, \n Functional Imaging Unit \n Department of Clinical Physiology, Nuclear Medicine and PET \n Rigshospitalet, Glostrup, Denmark \n July 2021.' 
    text_title = Label(popup,text=help_str,anchor="w", background='white')
    text_title.pack(side="top", fill="x", pady=10)
    text_name = Label(popup,text=name_str,justify="left",anchor="w")
    text_name.pack(side="top", fill="x", pady=10)
    
    B1 = Button(popup, text="Close", command = popup.destroy)
    B1.pack(side="top", pady=10)
    popup.mainloop()

analysismenu.add_command(label="About", command=popup_help, accelerator="Control+a")
window.bind_all('<Control-Key-a>', func=popup_help)


# For calculating pulsatility

def CalcPulsatility():
    PulsatilityWindow = Tk()
    GUI_width = window.winfo_screenwidth()*0.65
    GUI_height = window.winfo_screenheight()*0.65
    PulsatilityWindow.geometry( str(int(GUI_width)) +'x'+ str(int(GUI_height)) )
    PulsatilityWindow.resizable(True,True)
    PulsatilityWindow.wm_title("Calculate pulsatility")
    
    cardiac_phs=Img_data.img.header.general_info['max_cardiac_phases']
    cardiac_dyns=Img_data.img.header.general_info['max_dynamics']
    Flows_resh=Flow_output.Flows.reshape(cardiac_dyns,cardiac_phs)/60

     
    fig_puls = plt.figure(figsize=(5.5, 5.5), dpi=100)
    ax_puls = fig_puls.add_subplot(211)
    ax_cums = fig_puls.add_subplot(212)
    canvas_puls = FigureCanvasTkAgg(fig_puls, master=PulsatilityWindow) 
    canvas_puls.draw()
    canvas_puls.get_tk_widget().grid(row=0,rowspan=1,column=0,columnspan=1,padx=0,pady=0) # place in grid

    #ax_puls.set_position([0.2, 0.15, 0.7, 0.7])
    ax_flow.tick_params(labelsize=7)
    ax_puls.set_ylabel('Flow [ml/s]', fontsize = 8.0) # Y label
    ax_puls.set_xlabel('Index', fontsize = 8.0) # Y label
    Flows_mean=np.mean(Flows_resh,0)
    Qmean=np.mean(Flows_mean)
    puls_line_plot = ax_puls.plot(np.transpose(Flows_resh-Qmean), 'r',linewidth=0.5)

    mean_puls_line_plot = ax_puls.plot(Flows_mean-Qmean, '.k-')
    #puls_line_plot.set_ydata(Flows_resh)
    #puls_line_plot.set_xdata(range(cardiac_phs))
    #ax_puls.set_xlim([0,Img_data.Vel_image.shape[2]])
    #ax_puls.set_ylim([np.min(Flow_output.Flows)*0.9,np.max(Flow_output.Flows)*1.1])
    
    
    time=range(cardiac_phs)

    
    Vmean=scipy.integrate.cumtrapz(Flows_mean-Qmean, x=time, dx=1.0, axis=-1,initial=Flows_mean[0]-Qmean)/cardiac_phs
    ax_cums.tick_params(labelsize=7)
    ax_cums.set_ylabel('Cumulative flows', fontsize = 8.0) # Y label
    ax_cums.set_xlabel('Index', fontsize = 8.0) # Y label
    cums_line_plot = ax_cums.plot(time, Vmean, '.r-') 
 
    min_Vmean_line=ax_cums.plot( [time[0], time[-1]], [Vmean.min(), Vmean.min()], '.k:')
    max_Vmean_line=ax_cums.plot([time[0], time[-1]], [Vmean.max(),Vmean.max()], '.k:') 
    
    canvas_puls.draw()  
    
    DeltaV=Vmean.max()-Vmean.min()
    #cums_line_plot.set_ydata(Flow_output.Flows+10)
    #cums_line_plot.set_xdata(range(Img_data.Vel_image.shape[2]))
    #ax_cums.set_xlim([0,Img_data.Vel_image.shape[2]])
    #ax_cums.set_ylim([np.min(Flow_output.Flows)*0.9,np.max(Flow_output.Flows)*1.1])
    Puls_group = LabelFrame(PulsatilityWindow, text='Pulsatility', borderwidth=2,relief='solid')
    Puls_group.grid(row=0,rowspan=1,column=1,columnspan=1,sticky='nw',padx=10,pady=20)
    
    DeltaV_text=Label(Puls_group, text="\u0394 V = ",width = 15,padx=0)
    DeltaV_text.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=0)
    DeltaV_str="%5.5f" % DeltaV
    DeltaV_val=Label(Puls_group, text=DeltaV_str)
    DeltaV_val.grid(row=0,rowspan=1,column=0,columnspan=1,sticky='ne',pady=0,padx=0)
    

    HR_text = Label(Puls_group, text="Heart rate:",width = 15,padx=0)
    HR_text.grid(row=1,rowspan=1,column=0,columnspan=1,sticky='n',pady=0,padx=0)

# Entries for region growing algorithm 
    etr_HR = Entry(Puls_group,width=4)
    etr_HR.insert(END,'60')
    etr_HR.grid(row=1,rowspan=1,column=1,columnspan=1,sticky='nw',pady=0,padx=0)
    #sv.trace_add("write", calc_DeltaQ())
    class calc_DeltaQ:
        HR=etr_HR.get()
        etr_HR.delete(0, last=END )
        etr_HR.insert(END,HR)
        
        omega=2*3.14159*int(HR)/60
        DeltaQ=omega*DeltaV
    
        DeltaQ_text=Label(Puls_group, text="\u0394 Q = ",width = 15,padx=0)
        DeltaQ_text.grid(row=2,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=0)
        DeltaQ_str="%5.5f" % DeltaQ
        DeltaQ_val=Label(Puls_group, text=DeltaQ_str)
        DeltaQ_val.grid(row=2,rowspan=1,column=0,columnspan=1,sticky='ne',pady=0,padx=0)

        PP=DeltaQ/Qmean
    
        PP_text=Label(Puls_group, text="PP= ",width = 15,padx=0)
        PP_text.grid(row=3,rowspan=1,column=0,columnspan=1,sticky='nw',pady=0,padx=0)
        PP_str="%5.5f" % PP
        PP_val=Label(Puls_group, text=PP_str)
        PP_val.grid(row=3,rowspan=1,column=0,columnspan=1,sticky='ne',pady=0,padx=0)   
       
        def update_QDelta():
            HR=etr_HR.get()
            etr_HR.delete(0, last=END )
            etr_HR.insert(END,HR)
            calc_DeltaQ.omega=2*3.14159*int(HR)/60
            calc_DeltaQ.DeltaQ=calc_DeltaQ.omega*DeltaV
            print(calc_DeltaQ.DeltaQ)
            
            calc_DeltaQ.DeltaQ_str="%5.5f" % calc_DeltaQ.DeltaQ
            calc_DeltaQ.DeltaQ_val['text'] = calc_DeltaQ.DeltaQ_str
            calc_DeltaQ.PP=calc_DeltaQ.DeltaQ/Qmean
            calc_DeltaQ.PP_str="%5.5f" % calc_DeltaQ.PP
            calc_DeltaQ.PP_val['text'] = calc_DeltaQ.PP_str
          



    CalcPuls_button = Button(master=Puls_group,
                      height = 2,
                      width = 10,
                      text = "Calc Puls", command=calc_DeltaQ.update_QDelta)
    CalcPuls_button.grid(row=4,rowspan=1,column=0,columnspan=1,sticky='ne',pady=0,padx=0)

         




window.config(menu=menubar) # Insert topmenu
window.mainloop()