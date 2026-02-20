# PCMCalculator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18712355.svg)](https://doi.org/10.5281/zenodo.18712355)

**PCMCalculator** (**P**hase **C**ontrast **M**apping Calculator) is a Python-based GUI for analyzing 2D phase-contrast MRI data and quantifying blood flow in vessels (ml/min).
The software offers an intuitive graphical user interface that enables users to load, visualize, and analyze PCM data and calculate blood flow without any programming expertise. Images can be displayed with adjustable colormaps and intensity ranges to optimize vessel visibility. Region of interests (ROIs) covering the targeted vessels can be defined either through manual polygon delineation or generated semi-automatically using a built-in region-growing algorithm. From these ROIs, blood flow is calculated in quantitative units (ml/min), along with mean velocity (cm/s) and cross-sectional area (mm²). Results can be exported to CSV format for further statistical analysis, and ROI masks can be saved in NIfTI or NumPy archive formats for reproducibility. The software was developed so that personnel without a programming background could use it for their analyses, thus widening its usability across clinical and research settings.


---

## Features

- Supports Philips PAR/REC and NIfTI file formats
- Manual and semi-automatic (region-growing) ROI delineation
- Blood flow quantification (ml/min), mean velocity (cm/s), and cross-sectional area (mm²)
- Pulsatility analysis (PI and ΔV) for cardiac-gated acquisitions
- Export to CSV, NIfTI, NPZ, and animated GIF

---

## Installation

Requires Python 3.13+. Clone the repository and install dependencies:

```bash
git clone https://github.com/MarkVestergaard/PCMCalculator.git
cd PCMCalculator
pip install -r requirements.txt
```

The required packages are: tkinter (included with standard Python), matplotlib, nibabel, numpy, pandas, scipy, scikit-image, and Pillow.

---

## Usage

`PCMCalculator` is launched from the command line:

```bash
cd /path/to/PCMCalculator
python PCMCalculator.py
```

When launched without arguments, a file dialog opens for the user to select the input image file. Alternatively, input files can be specified directly via command-line arguments. For Philips PAR/REC files:

```bash
python PCMCalculator.py --img /path/to/file.PAR
```

For NIfTI files (converted from DICOM using dcm2niix), three files corresponding to the velocity (phase), modulus, and magnitude images must be provided:

```bash
python PCMCalculator.py \
    --img_nii_vel velocity.nii \
    --img_nii_mod modulus.nii \
    --img_nii_mag magnitude.nii
```

For NIfTI input, JSON sidecar files produced by [dcm2niix](https://github.com/rordenlab/dcm2niix) must be present alongside the NIfTI files.

The analysis workflow consists of:

1. Loading magnitude and phase (velocity) images from the phase-contrast acquisition in PAR/REC or NIfTI format
2. Adjusting visualization settings including colormap (jet, grayscale, or viridis) and intensity range to optimise vessel visibility
3. If needed, inverting the velocity image to ensure a positive flow direction in the target vessel
4. Defining ROIs around vessels of interest either by manual polygon delineation or by using the semi-automatic region-growing algorithm
5. Copying the ROI to all frames, then editing individual frames as needed to ensure accurate delineation throughout each frame
6. Computing mean velocity and cross-sectional area within each ROI for each frame
7. Calculating flow and reviewing the flow waveform and mean flow value
8. Optionally performing pulsatility analysis to obtain the Pulsatility Index and ΔV
9. Saving results in CSV format and ROIs in NIfTI, NPZ, and/or GIF format
---
 ## Interface 

 ![Interface](Interface.png)
`PCMCalculator` provides a graphical user interface that enables users to load and visualize the MRI data, delineate vessels using ROI tools, and calculate blood flow. The interface is organized into three main panels: a control panel on the left for ROI operations and flow calculation, a central image display panel showing the phase-contrast data with ROI overlays, and a right-hand panel displaying the calculated flow waveform across frames. The displayed image type can be switched between velocity, modulus, magnitude, ROI mask, and mean velocity map via the top menu (Image → Change Image Type) or by keyboard shortcuts (Ctrl+1 through Ctrl+5). The colormap can be changed between jet, grayscale, and viridis via the menu (Image → Change Colorbar) or by keyboard shortcuts (Ctrl+Q, Ctrl+W, and Ctrl+E). Colorbar intensity limits can also be manually adjusted. If the data contains multiple frames, these can be navigated using either the arrow keys or the slider. New data can be loaded via the top menu (File → Load New PAR/REC File or File → Load New NIfTI File).

`PCMCalculator` provides two methods for ROI delineation. Manual ROIs are drawn by selecting **"Add ROI"** in the ROI analysis button group, which activates the polygon drawing tool. Polygons are directly drawn on the image using an interactive polygon selector. After drawing, the ROI can be toggled between editable and locked states using the **"Edit ROI"** button. An ROI can be drawn on one frame and copied to all frames, and can be edited on each individual frame to ensure accurate delineation throughout the cardiac cycle. Alternatively, the semi-automatic region-growing algorithm allows users to place a seed point approximately in the center of the vessel, after which the algorithm grows outward to include neighbouring voxels whose velocity values exceed a user-defined threshold (default: 10 cm/s). The growth is additionally constrained by a maximum distance from the seed point (default: 8 pixels) to prevent the region from extending beyond the vessel boundary. For multi-frame data, the region-growing algorithm automatically tracks the vessel across all frames by updating the seed position to the voxel with the highest velocity within the current region.

Once the ROI has been defined, the flow in each frame can be calculated by pressing **"Calculate Flow"**. The resulting flow waveform is displayed in the rightmost panel. Results can be saved as a CSV file containing flow (ml/min), mean velocity (cm/s), and cross-sectional area (mm²) for each frame using **"Save Data"**. ROI masks can optionally be saved in NIfTI, NPZ, and/or animated GIF format.

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| Ctrl+1–5 | Switch image display (velocity, modulus, magnitude, ROI, mean velocity) |
| Ctrl+Q/W/E | Change colormap (jet, grayscale, viridis) |
| Ctrl+N | Load new PAR/REC file |
| Ctrl+M | Load new NIfTI file |
| Ctrl+R | Load saved ROI file |
| Ctrl+S | Save flow data |
| Ctrl+P | Open pulsatility analysis |
| Up/Down arrows | Navigate through frames |

## References
Please see the following references for research studies in which **PCMCalculator** has been used.

Reproducibility of cerebral blood flow, oxygen metabolism, and lactate and N-acetyl-aspartate concentrations measured using magnetic resonance imaging and spectroscopy
Signe Sloth Madsen, Ulrich Lindberg, Sohail Asghar, Karsten Skovgaard Olsen, Kirsten Møller, Henrik Bo Wiberg Larsson, Mark Bitsch Vestergaard. 
Frontiers in Physiology. 2023: 14:1213352.

Cerebral metabolic rate of oxygen is correlated to treatment effect of electroconvulsive therapy in patients with depression
Christoffer Cramer Lundsgaard, André Beyer Mathiassen, Henrik Bo Wiberg Larsson, Poul Videbech, Krzysztof Gbyl, Mark Bitsch Vestergaard
Brain Stimulation. 2025, 18(5): 1470-1478

Glucose-dependent insulinotropic polypeptide is involved in postprandial regulation of splanchnic blood supply. 
Rasmus S Rasmussen, Ludvig S Langberg, Frederikke Østergaard, Sophie W Nielsen, Mark B Vestergaard, Kirsa Skov-Jeppesen, Bolette Hartmann, Helle Hjorth Johannesen, Jens J Holst, Bryan Haddock, Henrik BW Larsson, Mette M Rosenkilde, Ali Asmar, Ulrik B Andersen, Lærke S Gasbjerg.
Diabetes. 2025, 74(8):1355-1366

## Author

**Mark B. Vestergaard**  
Functional Imaging Unit, Department of Clinical Physiology and Nuclear Medicine  
Copenhagen University Hospital Rigshospitalet, Glostrup, Denmark  
mark.bitsch.vestergaard@regionh.dk
markbvestergaard@gmail.com

---

