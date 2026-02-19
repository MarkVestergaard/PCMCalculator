# PCMCalculator

**PCMCalculator** (**P**hase **C**ontrast **M**apping Calculator) is a Python-based GUI for analyzing 2D phase-contrast MRI data and quantifying blood flow in vessels (ml/min).
The software offers an intuitive graphical user interface that enables users to load, visualize, and analyze PCM data and calculate blood flow without any programming expertise. Images can be displayed with adjustable colormaps and intensity ranges to optimize vessel visibility. Region of interests (ROIs) covering the targeted vessels can be defined either through manual polygon delineation or generated semi-automatically using a built-in region-growing algorithm. From these ROIs, blood flow is calculated in quantitative units (ml/min), along with mean velocity (cm/s) and cross-sectional area (mm²). Results can be exported to CSV format for further statistical analysis, and ROI masks can be saved in NIfTI or NumPy archive formats for reproducibility. The software was developed so that personnel without a programming background could use it for their analyses, thus widening its usability across clinical and research settings.

![Interface](Interface.png)

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

---

## Usage

```bash
# Launch with file dialog
python PCMCalculator.py

# Load a Philips PAR/REC file directly
python PCMCalculator.py --img /path/to/file.PAR

# Load NIfTI files directly
python PCMCalculator.py --img_nii_vel velocity.nii --img_nii_mod modulus.nii --img_nii_mag magnitude.nii
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
9. Saving results in CSV format and ROIs in NIfTI, NPZ, and/or animated GIF format
---

Please see following reference for research studies where **PCMCalculator** has been used.

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

