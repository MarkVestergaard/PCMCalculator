# PCMCalculator

**PCMCalculator** (**P**hase **C**ontrast **M**apping Calculator) is a Python-based GUI for analyzing 2D phase-contrast MRI data and quantifying blood flow in vessels (ml/min). It is designed for researchers and clinicians without programming expertise.

![Interface](paper/Figures/Figure1.png)

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

---


## Author

**Mark B. Vestergaard**  
Functional Imaging Unit, Department of Clinical Physiology and Nuclear Medicine  
Copenhagen University Hospital Rigshospitalet, Glostrup, Denmark  
mark.bitsch.vestergaard@regionh.dk
markbvestergaard@gmail.com

---

