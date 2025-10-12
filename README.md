# Reaction Path Refinement Algorithm (RPRA)

RPRA is a comprehensive computational tool designed to analyze, identify, and visualize reaction pathways from GAUSSIAN trajectories. The tool specializes in processing XYZ molecular trajectory files to detect energy peaks and valleys, calculate activation barriers, and identify low-barrier reaction pathways.

## General Description

RPRA performs automated analysis of chemical reaction pathways by:

- **Energy Profile Analysis**: Extracts energy data from XYZ trajectory files and identifies peaks (transition states) and valleys (intermediates) in the energy landscape
- **Barrier Calculation**: Computes activation barriers between consecutive energy minima and maxima
- **Molecular Structure Analysis**: Counts specific molecular species (e.g., HCN molecules) in the final structures to characterize reaction products
- **Low-Barrier Pathway Identification**: Filters and identifies reaction pathways with activation barriers below a specified threshold
- **Visualization**: Generates publication-quality plots showing energy profiles, peak/valley locations, and barrier heights
- **Data Export**: Outputs structured CSV files containing energy data and low-barrier pathway information

The tool is particularly useful for:
- Screening large numbers of reaction trajectories
- Identifying kinetically favorable reaction pathways
- Analyzing molecular dynamics or Monte Carlo simulation outputs
- Studying complex reaction networks with multiple pathways

## How to Use

### Prerequisites

Ensure you have the following Python packages installed:
```bash
pip install numpy matplotlib peakutils
```

### Directory Structure

Before running RPRA, set up the following directory structure:

```
project_directory/
├── peakDetect.py
├── calculations/         # Input XYZ files
├── energies/             # Output energy CSV files (created automatically)
├── plots/                # Output plots (created automatically)
└── utils/                # Utility modules
    ├── files.py
    ├── lists.py
    ├── peaks.py
    ├── pnds.py
    └── xyz.py
```

### Input Files

Place your XYZ trajectory files in the `./calculations/` directory. These files should contain:
- Multiple XYZ structures representing a reaction pathway
- Energy information in the format: `Energy: <value>` within the file

### Configuration

Edit the configuration variables in `peakDetect.py` as needed:

```python
XYZ_FOLDER = "./calculations/"    # Directory containing input XYZ files
ENERGY_DUMP = "./energies/"       # Directory for energy CSV output
PLOTS_DUMP = "./plots/"           # Directory for plot output
THRES = 0.005                     # Peak detection threshold
MIN_DIST = 1                      # Minimum distance between peaks
BARRIER_CUTOFF = 0.01             # Threshold for low-barrier identification
```

### Running the Analysis

Execute the main script:

```bash
python peakDetect.py
```

### Output Files

The tool generates several output files:

1. **Energy CSV files** (`./energies/<filename>.csv`): Individual energy profiles for each trajectory
2. **Peak plots** (`./plots/<filename>.png`): Visualization of energy profiles with identified peaks and valleys
3. **Barrier plots** (`./plots/<filename>_barriers.png`): Bar charts showing calculated activation barriers
4. **Low barrier summary** (`low_barrier_paths.csv`): Summary of all trajectories with barriers below the cutoff threshold

### Output Format

The `low_barrier_paths.csv` file contains:
- `filename`: Name of the trajectory file
- `n_hcn`: Number of HCN molecules identified in the final structure
- `steps/energies`: Complete energy profile as comma-separated values

### Customization

- **Peak Detection**: Adjust `THRES` and `MIN_DIST` parameters to fine-tune peak detection sensitivity
- **Barrier Cutoff**: Modify `BARRIER_CUTOFF` to change the threshold for low-barrier pathway identification
- **Molecular Analysis**: Edit the `find_hcn_molecules` function in `utils/xyz.py` to analyze different molecular species
- **Energy Extraction**: Modify the `REGEX` pattern to match different energy output formats

### References
- Lucas Hermann Negri, & Christophe Vestri. (2017). lucashn/peakutils: v1.1.0 (v1.1.0). Zenodo. <https://doi.org/10.5281/zenodo.887917>