import utils.files as fl
import utils.lists as ls
import utils.peaks as pk
import utils.xyz   as mol


# configs
XYZ_FOLDER  = "./calculations/"
ENERGY_DUMP = "./energies/"
PLOTS_DUMP  = "./plots/"
CWD         = fl.get_cwd()

THRES=0.005
MIN_DIST=1

REGEX = r"\bEnergy:\s+(.*)$"

BARRIER_CUTOFF = 0.01


files = fl.search_files(type="xyz", location=XYZ_FOLDER)
low_barriers_txt = 'filename,n_hcn,steps/energies\n'

for xyz in files:
    
    xyz_data = mol.read_xyz_file(xyz)
    first_xyz_data, last_xyz_data = mol.extract_xyz_files(xyz_data)
    
    # o_count, h_count = mol.find_closest_h_atoms(last_xyz_data)
    # h2o_count=f'{h_count/o_count} Hx{h_count} Ox{o_count}'
    
    hcn_count = mol.find_hcn_molecules(last_xyz_data, h_c_cutoff=1.25, c_n_cutoff=1.22)
    
    energy = fl.grep(filename=xyz, pattern=REGEX)
    ls.to_csv(data=energy, filename=CWD/ENERGY_DUMP/f"{xyz.name}.csv", orientation='col')
    
    # detect peaks
    x,y,x_peaks   = pk.detect_peak(energy, inverse=False, thres=THRES, min_dist=MIN_DIST)
    _,_,x_valleys = pk.detect_peak(energy, inverse=True,  thres=THRES, min_dist=MIN_DIST)
    y_peaks       = y[x_peaks]
    y_valleys     = y[x_valleys]
    
    print(xyz)
    print(y_peaks)
    print(y_valleys)
    
    barriers = pk.cal_barriers(x,y, x_peaks, x_valleys)
    print(barriers)
    
    # Identify large barriers. TRUE if barrier > BARRIER_CUTOFF
    l = pk.collect_low_barrier_structures(barriers, energies=y, filename=xyz.name, barrier_cutoff=BARRIER_CUTOFF, mol_count=hcn_count)
    if l is not None:
        low_barriers_txt=low_barriers_txt+l
    
    print('---')
    
    # plotting
    fig = pk.plot_peaks(x,y, x_peaks, x_valleys, xyz.name)
    fig.savefig(CWD/PLOTS_DUMP/f"{xyz.name}.png", format="png", dpi=300)
    
    fig = pk.plot_barriers(x, barriers, xyz.name)
    fig.savefig(CWD/PLOTS_DUMP/f"{xyz.name}_barriers.png", format="png", dpi=300)

# Dump data
with open('low_barrier_paths.csv', 'w') as f:
    f.write(low_barriers_txt)
    
del low_barriers_txt
    

