import peakutils
import numpy as np
import matplotlib.pyplot as plt


def detect_peak(data, inverse=False, thres=0.005, min_dist=2):
    y = np.array(data, dtype='float64')
    ly = len(y)
    x = np.linspace(0, ly-1, ly, dtype='int64')

    if inverse:
        x_peaks = peakutils.indexes(-y, thres=thres, min_dist=min_dist)
    else:
        x_peaks = peakutils.indexes(y, thres=thres, min_dist=min_dist)
    
    return x, y, x_peaks


def plot_peaks(x,y, x_peaks, x_valleys, filename='figure'):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(x, y, "--")
    
    mask_peaks = np.ones(y.shape, bool)
    mask_peaks[x_peaks] = False
    y1 = np.ma.MaskedArray(y, mask_peaks)
    
    mask_valleys = np.ones(y.shape, bool)
    mask_valleys[x_valleys] = False
    y2 = np.ma.MaskedArray(y, mask_valleys)
    
    ax.plot(x, y1, label="Peaks", color='r', marker='o', fillstyle='none')
    ax.plot(x, y2, label="Valleys", color='g', marker='o', fillstyle='none')
    
    plt.title(filename)
    plt.legend()
    fig.tight_layout()
    # fig.savefig(f"./{filename}.png", format="png", dpi=300)
    # plt.show()
    # plt.close(fig)
    
    return fig


def cal_barriers(x,y, x_peaks, x_valleys):
    mask_peaks = np.ones(y.shape, bool)
    mask_peaks[x_peaks] = False
    peaks = np.ma.MaskedArray(y, mask_peaks)
    
    mask_valleys = np.ones(y.shape, bool)
    mask_valleys[x_valleys] = False
    valleys = np.ma.MaskedArray(y, mask_valleys)
    
    barriers = []
    for idx, p in enumerate(peaks):
        print(idx,p)
        
        if np.ma.is_masked(p):
            barriers.append(0)
            print('pass')
        else:
            try:
                v = valleys[:idx]
                b = p - v[v.mask == False][-1]
            except:
                b = 0
            
            print(f'>> {v}, {b}')
            barriers.append(b)
        
    return barriers


def plot_barriers(x, y, filename='barriers'):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(x,y)
    
    plt.title(filename)
    fig.tight_layout()
    # fig.savefig(f"./{filename}.png", format="png", dpi=300)
    # plt.show()
    # plt.close(fig)
    
    return fig


def collect_low_barrier_structures(barriers, energies, filename, barrier_cutoff=0.01, mol_count=0):
    large_barrier = any(value > barrier_cutoff for value in barriers)
    if not large_barrier:
        array_str = ','.join(map(str, energies))
        line = f'{filename},{mol_count},{array_str}\n'
        return line
    
    
    
    