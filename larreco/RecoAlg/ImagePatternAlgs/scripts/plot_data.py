import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os,json
from scipy import ndimage
from skimage.measure import block_reduce

from utils import read_config, get_data
from prepare_data import get_data


def main(argv):

    print '#'*50,'\nPlot data'
    INPUT_DIR, OUTPUT_DIR, PATCH_SIZE = read_config()
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    files = [f for f in os.listdir(INPUT_DIR) if 'tpc_2.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        evt_no = fname.split('_')[2]
        print 'Process', fname, 'EVT', evt_no

        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname)
        total_learning_patches = int(np.sum(tracks) + np.sum(showers))
        print 'Tracks', np.sum(tracks), 'showers', np.sum(showers)
        print 'Total learning patches', total_learning_patches

        fig, ax = subplots(2,2, figsize=(15, 15))
        ax[0,0].imshow(-raw, cmap='gray', interpolation='none', aspect='auto')
        ax[0,0].set_title('ADC')
        ax[0,1].imshow(-deposit, cmap='gray', interpolation='none', aspect='auto')
        ax[0,1].set_title('DEPOSIT')
        ax[1,0].imshow(-tracks, cmap='gray', interpolation='none', aspect='auto')
        ax[1,0].set_title('TRACKS LABEL')
        ax[1,1].imshow(-showers, cmap='gray', interpolation='none', aspect='auto')
        ax[1,1].set_title('SHOWERS LABEL')
        plt.show()
        #break


if __name__ == "__main__":
    main(sys.argv)
