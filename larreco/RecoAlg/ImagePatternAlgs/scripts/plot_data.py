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

from utils import read_config, get_data, get_vertices

def main(argv):

    config = read_config('config.json')
    INPUT_DIR = config['prepare_data_em_track']['input_dir']
    PATCH_SIZE_D = config['prepare_data_em_track']['patch_size_d']
    crop_event = config['prepare_data_em_track']['crop_event']

    files = [f for f in os.listdir(INPUT_DIR) if 'raw_event_1_run_1_subrun_0_tpc_2_view_2.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        evt_no = fname.split('_')[2]
        print 'Process', fname, 'EVT', evt_no

        # get clipped data, margin depends on patch size in drift direction
        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname, PATCH_SIZE_D/2 + 2, crop_event)

        vtx = get_vertices(pdg)
        print 'found', vtx.shape[0], 'vertices'
        print vtx[:,0]
        print vtx[:,1]
        print vtx[:,2]
        
        fig, ax = subplots(2,2, figsize=(15, 15))
        ax[0,0].imshow(-raw.transpose(), cmap='gray', interpolation='none', aspect='auto')
        ax[0,0].scatter(vtx[:,0]+0.5, vtx[:,1]+0.5, s=50, c=32*vtx[:,2], cmap='rainbow')
        ax[0,0].set_title('ADC')
        ax[0,1].imshow(-deposit.transpose(), cmap='gray', interpolation='none', aspect='auto')
        ax[0,1].set_title('DEPOSIT')
        ax[1,0].imshow(-tracks.transpose(), cmap='gray', interpolation='none', aspect='auto')
        ax[1,0].set_title('TRACKS LABEL')
        ax[1,1].imshow(-showers.transpose(), cmap='gray', interpolation='none', aspect='auto')
        ax[1,1].set_title('SHOWERS LABEL')
        plt.show()
        #break


if __name__ == "__main__":
    main(sys.argv)

