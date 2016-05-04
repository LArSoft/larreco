import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os, json

from utils import read_config, get_data

def main(argv):

    print '#'*50,'\nPrepare data for CNN'
    INPUT_DIR, OUTPUT_DIR, PATCH_SIZE = read_config()
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50
    
    flat_txt_training = open(OUTPUT_DIR+'/flat_training.prn', 'w')
    flat_txt_testing = open(OUTPUT_DIR+'/flat_testing.prn', 'w')

    f_idx = 0
    files = [f for f in os.listdir(INPUT_DIR) if '.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        evt_no = fname.split('_')[2]
        ievt = int(evt_no)

        f_idx += 1
        
        print 'Process', fname, 'EVT', evt_no

        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname)
        total_learning_patches = int(np.sum(tracks) + np.sum(showers))
        print 'Tracks', np.sum(tracks), 'showers', np.sum(showers)
        print 'Total learning patches', total_learning_patches

        cnt_ind = 0
        cnt_trk = 0
        cnt_sh = 0
        for i in range(raw.shape[0]):
            for j in range(raw.shape[1]):
                target = -1
                if tracks[i,j] == 1:
                    target = 1
                elif showers[i,j] == 1:
                    target = 0
                if target == -1:
                    continue

                patch = np.zeros((PATCH_SIZE, PATCH_SIZE))
                x_start = np.max([0, i - PATCH_SIZE/2])
                x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE])

                y_start = np.max([0, j - PATCH_SIZE/2])
                y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE])

                if x_stop - x_start != PATCH_SIZE or y_stop - y_start != PATCH_SIZE:
                    continue

                # skip some fraction (1/2) of almost-track-only patches
                if target == 1 and np.count_nonzero(showers[x_start:x_stop, y_start:y_stop]) < 4 and np.random.randint(2) < 1:
                    continue
                        
                training_patch = raw[x_start:x_stop, y_start:y_stop]
                flat_patch = training_patch.flatten()
                
                if f_idx < 70:
                    flat_txt_training.write(evt_no + ' ' + str(i) + ' ' + str(j) + ' ')
                    for k in range(flat_patch.shape[0]):
                        flat_txt_training.write(str(flat_patch[k]) + ' ')
                    flat_txt_training.write(str(target) + '\n')
                else:
                    flat_txt_testing.write(evt_no + ' ' + str(i) + ' ' + str(j) + ' ')
                    for k in range(flat_patch.shape[0]):
                        flat_txt_testing.write(str(flat_patch[k]) + ' ')
                    flat_txt_testing.write(str(target) + '\n')

                if target == 1: cnt_trk += 1
                else: cnt_sh += 1

                cnt_ind += 1

        print 'Added', cnt_ind, 'of total expected', total_learning_patches
        print 'Track patches:', cnt_trk, 'shower patches:', cnt_sh

    flat_txt_training.close()
    flat_txt_testing.close()

if __name__ == "__main__":
    main(sys.argv)
