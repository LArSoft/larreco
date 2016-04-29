import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os, json
from scipy import ndimage
from skimage.measure import block_reduce

from utils import read_config, get_data

def main(argv):

    print '#'*50,'\nPrepare data for CNN'
    INPUT_DIR, OUTPUT_DIR, PATCH_SIZE = read_config()
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    files = [f for f in os.listdir(INPUT_DIR) if '.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        evt_no = fname.split('_')[2]
        print 'Process', fname, 'EVT', evt_no

        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname)
        total_learning_patches = int(np.sum(tracks) + np.sum(showers))
        print 'Tracks', np.sum(tracks), 'showers', np.sum(showers)
        print 'Total learning patches', total_learning_patches

        db = np.zeros((total_learning_patches, PATCH_SIZE, PATCH_SIZE))
        db_y = np.zeros((total_learning_patches, 1))

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

                #print '#',i,j, target, raw.shape
                #print x_start, x_stop
                #print y_start, y_stop
                db[cnt_ind] = raw[x_start:x_stop, y_start:y_stop]
                db_y[cnt_ind] = target
                
                if target == 1: cnt_trk += 1
                else: cnt_sh += 1
                cnt_ind += 1

        print 'Added', cnt_ind, 'of total expected', total_learning_patches
        print 'Tracks:', cnt_trk, 'showers:', cnt_sh

        np.save(OUTPUT_DIR+'/db_evt_' + evt_no, db[:cnt_ind])
        np.save(OUTPUT_DIR+'/db_y_evt_' + evt_no, db_y[:cnt_ind])



if __name__ == "__main__":
    main(sys.argv)
