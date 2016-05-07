import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os, json

from utils import read_config, get_data

def write_event(f, evt_no, i, j, data_row, label):
    f.write(evt_no + ' ' + str(i) + ' ' + str(j) + ' ')
    for k in range(data_row.shape[0]): f.write(str(data_row[k]) + ' ')
    for k in range(label.shape[0]): f.write(str(label[k]) + ' ')
    f.write('\n')
    return 0

def main(argv):

    print '#'*50,'\nPrepare data for CNN'
    INPUT_DIR, OUTPUT_DIR, PATCH_SIZE = read_config()
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    flat_txt_training = open(OUTPUT_DIR+'/flat_training_3class.prn', 'w')
    flat_txt_testing = open(OUTPUT_DIR+'/flat_testing_3class.prn', 'w')
    
    TEST_EVENTS_SPLIT = 1000
    
    cnt_ind = 0
    cnt_trk = 0
    cnt_sh = 0
    cnt_void = 0
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

        for i in range(raw.shape[0]):
            for j in range(raw.shape[1]):
                x_start = np.max([0, i - PATCH_SIZE/2])
                x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE])

                y_start = np.max([0, j - PATCH_SIZE/2])
                y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE])

                if x_stop - x_start != PATCH_SIZE or y_stop - y_start != PATCH_SIZE:
                    continue
                    
                target = np.zeros(3)
                if tracks[i,j] == 1:
                    # skip some fraction (1/2) of almost-track-only patches
                    if np.random.randint(2) < 1 and np.count_nonzero(showers[x_start:x_stop, y_start:y_stop]) < 4: continue
                    else:
                        target[0] = 1
                        cnt_trk += 1
                elif showers[i,j] == 1:
                    target[1] = 1
                    cnt_sh += 1
                else:
                    # use small fraction (1/130) of almost-empty patches as 3rd class
                    if np.random.randint(130) < 1:
                        npix = np.count_nonzero(showers[x_start:x_stop, y_start:y_stop])
                        npix += np.count_nonzero(tracks[x_start:x_stop, y_start:y_stop])
                        if npix > 6:
                            target[2] = 1
                            cnt_void += 1
                        else: continue
                    else: continue

                cnt_ind += 1
                        
                training_patch = raw[x_start:x_stop, y_start:y_stop]
                flat_patch = training_patch.flatten()
                
                if f_idx < TEST_EVENTS_SPLIT:
                    write_event(flat_txt_training, evt_no, i, j, flat_patch, target)
                else:
                    write_event(flat_txt_testing, evt_no, i, j, flat_patch, target)


    print '*** Total samples:', cnt_ind
    print '*** Track samples:', cnt_trk, 'shower samples:', cnt_sh, 'empty samples:', cnt_void

    flat_txt_training.close()
    flat_txt_testing.close()

if __name__ == "__main__":
    main(sys.argv)
