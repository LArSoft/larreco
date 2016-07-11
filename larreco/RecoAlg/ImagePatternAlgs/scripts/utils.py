import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os, json

def get_event_bounds(A):
    # get center with 99% of signal
    cum = np.cumsum(np.sum(A, axis=0))
    start_ind = np.where(cum > cum[-1]*0.005)[0][0]
    end_ind = np.where(cum > cum[-1]*0.995)[0][0]
    return start_ind, end_ind

def get_data(fname):
    print 'Reading', fname + '.raw'
    A_raw     = np.genfromtxt(fname + '.raw', delimiter=' ')
    A_deposit = np.genfromtxt(fname + '.deposit', delimiter=' ')
    A_pdg     = np.genfromtxt(fname + '.pdg', delimiter=' ')

    print np.sum(A_raw), np.sum(A_deposit), np.sum(A_pdg)
    #assert np.sum(A_deposit) > 0
    # get main event body (99% signal)
    evt_start_ind, evt_stop_ind = get_event_bounds(A_raw)
    print evt_start_ind, evt_stop_ind
    A_raw     = A_raw[:,evt_start_ind:evt_stop_ind]
    A_deposit = A_deposit[:,evt_start_ind:evt_stop_ind]
    A_pdg     = np.abs(A_pdg[:,evt_start_ind:evt_stop_ind]) # ABS

    # ugly fix
    deposit_th_ind = A_deposit < 2.0e-5
    A_pdg[deposit_th_ind] = 0
    tracks = A_pdg.copy()
    showers = A_pdg.copy()
    tracks[A_pdg == 11] = 0
    tracks[tracks > 0]   = 1

    showers[A_pdg != 11] = 0
    showers[showers > 0] = 1

    return A_raw, A_deposit, A_pdg, tracks, showers


def read_config():
    config = None
    with open('config.json', 'r') as fin:
        config = json.loads(fin.read());
    if config is None:
        print 'This script requires configuration file: config.json'
        exit(1)
    INPUT_DIR = config['prepare_data']['input_dir']
    OUTPUT_DIR = config['prepare_data']['output_dir']
    PATCH_SIZE = config['prepare_data']['patch_size']
    return INPUT_DIR, OUTPUT_DIR, PATCH_SIZE

def main(argv):

    print '#'*50,'\nPrepare data for CNN'
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

        db = np.zeros((total_learning_patches, PATCH_SIZE, PATCH_SIZE))
        db_y = np.zeros((total_learning_patches, 1))

        cnt_ind = 0
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
                #print '#',i,j, target, raw.shape
                #print x_start, x_stop
                #print y_start, y_stop
                db[cnt_ind] = raw[x_start:x_stop, y_start:y_stop]
                db_y[cnt_ind] = target
                cnt_ind += 1

        print 'Added', cnt_ind, 'of total expected', total_learning_patches

        np.save(OUTPUT_DIR+'/db_evt_' + evt_no, db[:cnt_ind])
        np.save(OUTPUT_DIR+'/db_y_evt_' + evt_no, db_y[:cnt_ind])



if __name__ == "__main__":
    main(sys.argv)
