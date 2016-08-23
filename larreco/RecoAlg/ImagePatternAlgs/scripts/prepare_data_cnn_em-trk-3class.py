import numpy as np
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os, json

from utils import read_config, get_data

def get_patch(a, wire, drift, wsize, dsize):
    halfSizeW = wsize / 2;
    halfSizeD = dsize / 2;

    w0 = wire - halfSizeW;
    w1 = wire + halfSizeW;

    d0 = drift - halfSizeD;
    d1 = drift + halfSizeD;

    patch = np.zeros((wsize, dsize))

    wpatch = 0
    for w in range(w0, w1):
        if w >= 0 and w < a.shape[0]:
            dpatch = 0
            for d in range(d0, d1):
                if d >= 0 and d < a.shape[1]:
                    patch[wpatch,dpatch] = a[w,d];
                dpatch += 1
        wpatch += 1
    
    return patch


def main(argv):

    print '#'*50,'\nPrepare data for CNN'
    INPUT_DIR, OUTPUT_DIR, PATCH_SIZE_W, PATCH_SIZE_D = read_config()
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    db = np.zeros((1800000, PATCH_SIZE_W, PATCH_SIZE_D))
    db_y = np.zeros((1800000, 3), dtype=np.int32)

    cnt_ind = 0
    cnt_trk = 0
    cnt_sh = 0
    cnt_void = 0

    files = [f for f in os.listdir(INPUT_DIR) if '.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        evt_no = fname.split('_')[2]
        print 'Process', fname, 'EVT', evt_no

        # get clipped data, margin depends on patch size in drift direction
        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname, PATCH_SIZE_D/2 + 2)
        print 'Tracks', np.sum(tracks), 'showers', np.sum(showers)

        for i in range(raw.shape[0]):
            for j in range(raw.shape[1]):
                # randomly skip (3/4) of patches
                if np.random.randint(100) < 25: continue
                
                x_start = np.max([0, i - PATCH_SIZE_W/2])
                x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE_W])

                y_start = np.max([0, j - PATCH_SIZE_D/2])
                y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE_D])

                if x_stop - x_start != PATCH_SIZE_W or y_stop - y_start != PATCH_SIZE_D:
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
                    # use small fraction (1/100) of empty-but-close-to-something patches
                    if np.random.randint(100) < 1:
                        nclose = np.count_nonzero(showers[i-2:i+3, j-2:j+3])
                        nclose += np.count_nonzero(tracks[i-2:i+3, j-2:j+3])
                        if nclose == 0:
                            npix = np.count_nonzero(showers[x_start:x_stop, y_start:y_stop])
                            npix += np.count_nonzero(tracks[x_start:x_stop, y_start:y_stop])
                            if npix > 6:
                                target[2] = 1
                                cnt_void += 1
                            else: continue # completely empty patch
                        else: continue # too close to trk/shower
                    else: continue # not selected randomly

                db[cnt_ind] = get_patch(raw, i, j, PATCH_SIZE_W, PATCH_SIZE_D)
                db_y[cnt_ind] = target
                
                cnt_ind += 1

    print 'Added', cnt_ind, 'Tracks:', cnt_trk, 'showers:', cnt_sh, 'empty:', cnt_void

    np.save(OUTPUT_DIR+'/db_x', db[:cnt_ind])
    np.save(OUTPUT_DIR+'/db_y', db_y[:cnt_ind])

if __name__ == "__main__":
    main(sys.argv)

