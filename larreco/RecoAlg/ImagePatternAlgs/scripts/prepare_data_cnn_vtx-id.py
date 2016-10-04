import numpy as np
from sys import argv
from os import listdir
from os.path import isfile, join
import os, json
import argparse

from utils import read_config, get_data, get_patch, get_vertices

def main(argv):

    parser = argparse.ArgumentParser(description='Makes training data set for various vertex/decay ID')
    parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
    args = parser.parse_args()

    config = read_config(args.config)

    print '#'*50,'\nPrepare data for CNN'
    INPUT_DIR = config['prepare_data_vtx_id']['input_dir']
    OUTPUT_DIR = config['prepare_data_vtx_id']['output_dir']
    PATCH_SIZE_W = config['prepare_data_vtx_id']['patch_size_w']
    PATCH_SIZE_D = config['prepare_data_vtx_id']['patch_size_d']
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    selected_view_idx = config['prepare_data_vtx_id']['selected_view_idx']   # set the view id
    nearby_empty = config['prepare_data_vtx_id']['nearby_empty']             # number of patches near each vtx, but with empty area in the central pixel
    nearby_on_track = config['prepare_data_vtx_id']['nearby_on_track']       # number of patches on tracks, somewhere close to each vtx
    crop_event = config['prepare_data_vtx_id']['crop_event']                 # use true only if no crop on LArSoft level and not a noise dump

    print 'Using', nearby_empty, 'empty and', nearby_on_track, 'on track patches per each verex in view', selected_view_idx

    max_capacity = 500000
    db = np.zeros((max_capacity, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
    db_y = np.zeros((max_capacity, 3), dtype=np.int32)

    kHadr  = 0x1   # hadronic inelastic scattering
    kPi0   = 0x2   # pi0 produced in this vertex
    kDecay = 0x4   # point of particle decay
    kConv  = 0x8   # gamma conversion

    cnt_ind = 0
    cnt_vtx = 0
    cnt_decay = 0
    cnt_trk = 0
    cnt_void = 0

    fcount = 0

    files = [f for f in os.listdir(INPUT_DIR) if '.raw' in f]
    for fname in files:
        fname = fname[:-4] # only main part of file, without extension
        finfo = fname.split('_')
        evt_no = finfo[2]
        tpc_idx = int(finfo[8])
        view_idx = int(finfo[10])

        if view_idx != selected_view_idx: continue
        fcount += 1

        print 'Process file', fcount, fname, 'EVT', evt_no

        # get clipped data, margin depends on patch size in drift direction
        raw, deposit, pdg, tracks, showers = get_data(INPUT_DIR+'/'+fname, PATCH_SIZE_D/2 + 2, crop_event)
        if raw is None:
            print 'Skip empty event...'
            continue

        vtx = get_vertices(pdg)
        print 'Found', vtx.shape[0], 'vertices'

        for v in range(vtx.shape[0]):
            flags = vtx[v,2]
            if (flags & kHadr) > 0 or (flags & kDecay) > 0 or ((flags & kPi0) > 0 and (flags & kConv) > 0):

                wire = vtx[v,0]
                drif = vtx[v,1]

                x_start = np.max([0, wire - PATCH_SIZE_W/2])
                x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE_W])

                y_start = np.max([0, drif - PATCH_SIZE_D/2])
                y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE_D])

                if x_stop - x_start != PATCH_SIZE_W or y_stop - y_start != PATCH_SIZE_D:
                    continue

                target = np.zeros(3)
                if (flags & kDecay) > 0:
                    target[0] = 1
                    cnt_vtx += 1
                else:
                    target[1] = 1
                    cnt_decay += 1

                patch = get_patch(raw, wire, drif, PATCH_SIZE_W, PATCH_SIZE_D)
                if cnt_ind < max_capacity:
                    db[cnt_ind] = patch
                    db_y[cnt_ind] = target
                    cnt_ind += 1
                else: break

                n_empty = 0
                n_trials = 0
                while n_empty < nearby_empty and n_trials < 500:
                    wi = np.random.randint(x_start+1, x_stop-1)
                    di = np.random.randint(y_start+1, y_stop-1)
                    if (wi < wire-1 or wi > wire+1) and (di < drif-2 or di > drif+2):
                        if tracks[wi,di] == 0 and showers[wi,di] == 0:
                            if cnt_ind < max_capacity:
                                patch = get_patch(raw, wi, di, PATCH_SIZE_W, PATCH_SIZE_D)
                                target = np.zeros(3)
                                target[2] = 1
                                db[cnt_ind] = patch
                                db_y[cnt_ind] = target
                                cnt_void += 1
                                cnt_ind += 1
                                n_empty += 1
                            else: break
                    n_trials += 1

                n_track = 0
                n_trials = 0
                while n_track < nearby_on_track and n_trials < 500:
                    wi = np.random.randint(x_start+1, x_stop-1)
                    di = np.random.randint(y_start+1, y_stop-1)
                    if (wi < wire-1 or wi > wire+1) and (di < drif-2 or di > drif+2):
                        if tracks[wi,di] == 1 or showers[wi,di] == 1:
                            if cnt_ind < max_capacity:
                                patch = get_patch(raw, wi, di, PATCH_SIZE_W, PATCH_SIZE_D)
                                target = np.zeros(3)
                                target[2] = 1
                                db[cnt_ind] = patch
                                db_y[cnt_ind] = target
                                cnt_trk += 1
                                cnt_ind += 1
                                n_track += 1
                            else: break
                    n_trials += 1

    print 'Total size', cnt_ind, ':: interactions:', cnt_vtx, 'decays:', cnt_decay, 'empty:', cnt_void, 'on-track', cnt_trk

    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_x', db[:cnt_ind])
    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_y', db_y[:cnt_ind])

if __name__ == "__main__":
    main(argv)

