from ROOT import TFile
import numpy as np
from sys import argv
from os import listdir
from os.path import isfile, join
import os, json
import argparse

from utils import read_config, get_data, get_patch, get_vertices, get_nu_vertices

def main(argv):

    parser = argparse.ArgumentParser(description='Makes training data set for various vertex/decay ID')
    parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
    parser.add_argument('-t', '--type', help="Input file format")
    parser.add_argument('-i', '--input', help="Input directory")
    parser.add_argument('-o', '--output', help="Output directory")
    parser.add_argument('-v', '--view', help="view")
    args = parser.parse_args()

    config = read_config(args.config)

    print '#'*50,'\nPrepare data for CNN'
    INPUT_TYPE   = config['prepare_data_vtx_id']['input_type']
    INPUT_DIR = config['prepare_data_vtx_id']['input_dir']
    OUTPUT_DIR = config['prepare_data_vtx_id']['output_dir']
    PATCH_SIZE_W = config['prepare_data_vtx_id']['patch_size_w']
    PATCH_SIZE_D = config['prepare_data_vtx_id']['patch_size_d']
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    print '#'*50

    rootModule = config['prepare_data_vtx_id']['module_name']                # larsoft module name used for data dumps in ROOT format
    selected_view_idx = config['prepare_data_vtx_id']['selected_view_idx']   # set the view id
    nearby_empty = config['prepare_data_vtx_id']['nearby_empty']             # number of patches near each vtx, but with empty area in the central pixel
    nearby_on_track = config['prepare_data_vtx_id']['nearby_on_track']       # number of patches on tracks or showers, somewhere close to each vtx
    crop_event = config['prepare_data_vtx_id']['crop_event']                 # use true only if no crop on LArSoft level and not a noise dump

    print 'Using', nearby_empty, 'empty and', nearby_on_track, 'on track patches per each verex in view', selected_view_idx

    max_capacity = 300000
    db = np.zeros((max_capacity, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
    db_y = np.zeros((max_capacity, 5), dtype=np.int32)

    kHadr  = 0x1   # hadronic inelastic scattering
    kPi0   = 0x2   # pi0 produced in this vertex
    kDecay = 0x4   # point of particle decay (except pi0 decays)
    kConv  = 0x8   # gamma conversion

    cnt_ind = 0
    cnt_vtx = 0
    cnt_decay = 0
    cnt_gamma = 0
    cnt_nu = 0
    cnt_trk = 0
    cnt_void = 0

    fcount = 0

    event_list = []
    if INPUT_TYPE == "root":
        fnames = [f for f in os.listdir(INPUT_DIR) if '.root' in f]
        for n in fnames:
            rootFile = TFile(INPUT_DIR+'/'+n)
            keys = [rootModule+'/'+k.GetName()[:-4] for k in rootFile.Get(rootModule).GetListOfKeys() if '_raw' in k.GetName()]
            event_list.append((rootFile, keys))
    else:
        keys = [f[:-4] for f in os.listdir(INPUT_DIR) if '.raw' in f] # only main part of file name, without extension
        event_list.append((INPUT_DIR, keys)) # single entry in the list of txt files

    for entry in event_list:
        folder = entry[0]
        event_names = entry[1]

        for evname in event_names:
            finfo = evname.split('_')
            evt_no = finfo[2]
            tpc_idx = int(finfo[8])
            view_idx = int(finfo[10])

            if view_idx != selected_view_idx: continue
            fcount += 1

            print 'Process event', fcount, evname, 'NO.', evt_no

            # get clipped data, margin depends on patch size in drift direction
            raw, deposit, pdg, tracks, showers = get_data(folder, evname, PATCH_SIZE_D/2 + 2, crop_event)
            if raw is None:
                print 'Skip empty event...'
                continue

            vtx = get_vertices(pdg)
            nuvtx = get_nu_vertices(pdg)
            print 'Found', vtx.shape[0], 'hadronic vertices/decay', nuvtx.shape[0], 'neutrino vertices'

            for v in range(vtx.shape[0]):
                flags = 0
                if vtx.shape[0] > 0:
                    flags = vtx[v,2]

                nuflags = 0
                if nuvtx.shape[0] > 0:
                    nuflags = nuvtx[v,2]

                if (flags & kHadr) > 0 or (flags & kDecay) > 0 or ((flags & kPi0) > 0 and (flags & kConv) > 0):

                    wire = vtx[v,0]
                    drif = vtx[v,1]

                    x_start = np.max([0, wire - PATCH_SIZE_W/2])
                    x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE_W])

                    y_start = np.max([0, drif - PATCH_SIZE_D/2])
                    y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE_D])

                    if x_stop - x_start != PATCH_SIZE_W or y_stop - y_start != PATCH_SIZE_D:
                        continue

                    target = np.zeros(5, dtype=np.int32)     # [decay, hadronic_vtx, g_conversion, nu_primary, not_vtx]
                    if nuflags > 0:
                        target[3] = 1
                        cnt_nu += 1
                    elif (flags & kDecay) > 0:
                        target[0] = 1
                        cnt_decay += 1
                    elif (flags & kHadr) > 0:
                        target[1] = 1
                        cnt_vtx += 1
                    elif (flags & kConv) > 0:
                        target[2] = 1
                        cnt_gamma += 1

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
                                    target = np.zeros(5, dtype=np.int32)
                                    target[4] = 1
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
                                    target = np.zeros(5, dtype=np.int32)
                                    target[4] = 1
                                    db[cnt_ind] = patch
                                    db_y[cnt_ind] = target
                                    cnt_trk += 1
                                    cnt_ind += 1
                                    n_track += 1
                                else: break
                        n_trials += 1

    print 'Total size', cnt_ind, ':: hadronic:', cnt_vtx, 'decays:', cnt_decay, 'nu-primary:', cnt_nu, 'g-conv:', cnt_gamma, 'empty:', cnt_void, 'on-track:', cnt_trk

    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_x', db[:cnt_ind])
    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_y', db_y[:cnt_ind])

if __name__ == "__main__":
    main(argv)

