from ROOT import TFile
import numpy as np
from sys import argv
from os import listdir
from os.path import isfile, join
import os, json
import argparse

from utils import read_config, get_data, get_patch

def main(argv):

    parser = argparse.ArgumentParser(description='Makes training data set for EM vs track separation')
    parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
    args = parser.parse_args()

    config = read_config(args.config)

    print '#'*50,'\nPrepare data for CNN'
    INPUT_TYPE   = config['prepare_data_em_track']['input_type']
    INPUT_DIR    = config['prepare_data_em_track']['input_dir']
    OUTPUT_DIR   = config['prepare_data_em_track']['output_dir']
    PATCH_SIZE_W = config['prepare_data_em_track']['patch_size_w']
    PATCH_SIZE_D = config['prepare_data_em_track']['patch_size_d']
    print 'Using %s as input dir, and %s as output dir' % (INPUT_DIR, OUTPUT_DIR)
    if INPUT_TYPE == 'root': print 'Reading from ROOT file'
    else: print 'Reading from TEXT files'
    print '#'*50

    doing_nue = config['prepare_data_em_track']['doing_nue']                       # set to true for nu_e events (will skip more showers)
    selected_view_idx = config['prepare_data_em_track']['selected_view_idx']       # set the view id
    patch_fraction = config['prepare_data_em_track']['patch_fraction']             # percent of used patches
    empty_fraction = config['prepare_data_em_track']['empty_fraction']             # percent of "empty background" patches
    clean_track_fraction = config['prepare_data_em_track']['clean_track_fraction'] # percent of selected patches, where only a clean track is present
    muon_track_fraction = config['prepare_data_em_track']['muon_track_fraction']   # ***** new: preselect muos, they are many *****
    crop_event = config['prepare_data_em_track']['crop_event']                     # use true only if no crop on LArSoft level and not a noise dump

    blur_kernel = np.asarray(config['prepare_data_em_track']['blur'])              # add blur in wire direction with given kernel if it is not empty (only for tests)
    white_noise = config['prepare_data_em_track']['noise']                         # add gauss noise with given sigma if value > 0 (only for tests)
    coherent_noise = config['prepare_data_em_track']['coherent']                   # add coherent (groups of 32 wires) gauss noise with given sigma if value > 0 (only for tests)

    print 'Using', patch_fraction, '% of data from view', selected_view_idx
    print 'Using', muon_track_fraction, '% of muon points'
    if doing_nue: print 'Neutrino mode, will skip more showers.'

    print 'Blur kernel', blur_kernel, 'noise RMS', white_noise

    max_capacity = 1700000
    db = np.zeros((max_capacity, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
    db_y = np.zeros((max_capacity, 4), dtype=np.int32)

    patch_area = PATCH_SIZE_W * PATCH_SIZE_D

    cnt_ind = 0
    cnt_trk = 0
    cnt_sh = 0
    cnt_michel = 0
    cnt_void = 0

    fcount = 0

    rootFile = None
    rootModule = 'datadump'
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
            raw, deposit, pdg, tracks, showers = get_data(folder, evname, PATCH_SIZE_D/2 + 2, crop_event, blur_kernel, white_noise, coherent_noise)
            if raw is None:
                print 'Skip empty event...'
                continue

            pdg_michel = ((pdg & 0xF000) == 0x2000)
            vtx_map = (pdg >> 24) & 0x3 # only 1 or 2: interaction or pi0 decay

            print 'Tracks', np.sum(tracks), 'showers', np.sum(showers), 'michels', np.sum(pdg_michel)

            sel_trk = 0
            sel_sh = 0
            sel_muon = 0
            sel_mu_near_stop = 0
            sel_michel = 0
            sel_empty = 0
            sel_clean_trk = 0
            for i in range(raw.shape[0]):
                for j in range(raw.shape[1]):
                    is_raw_zero = (raw[i,j] < 0.01)
                    is_michel = (pdg[i,j] & 0xF000 == 0x2000) # has michel flag set, wont skip it
                    is_muon = (pdg[i,j] & 0xFFF == 13)
                
                    is_vtx = (vtx_map[i,j] > 0)
                
                    x_start = np.max([0, i - PATCH_SIZE_W/2])
                    x_stop  = np.min([raw.shape[0], x_start + PATCH_SIZE_W])

                    y_start = np.max([0, j - PATCH_SIZE_D/2])
                    y_stop  = np.min([raw.shape[1], y_start + PATCH_SIZE_D])

                    if x_stop - x_start != PATCH_SIZE_W or y_stop - y_start != PATCH_SIZE_D:
                        continue

                    is_mu_near_stop = False
                    if is_muon:
                        pdg_patch = pdg_michel[x_start+2:x_stop-2, y_start+2:y_stop-2]
                        if np.count_nonzero(pdg_patch) > 2:
                            is_mu_near_stop = True
                            sel_mu_near_stop += 1
                
                    eff_patch_fraction = patch_fraction
                    near_vtx_count = np.count_nonzero(vtx_map[x_start+2:x_stop-2, y_start+2:y_stop-2])
                    if near_vtx_count > 0 and eff_patch_fraction < 50:
                        eff_patch_fraction = 50 # min half of points near vertices

                    # randomly skip fraction of patches
                    if not(is_michel | is_mu_near_stop | is_vtx) and np.random.randint(10000) > int(100*eff_patch_fraction): continue

                    track_pixels = np.count_nonzero(tracks[x_start:x_stop, y_start:y_stop])
                    shower_pixels = np.count_nonzero(showers[x_start:x_stop, y_start:y_stop])

                    target = np.zeros(4, dtype=np.int32)
                    if tracks[i,j] == 1:
                        target[0] = 1 # label as a track-like
                        if is_raw_zero: continue
                        # skip fraction of almost-track-only patches
                        if shower_pixels < 8 and near_vtx_count == 0:
                            if np.random.randint(10000) > int(100*clean_track_fraction): continue
                            else: sel_clean_trk += 1
                        else:
                            if is_muon:
                                if not(is_mu_near_stop) and np.random.randint(10000) > int(100*muon_track_fraction): continue
                                sel_muon += 1
                        cnt_trk += 1
                        sel_trk += 1
                    elif showers[i,j] == 1:
                        target[1] = 1 # label as a em-like
                        if doing_nue: # for nu_e events (lots of showers) skip some fraction of shower patches
                            if near_vtx_count == 0 and np.random.randint(100) < 40: continue  # skip 40% of any shower
                            if shower_pixels > 0.05*patch_area and np.random.randint(100) < 50: continue
                            if shower_pixels > 0.20*patch_area and np.random.randint(100) < 90: continue
                        if is_raw_zero: continue
                        if is_michel:
                            target[2] = 1 # additionally label as a michel electron
                            cnt_michel += 1
                            sel_michel += 1
                        cnt_sh += 1
                        sel_sh += 1
                    else: # use small fraction of empty-but-close-to-something patches
                        target[3] = 1 # label an empty pixel
                        if np.random.randint(10000) < int(100*empty_fraction):
                            nclose = np.count_nonzero(showers[i-2:i+3, j-2:j+3])
                            nclose += np.count_nonzero(tracks[i-2:i+3, j-2:j+3])
                            if nclose == 0:
                                npix = shower_pixels + track_pixels
                                if npix > 6:
                                    cnt_void += 1
                                    sel_empty += 1
                                else: continue # completely empty patch
                            else: continue # too close to trk/shower
                        else: continue # not selected randomly

                    if np.count_nonzero(target) == 0:
                        print 'NEED A LABEL IN THE TARGET!!!'
                        continue

                    if cnt_ind < max_capacity:
                        db[cnt_ind] = get_patch(raw, i, j, PATCH_SIZE_W, PATCH_SIZE_D)
                        db_y[cnt_ind] = target
                        cnt_ind += 1
                    else:
                        print 'MAX CAPACITY REACHED!!!'
                        break

            print 'Selected: tracks', sel_trk, 'showers', sel_sh, 'empty', sel_empty, '/// muons', sel_muon, 'michel', sel_michel, 'clean trk', sel_clean_trk

    print 'Added', cnt_ind, 'tracks:', cnt_trk, 'showers:', cnt_sh, 'michels:', cnt_michel, 'empty:', cnt_void

    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_x', db[:cnt_ind])
    np.save(OUTPUT_DIR+'/db_view_'+str(selected_view_idx)+'_y', db_y[:cnt_ind])

if __name__ == "__main__":
    main(argv)

