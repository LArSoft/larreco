import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from os import listdir
from os.path import isfile, join
import os, json

def get_event_bounds(A, drift_margin = 0):
    # get center with 99% of signal
    cum = np.cumsum(np.sum(A, axis=0))
    start_ind = np.max([0, np.where(cum > cum[-1]*0.005)[0][0] - drift_margin])
    end_ind = np.min([A.shape[1], np.where(cum > cum[-1]*0.995)[0][0] + drift_margin])
    return start_ind, end_ind

def get_data(fname, drift_margin = 0):
    print 'Reading', fname + '.raw'
    A_raw     = np.genfromtxt(fname + '.raw', delimiter=' ')
    A_deposit = np.genfromtxt(fname + '.deposit', delimiter=' ')
    A_pdg     = np.genfromtxt(fname + '.pdg', delimiter=' ', dtype=int32)

    print np.sum(A_raw), np.sum(A_deposit), np.sum(A_pdg)
    #assert np.sum(A_deposit) > 0
    # get main event body (99% signal)
    evt_start_ind, evt_stop_ind = get_event_bounds(A_deposit, drift_margin)
    print evt_start_ind, evt_stop_ind
    A_raw     = A_raw[:,evt_start_ind:evt_stop_ind]
    A_deposit = A_deposit[:,evt_start_ind:evt_stop_ind]
    A_pdg     = A_pdg[:,evt_start_ind:evt_stop_ind]

    deposit_th_ind = A_deposit < 2.0e-5
    A_pdg[deposit_th_ind] = 0
    tracks = A_pdg.copy()
    showers = A_pdg.copy()
    tracks[(A_pdg & 0xFFFF) == 11] = 0
    tracks[tracks > 0]   = 1
    showers[(A_pdg & 0xFFFF) != 11] = 0
    showers[showers > 0] = 1
    return A_raw, A_deposit, A_pdg, tracks, showers

def get_vertices(A, mask):
    vtx = np.zeros((A.shape[0]*A.shape[1],3))
    nvtx = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if (A[i,j] & mask) > 0:
                t = A[i,j] >> 24
                v = np.zeros(3)
                v[0] = i
                v[1] = j
                v[2] = t
                vtx[nvtx] = v
                nvtx += 1
    return vtx[:nvtx]

def read_config():
    config = None
    with open('config.json', 'r') as fin:
        config = json.loads(fin.read());
    if config is None:
        print 'This script requires configuration file: config.json'
        exit(1)
    INPUT_DIR = config['prepare_data']['input_dir']
    OUTPUT_DIR = config['prepare_data']['output_dir']
    PATCH_SIZE_W = config['prepare_data']['patch_size_w']
    PATCH_SIZE_D = config['prepare_data']['patch_size_d']
    return INPUT_DIR, OUTPUT_DIR, PATCH_SIZE_W, PATCH_SIZE_D

