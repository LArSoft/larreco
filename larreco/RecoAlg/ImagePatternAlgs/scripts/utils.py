import numpy as np
from os import listdir
from os.path import isfile, join
import os, json

def get_event_bounds(A, drift_margin = 0):
    # get center with 99% of signal
    cum = np.cumsum(np.sum(A, axis=0))
    start_ind = np.max([0, np.where(cum > cum[-1]*0.005)[0][0] - drift_margin])
    end_ind = np.min([A.shape[1], np.where(cum > cum[-1]*0.995)[0][0] + drift_margin])
    return start_ind, end_ind

def get_data(fname, drift_margin = 0, crop = True):
    print 'Reading', fname + '.raw'
    try:
        A_raw     = np.genfromtxt(fname + '.raw', delimiter=' ', dtype=np.float32)
        A_deposit = np.genfromtxt(fname + '.deposit', delimiter=' ', dtype=np.float32)
        A_pdg     = np.genfromtxt(fname + '.pdg', delimiter=' ', dtype=np.int32)
    except Exception as exc:
        print 'Bad event, return empty arrays'
        return None, None, None, None, None

    if A_raw.shape[0] < 8 or A_raw.shape[1] < 8: return None, None, None, None, None

    test_pdg = np.sum(A_pdg)
    test_dep = np.sum(A_deposit)
    test_raw = np.sum(A_raw)
    if test_raw == 0.0 or test_dep == 0.0 or test_pdg == 0: return None, None, None, None, None

    print test_raw, test_dep, test_pdg
    #assert np.sum(A_deposit) > 0
    # get main event body (99% signal)
    if crop:
        evt_start_ind, evt_stop_ind = get_event_bounds(A_deposit, drift_margin)
        A_raw     = A_raw[:,evt_start_ind:evt_stop_ind]
        A_deposit = A_deposit[:,evt_start_ind:evt_stop_ind]
        A_pdg     = A_pdg[:,evt_start_ind:evt_stop_ind]
    else:
        evt_start_ind = 0
        evt_stop_ind = A_raw.shape[1]
    print evt_start_ind, evt_stop_ind

    deposit_th_ind = A_deposit < 2.0e-5
    A_pdg[deposit_th_ind] = 0
    tracks = A_pdg.copy()
    showers = A_pdg.copy()
    tracks[(A_pdg & 0x0FFF) == 11] = 0
    tracks[tracks > 0]   = 1
    showers[(A_pdg & 0x0FFF) != 11] = 0
    showers[showers > 0] = 1
    return A_raw, A_deposit, A_pdg, tracks, showers

def get_patch(a, wire, drift, wsize, dsize):
    halfSizeW = wsize / 2;
    halfSizeD = dsize / 2;

    w0 = wire - halfSizeW;
    w1 = wire + halfSizeW;

    d0 = drift - halfSizeD;
    d1 = drift + halfSizeD;

    patch = np.zeros((wsize, dsize), dtype=np.float32)

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

def get_vertices(A):
    max_count = A.shape[0]*A.shape[1] / 4 # rather not more than 25% of plane filled with vertices
    vtx = np.zeros((max_count, 3), dtype=np.int32)
    nvtx = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if nvtx >= max_count: break
            if (A[i,j] & 0xFF000000) > 0:
                t = A[i,j] >> 24
                v = np.zeros(3)
                v[0] = i
                v[1] = j
                v[2] = t
                vtx[nvtx] = v
                nvtx += 1
    return vtx[:nvtx]

def get_nu_vertices(A):
    max_count = 10 # 10 vertices per view shoud be enough...
    vtx = np.zeros((max_count, 3), dtype=np.int32)
    nvtx = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if nvtx >= max_count: break
            if (A[i,j] & 0xFF0000) > 0:
                t = (A[i,j] >> 16) & 0xFF
                v = np.zeros(3)
                v[0] = i
                v[1] = j
                v[2] = t
                vtx[nvtx] = v
                nvtx += 1
    return vtx[:nvtx]

def read_config(cfgname):
    config = None
    with open(cfgname, 'r') as fin:
        config = json.loads(fin.read());
    if config is None:
        print 'This script requires configuration file: config.json'
        exit(1)
    return config

