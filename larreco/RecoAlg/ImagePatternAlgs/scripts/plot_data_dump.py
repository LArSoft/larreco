import argparse

from sys import argv
import numpy as np
import matplotlib.pyplot as plt

import ROOT
from ROOT import TFile
from root_numpy import hist2array

from utils import get_vertices

def main(argv):

    parser = argparse.ArgumentParser(description='Plot data dumped to ROOT file.')
    parser.add_argument('-i', '--input', help="Input file", default='datadump_hist.root')
    parser.add_argument('-e', '--event', help="Event index", default='0')
    parser.add_argument('-v', '--view', help="view index", default='')
    parser.add_argument('-m', '--module', help="LArSoft module name", default='datadump')
    args = parser.parse_args()

    if args.view != '': view_base = 'view_' + args.view
    else: view_base = ''

    file0 = TFile(args.input)
    raw_keys = [k.GetName() for k in file0.Get(args.module).GetListOfKeys() if view_base + '_raw' in k.GetName()]
    pdg_keys = [k.GetName() for k in file0.Get(args.module).GetListOfKeys() if view_base + '_pdg' in k.GetName()]
    dep_keys = [k.GetName() for k in file0.Get(args.module).GetListOfKeys() if view_base + '_deposit' in k.GetName()]
    
    ev_idx = int(args.event)
    if ev_idx >= len(raw_keys):
        print 'event index out of range'
        return
    print len(raw_keys), 'keys in the file, reading event key:', raw_keys[ev_idx]

    dep = hist2array(file0.Get(args.module + '/' + dep_keys[ev_idx]))
    raw = hist2array(file0.Get(args.module + '/' + raw_keys[ev_idx]))
    mc  = hist2array(file0.Get(args.module + '/' + pdg_keys[ev_idx]))
    pdg = mc & 0xFF

    vtx_list = get_vertices(mc)

    vtx_hadr = vtx_list[:,2] & 1
    vtx_pi0 = (vtx_list[:,2] >> 1) & 1
    vtx_decay = (vtx_list[:,2] >> 2) & 1
    vtx_conv = (vtx_list[:,2] >> 3) & 1
    vtx_eend = (vtx_list[:,2] >> 4) & 1
    vtx_list[:,2] = vtx_hadr + 2*vtx_pi0 + 3*vtx_decay + 4*vtx_conv + 5*vtx_eend
    print 'all vtx:', np.count_nonzero(vtx_list[:,2]), \
          'hadr:', np.count_nonzero(vtx_hadr), \
          'pi0:', np.count_nonzero(vtx_pi0), \
          'decay:', np.count_nonzero(vtx_decay), \
          'conv:', np.count_nonzero(vtx_conv), \
          'eend:', np.count_nonzero(vtx_eend)

    fig, ax = plt.subplots(1, 3, figsize=(36, 10))

    cs = ax[0].pcolor(np.transpose(pdg), cmap='gist_ncar')
    ax[0].scatter(vtx_list[:,0]+0.5, vtx_list[:,1]+0.5, s=50, c=vtx_list[:,2], cmap='rainbow', alpha=0.75)
    ax[0].set_title('PDG')
    fig.colorbar(cs, ax=ax[0])

    cs = ax[1].pcolor(np.transpose(dep), cmap='jet')
    ax[1].set_title('MC truth deposit')
    fig.colorbar(cs, ax=ax[1])

    cs = ax[2].pcolor(np.transpose(raw), cmap='jet')
    ax[2].scatter(vtx_list[:,0]+0.5, vtx_list[:,1]+0.5, s=50, c=vtx_list[:,2], cmap='rainbow', alpha=0.75)
    ax[2].set_title('ADC')
    fig.colorbar(cs, ax=ax[2])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main(argv)

