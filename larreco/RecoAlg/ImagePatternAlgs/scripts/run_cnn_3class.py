import argparse
parser = argparse.ArgumentParser(description='Run CNN over a full 2D projection.')
parser.add_argument('-i', '--input', help="Input file", default='datadump_hist.root') # '/eos/user/r/rosulej/ProtoDUNE/datadump/datadump_hist.root'
parser.add_argument('-e', '--event', help="Event index", default='0')
parser.add_argument('-m', '--module', help="LArSoft module name", default='datadump')
parser.add_argument('-f', '--full', help="Full 2D plane (1), or not-empty pixels only (0)", default='1')
parser.add_argument('-n', '--net', help="Network model name (json + h5 files)", default='model') # /eos/user/r/rosulej/models/pdune_em-trk-michel_clean_iter350
parser.add_argument('-g', '--gpu', help="GPU index to use (default is CPU)", default='-1')
parser.add_argument('-r', '--rows', help="Patch rows (wires)", default='44')
parser.add_argument('-c', '--cols', help="Patch cols (ticks)", default='48')
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
from ROOT import TFile
from utils import get_data, get_patch

import theano
import theano.sandbox.cuda
if args.gpu == '-1':
    theano.sandbox.cuda.use('cpu')
    print 'Running on CPU, use -g option to say which device index should be used.'
else:
    theano.sandbox.cuda.use('gpu' + args.gpu)

import os
os.environ['KERAS_BACKEND'] = "theano"

import keras
from keras.models import model_from_json

print 'Software versions: Theano ', theano.__version__, ', Keras ', keras.__version__
if keras.backend.backend() != 'theano':
    print '**** You should be using Theano backend now...****'
    quit()
keras.backend.set_image_dim_ordering('th')

def load_model(name):
    with open(name + '_architecture.json') as f:
        model = model_from_json(f.read())
    model.load_weights(name + '_weights.h5')
    return model


PATCH_SIZE_W = int(args.rows) # wires
PATCH_SIZE_D = int(args.cols) # ticks (downsampled)
crop_event = False

rootModule = args.module
rootFile = TFile(args.input)
keys = [rootModule+'/'+k.GetName()[:-4] for k in rootFile.Get(rootModule).GetListOfKeys() if '_raw' in k.GetName()]
evname = keys[int(args.event)]

raw, deposit, pdg, tracks, showers = get_data(rootFile, evname, PATCH_SIZE_D/2 + 2, crop_event)
full2d = int(args.full)
if full2d == 1: total_patches = raw.size
else: total_patches = int(np.sum(tracks) + np.sum(showers))
print 'Number of pixels:', total_patches

inputs = np.zeros((total_patches, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)

cnt_ind = 0
for r in range(raw.shape[0]):
    for c in range(raw.shape[1]):
        if full2d == 0 and not(tracks[r, c] == 1 or showers[r, c] == 1):
            continue

        inputs[cnt_ind] = get_patch(raw, r, c, PATCH_SIZE_W, PATCH_SIZE_D)
        cnt_ind += 1

inputs = inputs[:(cnt_ind)]
print inputs.shape, cnt_ind

model_name = args.net
m = load_model(model_name)
m.compile(loss='mean_squared_error', optimizer='sgd')

print 'running CNN...'
pred = m.predict(inputs.reshape(inputs.shape[0], 1, PATCH_SIZE_W, PATCH_SIZE_D)) # nsamples, channel, rows, cols
print '...done', pred.shape

outputs = np.zeros((pred.shape[1], raw.shape[0], raw.shape[1]), dtype=np.float32)

mask = np.zeros(raw.shape, dtype=np.int32)
mask_thr = 0.67

none_idx = pred.shape[1] - 1 # here is "none" label usually
pnorm = pred[:, 0] + pred[:, 1] + pred[:, none_idx]

cnt_ind = 0
for r in range(outputs.shape[1]):
    for c in range(outputs.shape[2]):
        if full2d == 0 and not(tracks[r, c] == 1 or showers[r, c] == 1):
            continue

        if pnorm[cnt_ind] > 0:
            pn = 1.0 / pnorm[cnt_ind]
            outputs[0, r, c] = pred[cnt_ind, 0] * pn
            outputs[1, r, c] = pred[cnt_ind, 1] * pn

            outputs[none_idx, r, c] = 1 - pred[cnt_ind, none_idx] * pn
            if outputs[none_idx, r, c] > mask_thr:
                if outputs[1, r, c] > outputs[0, r, c]: mask[r, c] = 1
                else: mask[r, c] = -1

        cnt_ind += 1


fig, ax = plt.subplots(2, 3, figsize=(28, 14))

cs = ax[0,0].pcolor(np.transpose(pdg & 0xFF), cmap='gist_ncar')
ax[0,0].set_title('PDG')
fig.colorbar(cs, ax=ax[0,0])

cs = ax[0,1].pcolor(np.transpose(deposit), cmap='jet')
ax[0,1].set_title('MC truth')
fig.colorbar(cs, ax=ax[0,1])

cs = ax[0,2].pcolor(np.transpose(raw), cmap='jet')
ax[0,2].set_title('ADC')
fig.colorbar(cs, ax=ax[0,2])

cs = ax[1,0].pcolor(-np.transpose(outputs[0]), cmap='CMRmap')
ax[1,0].set_title('P(track-like)')
fig.colorbar(cs, ax=ax[1,0])

cs = ax[1,1].pcolor(-np.transpose(outputs[1]), cmap='CMRmap')
ax[1,1].set_title('P(EM-like)')
fig.colorbar(cs, ax=ax[1,1])

#cs = ax[1,2].contour(np.transpose(outputs[none_idx]), levels=[0, 0.33, 0.67, 1], colors=('black' , (0.6, 0.6, 1), (0.9, 0, 0), 'green'))
#ax[1,2].set_title('1 - P(empty)')
cs = ax[1,2].pcolor(np.transpose(mask), cmap='RdYlGn')
ax[1,2].set_title('ROI')
fig.colorbar(cs, ax=ax[1,2])

#plt.subplots_adjust(left=0.03, right=0.999, bottom=0.07, top=0.93)
plt.tight_layout()
plt.show()
