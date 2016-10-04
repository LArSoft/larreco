import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pylab import *
from scipy import ndimage

from utils import get_data, get_patch

from keras.models import model_from_json
from keras.optimizers import SGD

def load_model(name):
    with open(name + '_architecture.json') as f:
        model = model_from_json(f.read())
    model.load_weights(name + '_weights.h5')
    return model

raw, deposit, pdg, tracks, showers = get_data('your_dir/raw_data_filename_no_extension')
total_learning_patches = int(np.sum(tracks) + np.sum(showers))

PATCH_SIZE_W = 32
PATCH_SIZE_D = 44

# input image dimensions
img_rows, img_cols = PATCH_SIZE_W, PATCH_SIZE_D

db = np.zeros((total_learning_patches, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
db_y = np.zeros((total_learning_patches, 3), dtype=np.int32)

cnt_ind = 0
for i in range(raw.shape[0]):
    for j in range(raw.shape[1]):
        target = np.zeros(3)
        if tracks[i,j] == 1:
            target[0] = 1
        elif showers[i,j] == 1:
            target[1] = 1
        else: continue

        db[cnt_ind] = get_patch(raw, i, j, PATCH_SIZE_W, PATCH_SIZE_D)
        db_y[cnt_ind] = target
        cnt_ind += 1
        
db = db[:(cnt_ind)]
print db.shape, cnt_ind

sgd = SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
m = load_model('your_dir/cnn_model_filename_no_extension')
m.compile(loss='categorical_crossentropy', optimizer=sgd)

pred = m.predict(db.reshape(db.shape[0], 1, img_rows, img_cols))

A_pred_inter = np.zeros(raw.shape)
A_pred_decay = np.zeros(raw.shape)

pred_thr = 0.9
cnt_ind = 0
for i in range(raw.shape[0]):
    for j in range(raw.shape[1]):
        if not (tracks[i,j] == 1 or showers[i,j] == 1):
            continue
        
        if cnt_ind < pred.shape[0]:
            A_pred_inter[i,j] = pred[cnt_ind,0]
            A_pred_decay[i,j] = pred[cnt_ind,1]
        
        cnt_ind += 1
        
print cnt_ind, A_pred_inter.shape, A_pred_decay.shape
        
fig, ax = subplots(2,2, figsize=(15, 15))
ax[0,0].imshow(-np.transpose(raw), cmap='gray', interpolation='none', aspect='auto')
ax[0,0].set_title('ADC')
ax[0,1].imshow(-np.transpose(showers), cmap='gray', interpolation='none', aspect='auto')
ax[0,1].set_title('SHOWERS LABEL')
ax[1,0].imshow(-np.transpose(A_pred_inter), cmap='gray', interpolation='none', aspect='auto')
ax[1,0].set_title('CNN out[0]')
ax[1,1].imshow(-np.transpose(A_pred_decay), cmap='gray', interpolation='none', aspect='auto')
ax[1,1].set_title('CNN out[1]')
plt.show()
