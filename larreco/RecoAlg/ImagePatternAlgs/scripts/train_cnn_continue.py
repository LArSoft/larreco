import argparse
parser = argparse.ArgumentParser(description='Run CNN training on patches with a few different hyperparameter sets.')
parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
parser.add_argument('-m', '--model', help="input CNN model name (saved in JSON and h5 files)", default='cnn_model')
parser.add_argument('-o', '--output', help="output CNN model name (saved in JSON and h5 files)", default='cnn_model_out')
parser.add_argument('-g', '--gpu', help="Which GPU index", default='0')
args = parser.parse_args()

import theano.sandbox.cuda
theano.sandbox.cuda.use('gpu'+args.gpu)

import os
os.environ['KERAS_BACKEND'] = "theano"

import keras
if keras.__version__[0] != '1':
    print 'Please use Keras 1.x.x API due to matrix shape constraints in LArSoft interface'
    quit()
keras.backend.set_image_dim_ordering('th')

import numpy as np
np.random.seed(2017)  # for reproducibility
from keras.datasets import mnist
from keras.models import model_from_json
from keras.optimizers import SGD
from keras.utils import np_utils
from os.path import exists, isfile, join
import json

from utils import read_config, get_patch_size, count_events, shuffle_in_place

def load_model(name):
    with open(name + '_architecture.json') as f:
        model = model_from_json(f.read())
    model.load_weights(name + '_weights.h5')
    return model

def save_model(model, name):
    try:
        with open(name + '_architecture.json', 'w') as f:
            f.write(model.to_json())
        model.save_weights(name + '_weights.h5', overwrite=True)
        return True   # Save successful
    except:
        return False  # Save failed

#######################  configuration  #############################
print 'Reading configuration...'

config = read_config(args.config)
cfg_name = args.model
out_name = args.output

CNN_INPUT_DIR = config['training_on_patches']['input_dir']
# input image dimensions
PATCH_SIZE_W, PATCH_SIZE_D = get_patch_size(CNN_INPUT_DIR)
img_rows, img_cols = PATCH_SIZE_W, PATCH_SIZE_D

batch_size = config['training_on_patches']['batch_size']
nb_epoch = config['training_on_patches']['nb_epoch']
# nb_classes = config['training_on_patches']['nb_classes']

######################  CNN commpilation  ###########################
print 'Compiling CNN model...'

model = load_model(cfg_name)

#model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])

sgd = SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy', optimizer=sgd)

#######################  read data sets  ############################
n_training = count_events(CNN_INPUT_DIR, 'training')
X_train = np.zeros((n_training, 1, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
Y_train = np.zeros((n_training, nb_classes), dtype=np.int32)
print 'Training data size:', n_training, 'events; patch size:', PATCH_SIZE_W, 'x', PATCH_SIZE_D

ntot = 0
subdirs = [f for f in os.listdir(CNN_INPUT_DIR) if 'training' in f]
subdirs.sort()
for dirname in subdirs:
    print 'Reading data in', dirname
    filesX = [f for f in os.listdir(CNN_INPUT_DIR + '/' + dirname) if '_x.npy' in f]
    for fnameX in filesX:
        print '...training data', fnameX
        fnameY = fnameX.replace('_x.npy', '_y.npy')
        dataX = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameX)
        if dataX.dtype != np.dtype('float32'):
            dataX = dataX.astype("float32")
        dataY = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameY)
        n = dataY.shape[0]
        X_train[ntot:ntot+n] = dataX.reshape(n, 1, img_rows, img_cols)
        Y_train[ntot:ntot+n] = dataY
        ntot += n
print ntot, 'events ready'

n_testing = count_events(CNN_INPUT_DIR, 'testing')
X_test = np.zeros((n_testing, 1, PATCH_SIZE_W, PATCH_SIZE_D), dtype=np.float32)
Y_test = np.zeros((n_testing, nb_classes), dtype=np.int32)
print 'Testing data size:', n_testing, 'events'

ntot = 0
subdirs = [f for f in os.listdir(CNN_INPUT_DIR) if 'testing' in f]
subdirs.sort()
for dirname in subdirs:
    print 'Reading data in', dirname
    filesX = [f for f in os.listdir(CNN_INPUT_DIR + '/' + dirname) if '_x.npy' in f]
    for fnameX in filesX:
        print '...testing data', fnameX
        fnameY = fnameX.replace('_x.npy', '_y.npy')
        dataX = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameX)
        if dataX.dtype != np.dtype('float32'):
            dataX = dataX.astype("float32")
        dataY = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameY)
        n = dataY.shape[0]
        X_test[ntot:ntot+n] = dataX.reshape(n, 1, img_rows, img_cols)
        Y_test[ntot:ntot+n] = dataY
        ntot += n
print ntot, 'events ready'

dataX = None
dataY = None

print 'Shuffle training set...'
shuffle_in_place(X_train, Y_train)

print 'Training', X_train.shape, 'testing', X_test.shape

##########################  training  ###############################
print 'Fit config:', out_name
model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
            verbose=1, validation_data=(X_test, Y_test))

X_train = None
Y_train = None

score = model.evaluate(X_test, Y_test, verbose=0)
print('Test score:', score)

X_test = None
Y_test = None
#####################################################################

if save_model(model, out_name):
    print('All done!')
else:
    print('Error: model not saved.')

