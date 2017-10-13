import argparse
parser = argparse.ArgumentParser(description='Run CNN training on patches with a few different hyperparameter sets.')
parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
parser.add_argument('-m', '--model', help="input CNN model name (saved in JSON and h5 files)", default='cnn_model')
parser.add_argument('-o', '--output', help="output CNN model name (saved in JSON and h5 files)", default='cnn_model_out')
parser.add_argument('-g', '--gpu', help="Which GPU index", default='0')
args = parser.parse_args()

import os
os.environ['KERAS_BACKEND'] = "tensorflow"
os.environ["CUDA_VISIBLE_DEVICES"]=args.gpu

import tensorflow as tf
import keras
if keras.__version__[0] != '2':
    print 'Please use the newest Keras 2.x.x API with the Tensorflow backend'
    quit()
keras.backend.set_image_data_format('channels_last')
keras.backend.set_image_dim_ordering('tf')

import numpy as np
np.random.seed(2017)  # for reproducibility
from keras.preprocessing.image import ImageDataGenerator
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
nb_classes = config['training_on_patches']['nb_classes']

######################  CNN commpilation  ###########################
print 'Compiling CNN model...'
with tf.device('/gpu:' + args.gpu):
    model = load_model(cfg_name)

    sgd = SGD(lr=0.002, decay=1e-5, momentum=0.9, nesterov=True)
    model.compile(optimizer=sgd,
                  loss={'em_trk_none_netout': 'categorical_crossentropy', 'michel_netout': 'mean_squared_error'},
                  loss_weights={'em_trk_none_netout': 0.1, 'michel_netout': 1.0})

#######################  read data sets  ############################
n_training = count_events(CNN_INPUT_DIR, 'training')
X_train = np.zeros((n_training, PATCH_SIZE_W, PATCH_SIZE_D, 1), dtype=np.float32)
EmTrkNone_train = np.zeros((n_training, 3), dtype=np.int32)
Michel_train = np.zeros((n_training, 1), dtype=np.int32)
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
        X_train[ntot:ntot+n] = dataX.reshape(n, img_rows, img_cols, 1)
        EmTrkNone_train[ntot:ntot+n] = dataY[:,[0, 1, 3]]
        Michel_train[ntot:ntot+n] = dataY[:,[2]]
        ntot += n
print ntot, 'events ready'

n_testing = count_events(CNN_INPUT_DIR, 'testing')
X_test = np.zeros((n_testing, PATCH_SIZE_W, PATCH_SIZE_D, 1), dtype=np.float32)
EmTrkNone_test = np.zeros((n_testing, 3), dtype=np.int32)
Michel_test = np.zeros((n_testing, 1), dtype=np.int32)
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
        X_test[ntot:ntot+n] = dataX.reshape(n, img_rows, img_cols, 1)
        EmTrkNone_test[ntot:ntot+n] = dataY[:,[0, 1, 3]]
        Michel_test[ntot:ntot+n] = dataY[:,[2]]
        ntot += n
print ntot, 'events ready'

dataX = None
dataY = None

print 'Training', X_train.shape, 'testing', X_test.shape

##########################  training  ###############################
datagen = ImageDataGenerator(
                featurewise_center=False, samplewise_center=False,
                featurewise_std_normalization=False,
                samplewise_std_normalization=False,
                zca_whitening=False,
                rotation_range=0, width_shift_range=0, height_shift_range=0,
                horizontal_flip=True, # randomly flip images
                vertical_flip=True)  # randomly flip images
datagen.fit(X_train)

def generate_data_generator(generator, X, Y1, Y2, b):
    genY1 = generator.flow(X, Y1, batch_size=b, seed=7)
    genY2 = generator.flow(X, Y2, batch_size=b, seed=7)
    while True:
            g1 = genY1.next()
            g2 = genY2.next()
            yield {'main_input': g1[0]}, {'em_trk_none_netout': g1[1], 'michel_netout': g2[1]}

print 'Fit config:', cfg_name
h = model.fit_generator(
              generate_data_generator(datagen, X_train, EmTrkNone_train, Michel_train, b=batch_size),
              validation_data=(
                  {'main_input': X_test},
                  {'em_trk_none_netout': EmTrkNone_test, 'michel_netout': Michel_test}),
              steps_per_epoch=X_train.shape[0]/batch_size, epochs=nb_epoch,
              verbose=1)

X_train = None
EmTrkNone_train = None
Michel_train = None

score = model.evaluate({'main_input': X_test},
                       {'em_trk_none_netout': EmTrkNone_test, 'michel_netout': Michel_test},
                       verbose=0)
print('Test score:', score)

X_test = None
EmTrkNone_test = None
Michel_test = None
#####################################################################

print h.history['loss']
print h.history['val_loss']

if save_model(model, args.output + cfg_name):
    print('All done!')
else:
    print('Error: model not saved.')

