import argparse
parser = argparse.ArgumentParser(description='Run CNN training on patches with a few different hyperparameter sets.')
parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
parser.add_argument('-o', '--output', help="Output model file name", default='model')
parser.add_argument('-g', '--gpu', help="Which GPU index", default='0')
args = parser.parse_args()

import theano.sandbox.cuda
theano.sandbox.cuda.use('gpu'+args.gpu)

import keras
if keras.__version__[0] > 1:
    print 'Please use Keras 1.x.x API due to matrix shape constraints in LArSoft interface'
    quit()

import numpy as np
np.random.seed(2017)  # for reproducibility
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers.advanced_activations import LeakyReLU
from keras.optimizers import SGD
from keras.utils import np_utils
from os.path import exists, isfile, join
import os, json

from utils import read_config

config = read_config(args.config)

CNN_INPUT_DIR = config['training_on_patches']['input_dir']

batch_size = config['training_on_patches']['batch_size']
nb_classes = config['training_on_patches']['nb_classes']
nb_epoch = config['training_on_patches']['nb_epoch']

#       name,                 nflt1, kernel1, maxpool, nflt2, kernel2, convactfn, drop1, densesize,  actfn, drop2 -> softmax
parameters = np.array([
#      ['deep',                 64,     5,        0,     32,      7,     'relu',   0.2,      64,     'tanh', 0.2]
       ['small1_sgd_lorate',    48,     5,        0,      0,      7,     'relu',   0.2,     128,     'tanh', 0.2]
#      ['small1_ada_def',       64,     5,        0,      0,      7,     'relu',   0.2,     128,     'relu', 0.2]
#    , ['large',                64,     5,        1,     64,      7,     'relu',   0.2,     128,     'relu', 0.2]
#    , ['leakytanh',            64,     5,        1,     64,      5,     'leaky',  0.2,      64,     'tanh', 0.1]
])


def save_model(model, name):
    try:
        with open(name + '_architecture.json', 'w') as f:
            f.write(model.to_json())
        model.save_weights(name + '_weights.h5', overwrite=True)
        return True  # Save successful
    except:
        print 'save failed' #sys.exc_info()  # Prints exceptions
        return False  # Save failed

def shuffle_in_place(a, b):
    assert len(a) == len(b)
    rng_state = np.random.get_state()
    np.random.shuffle(a)
    np.random.set_state(rng_state)
    np.random.shuffle(b)


# read train and test sets
X_train = None
Y_train = None
X_test = None
Y_test = None

subdirs = [f for f in os.listdir(CNN_INPUT_DIR) if 'training' in f]
subdirs.sort()
for dirname in subdirs:
    print 'Reading data in', dirname
    filesX = [f for f in os.listdir(CNN_INPUT_DIR + '/' + dirname) if '_x.npy' in f]
    for fnameX in filesX:
        print '...training data', fnameX
        fnameY = fnameX.replace('_x.npy', '_y.npy')
        dataX = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameX)
        dataY = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameY)
        if nb_classes == 3 and dataY.shape[1] == 4:
            dataY = dataY[:, [0, 1, 3]] #skip column with Michel labels
        if X_train is None:
            X_train = dataX
            Y_train = dataY
        else:
            X_train = np.concatenate((X_train, dataX))
            Y_train = np.concatenate((Y_train, dataY))

subdirs = [f for f in os.listdir(CNN_INPUT_DIR) if 'testing' in f]
subdirs.sort()
for dirname in subdirs:
    print 'Reading data in', dirname
    filesX = [f for f in os.listdir(CNN_INPUT_DIR + '/' + dirname) if '_x.npy' in f]
    for fnameX in filesX:
        print '...testing data', fnameX
        fnameY = fnameX.replace('_x.npy', '_y.npy')
        dataX = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameX)
        dataY = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameY)
        if nb_classes == 3 and dataY.shape[1] == 4:
            dataY = dataY[:, [0, 1, 3]] #skip column with Michel labels
        if X_test is None:
            X_test = dataX
            Y_test = dataY
        else:
            X_test = np.concatenate((X_test, dataX))
            Y_test = np.concatenate((Y_test, dataY))

dataX = None
dataY = None

print 'Shuffle training set...'
shuffle_in_place(X_train, Y_train)

print 'Train', X_train.shape, 'test', X_test.shape

# input image dimensions
PATCH_SIZE_W = X_train.shape[1]
PATCH_SIZE_D = X_train.shape[2]
img_rows, img_cols = PATCH_SIZE_W, PATCH_SIZE_D

X_train = X_train.reshape(X_train.shape[0], 1, img_rows, img_cols)
if X_train.dtype != np.dtype('float32'):
    X_train = X_train.astype("float32")

X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
if X_test.dtype != np.dtype('float32'):
    X_test = X_test.astype("float32")

print('X_train shape:', X_train.shape)
print(X_train.shape[0], 'train samples')
print(X_test.shape[0], 'test samples')


for p in range(parameters.shape[0]):

    nb_pool = 2 # size of pooling area for max pooling

    cfg_name = parameters[p,0]

    nb_filters1 = int(parameters[p,1]) # number of convolutional filters to use
    nb_conv1 = int(parameters[p,2]) # convolution kernel size

    maxpool = int(parameters[p,3]) # max pooling between conv. layers

    nb_filters2 = int(parameters[p,4]) # number of convolutional filters to use
    nb_conv2 = int(parameters[p,5]) # convolution kernel size
    
    convactfn = parameters[p,6]
    
    drop1 = float(parameters[p,7])
    
    densesize = int(parameters[p,8])
    
    actfn = parameters[p,9]
    
    drop2 = float(parameters[p,10])

    # CNN definition
    model = Sequential()
    model.add(Convolution2D(nb_filters1, nb_conv1, nb_conv1,
                            border_mode='valid',
                            input_shape=(1, img_rows, img_cols)))
                            
    if convactfn == 'leaky':
        model.add(LeakyReLU())
    else:
        model.add(Activation(convactfn))
    
    if nb_conv2 > 0:
        if maxpool == 1:
            model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
        model.add(Convolution2D(nb_filters2, nb_conv2, nb_conv2))
        if convactfn == 'leaky':
            model.add(LeakyReLU())
        else:
            model.add(Activation(convactfn))
    
    # model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(drop1))
    model.add(Flatten())
    
    # dense layers
    model.add(Dense(densesize))
    model.add(Activation(actfn))
    model.add(Dropout(drop2))

    model.add(Dense(32))
    model.add(Activation(actfn))
    model.add(Dropout(drop2))

    # output
    model.add(Dense(nb_classes))

    sgd = SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
    if nb_classes > 3:
        model.add(Activation('sigmoid'))
        model.compile(loss='mean_squared_error', optimizer=sgd)
    else:
        model.add(Activation('softmax'))
        #model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])
        model.compile(loss='categorical_crossentropy', optimizer=sgd)

    print('Fit config:', cfg_name)

    model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
                verbose=1, validation_data=(X_test, Y_test))
    score = model.evaluate(X_test, Y_test, verbose=0)
    print('Test score:', score)

    save_model(model, args.output + cfg_name)

print('All done!')
