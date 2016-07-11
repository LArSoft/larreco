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

batch_size = 256
nb_classes = 3
nb_epoch = 1000 #1000

#       name,      nflt1, kernel1, maxpool, nflt2, kernel2, convactfn, drop1, densesize, actfn, drop2 -> softmax
parameters = np.array([
#      ['deep',      64,     5,        0,     32,      7,     'relu',   0.2,      64,     'tanh', 0.2]
       ['small1_sgd_lorate_8k_coll',    32,     5,        0,      0,      7,     'relu',   0.2,     128,     'tanh', 0.2]
#      ['small1_ada_def',    64,     5,        0,      0,      7,     'relu',   0.2,     128,     'relu', 0.2]
#    , ['large',     64,     5,        1,     64,      7,     'relu',   0.2,     128,     'relu', 0.2]
#    , ['leakytanh', 64,     5,        1,     64,      5,     'leaky',  0.2,      64,     'tanh', 0.1]
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

_, CNN_INPUT_DIR, PATCH_SIZE = read_config()

# read train and test sets
X_train = None
Y_train = None
X_test = None
Y_test = None
filesX = [f for f in os.listdir(CNN_INPUT_DIR) if 'db_x' in f]
for fnameX in filesX:
    fnameY = fnameX.replace('_x_', '_y_')
    dataX = np.load(CNN_INPUT_DIR + '/' + fnameX)
    dataY = np.load(CNN_INPUT_DIR + '/' + fnameY)
    
    if '5000' in fnameX:
        X_test = dataX
        Y_test = dataY        
        continue
    
    if X_train is None:
        X_train = dataX
        Y_train = dataY
    else:
        X_train = np.concatenate((X_train, dataX))
        Y_train = np.concatenate((Y_train, dataY))

print 'Train', X_train.shape, 'test', X_test.shape

# input image dimensions
img_rows, img_cols = PATCH_SIZE, PATCH_SIZE

X_train = X_train.reshape(X_train.shape[0], 1, img_rows, img_cols)
X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
X_train = X_train.astype("float32")
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
    
    #if nb_conv2 > 0:
    #    if maxpool == 1:
    #        model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    #    model.add(Convolution2D(nb_filters2, nb_conv2, nb_conv2))
    #    if convactfn == 'leaky':
    #        model.add(LeakyReLU())
    #    else:
    #        model.add(Activation(convactfn))
    
    # model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
    model.add(Dropout(drop1))
    model.add(Flatten())
    
    # dense layers
    model.add(Dense(densesize))
    model.add(Activation(actfn))
    model.add(Dropout(drop2))
    model.add(Dense(nb_classes))
    model.add(Activation('softmax'))
    #model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])

    sgd = SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='categorical_crossentropy', optimizer=sgd)

    print('Fit config:', cfg_name)

    model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
                verbose=1, validation_data=(X_test, Y_test))
    score = model.evaluate(X_test, Y_test, verbose=0)
    print('Test score:', score)

    save_model(model, cfg_name)

print('All done!')
