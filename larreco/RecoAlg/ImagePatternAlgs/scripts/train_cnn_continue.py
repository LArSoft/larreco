import numpy as np
np.random.seed(2017)  # for reproducibility
from keras.datasets import mnist
from keras.models import model_from_json
from keras.optimizers import SGD
from keras.utils import np_utils
from os.path import exists, isfile, join
import os, json
import argparse

from utils import read_config

parser = argparse.ArgumentParser(description='Run CNN training on patches with a few different hyperparameter sets.')
parser.add_argument('-c', '--config', help="JSON with script configuration", default='config.json')
parser.add_argument('-m', '--model', help="input CNN model name (saved in JSON and h5 files)", default='cnn_model')
parser.add_argument('-o', '--output', help="output CNN model name (saved in JSON and h5 files)", default='cnn_model_out')
args = parser.parse_args()

config = read_config(args.config)
cfg_name = args.model
out_name = args.output

CNN_INPUT_DIR = config['training_on_patches']['input_dir']

batch_size = config['training_on_patches']['batch_size']
nb_classes = config['training_on_patches']['nb_classes']
nb_epoch = config['training_on_patches']['nb_epoch']

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
        return True  # Save successful
    except:
        print 'save failed' #sys.exc_info()  # Prints exceptions
        return False  # Save failed

# read train and test sets
X_train = None
Y_train = None
X_test = None
Y_test = None

n_test_dir = 1 # how many dirs used as testing data

subdirs = [f for f in os.listdir(CNN_INPUT_DIR) if '000' in f]
subdirs.sort()
for dirname in subdirs:
    print 'Reading data in', dirname
    filesX = [f for f in os.listdir(CNN_INPUT_DIR + '/' + dirname) if '_x.npy' in f]
    for fnameX in filesX:
        fnameY = fnameX.replace('_x.npy', '_y.npy')
        dataX = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameX)
        dataY = np.load(CNN_INPUT_DIR + '/' + dirname + '/' + fnameY)

        if n_test_dir > 0:
            print '...testing data', fnameX
            if X_test is None:
                X_test = dataX
                Y_test = dataY
            else:
                X_test = np.concatenate((X_test, dataX))
                Y_test = np.concatenate((Y_test, dataY))
        else:
            print '...training data', fnameX
            if X_train is None:
                X_train = dataX
                Y_train = dataY
            else:
                X_train = np.concatenate((X_train, dataX))
                Y_train = np.concatenate((Y_train, dataY))

    n_test_dir -= 1

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

model = load_model(cfg_name)

#model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])

sgd = SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
model.compile(loss='categorical_crossentropy', optimizer=sgd)

model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
            verbose=1, validation_data=(X_test, Y_test))
score = model.evaluate(X_test, Y_test, verbose=0)
print('Test score:', score)

save_model(model, out_name)

print('All done!')

