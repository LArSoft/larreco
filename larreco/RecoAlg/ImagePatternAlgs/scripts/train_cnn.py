import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(2017)  # for reproducibility
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.utils import np_utils

from utils import get_data, read_config, prepare_signal

def save_model(model, name):
    try:
        with open(name + '_architecture.json', 'w') as f:
            f.write(model.to_json())
        model.save_weights(name + '_weights.h5', overwrite=True)
        return True  # Save successful
    except:
        print 'save failed' #sys.exc_info()  # Prints exceptions
        return False  # Save failed

TOTAL_EVENTS = 100
TEST_EVENTS_SPLIT = 70

_, CNN_INPUT_DIR, PATCH_SIZE = read_config()

db, db_y = None, None
test_split = 0
for evt_no in range(1,TOTAL_EVENTS+1):
    tmp_db   = np.load(CNN_INPUT_DIR+'/db_evt_'+str(evt_no)+'.npy')
    tmp_db_y = np.load(CNN_INPUT_DIR+'/db_y_evt_'+str(evt_no)+'.npy')

    if db is None:
        db = tmp_db
        db_y = tmp_db_y
    else:
        db   = np.concatenate((db, tmp_db))
        db_y = np.concatenate((db_y, tmp_db_y))
    if evt_no == TEST_EVENTS_SPLIT:
        test_split = db.shape[0]
print 'Total data', db.shape, db_y.shape, 'test_split', test_split
print 'Tracks', np.sum(db_y == 1), 'showers', np.sum(db_y == 0)

# split between tran and test sets
X_train = db[:test_split]
X_test  = db[test_split:]
y_train = db_y[:test_split]
y_test  = db_y[test_split:]
print 'Train', X_train.shape, 'test', X_test.shape


batch_size = 256 #128
nb_classes = int(2)
nb_epoch = 100 ## 1000
# input image dimensions
img_rows, img_cols = PATCH_SIZE, PATCH_SIZE
# number of convolutional filters to use
nb_filters = 64
# size of pooling area for max pooling
nb_pool = 2
# convolution kernel size
nb_conv = 5

X_train = X_train.reshape(X_train.shape[0], 1, img_rows, img_cols)
X_test = X_test.reshape(X_test.shape[0], 1, img_rows, img_cols)
X_train = X_train.astype("float32")
X_test = X_test.astype("float32")
print('X_train shape:', X_train.shape)
print(X_train.shape[0], 'train samples')
print(X_test.shape[0], 'test samples')
Y_train = np_utils.to_categorical(y_train, nb_classes)
Y_test = np_utils.to_categorical(y_test, nb_classes)

# CNN definition
model = Sequential()
model.add(Convolution2D(nb_filters, nb_conv, nb_conv,
                        border_mode='valid',
                        input_shape=(1, img_rows, img_cols)))
model.add(Activation('relu'))
model.add(Convolution2D(nb_filters, nb_conv, nb_conv))
model.add(Activation('relu'))
model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
model.add(Dropout(0.25))
model.add(Flatten())
# dense layers
model.add(Dense(64))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(nb_classes))
model.add(Activation('softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adadelta', metrics=['accuracy'])


model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
              verbose=1, validation_data=(X_test, Y_test))
score = model.evaluate(X_test, Y_test, verbose=0)
print('Test score:', score[0])
print('Test accuracy:', score[1])

save_model(model, '../cnn_models/test')

