import argparse
parser = argparse.ArgumentParser(description='Save model to TF protobuf')
parser.add_argument('-m', '--model', help="Keras TF model", default='model')
parser.add_argument('-o', '--output', help="TF graph", default='tf_graph.proto')
parser.add_argument('-n', '--nname', help="Use node as output if name contains", default='cnn_output')
parser.add_argument('-g', '--gpu', help="Which GPU index", default='0')
args = parser.parse_args()

from keras import backend as K
from keras.models import model_from_json
from keras.optimizers import SGD

from tensorflow.python.framework.graph_util import convert_variables_to_constants
import tensorflow as tf

import os
os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu

def load_model(name):
    with open(name + '_architecture.json') as f:
        model = model_from_json(f.read())
    model.load_weights(name + '_weights.h5')
    return model

K.set_learning_phase(0)
m = load_model(args.model)

nnames = [n.name for n in K.get_session().graph.as_graph_def().node]

lastnode = None
currentname = ''
output_nodes = []
for i, n in enumerate(nnames):
    basename = n[:n.find('/')]
    print i, n, basename
    if ('_netout' in n) or (args.nname in n):
        if (not lastnode is None) and (basename != currentname):
            output_nodes.append(lastnode)
        currentname = basename
        lastnode = n

if not lastnode is None:
    output_nodes.append(lastnode) # last of last output

if len(output_nodes) == 0:
    print 'Cannot find output node'
    exit(1)
print 'Output node names found:', output_nodes

minimal_graph = convert_variables_to_constants(K.get_session(), K.get_session().graph.as_graph_def(), output_nodes)

nnames = [n.name for n in minimal_graph.node]
print nnames

tf.train.write_graph(minimal_graph, '.', args.output + '.pb', as_text=False)

print 'all done!'
