export PATH=/usr/local/python2/bin/:$PATH
export LD_LIBRARY_PATH=/usr/local/python2/lib/:$LD_LIBRARY_PATH
export CUDA_VISIBLE_DEVICES=`cat /tmp/pbs.prologue.$PBS_JOBID`
export CUDA_ROOT=/usr/local/cuda-7.5
export THEANO_FLAGS='cuda.root=/usr/local/cuda-7.5,device=gpu,floatX=float32'

cd /home/rsulej/CNN/scripts
python train_cnn.py

