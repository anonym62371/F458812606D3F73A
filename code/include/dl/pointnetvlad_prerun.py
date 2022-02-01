import argparse
import math
import numpy as np
import tensorflow as tf
import socket
import importlib
import os
import sys
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

from pointnetvlad_cls import *
from loading_pointclouds import *
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KDTree

#params
parser = argparse.ArgumentParser()
parser.add_argument("--gpu", type=int, default=0, help="GPU to use [default: GPU 1]")
parser.add_argument("--log_dir", default="log/", help="Log dir [default: log]")
parser.add_argument("--positives_per_query", type=int, default=1, help="Number of potential positives in each training tuple [default: 2]")
parser.add_argument("--negatives_per_query", type=int, default=9, help="Number of definite negatives in each training tuple [default: 20]")
parser.add_argument("--batch_num_queries", type=int, default=1, help="Batch Size during training [default: 1]")
parser.add_argument("--dimension", type=int, default=256)
parser.add_argument("--decay_step", type=int, default=200000, help="Decay step for lr decay [default: 200000]")
parser.add_argument("--decay_rate", type=float, default=0.7, help="Decay rate for lr decay [default: 0.8]")

parser.add_argument("--dat_path", default=".")
parser.add_argument("--run_arg", default="15")
parser.add_argument("--num_points", type=int, default=128)
parser.add_argument("--model_filename", default="model")
parser.add_argument("--output_dir", default="")

FLAGS = parser.parse_args()

#BATCH_SIZE = FLAGS.batch_size
BATCH_NUM_QUERIES = FLAGS.batch_num_queries
EVAL_BATCH_SIZE = 1
# NUM_POINTS = 4096
NUM_POINTS = FLAGS.num_points
POSITIVES_PER_QUERY= FLAGS.positives_per_query
NEGATIVES_PER_QUERY= FLAGS.negatives_per_query
GPU_INDEX = FLAGS.gpu
DECAY_STEP = FLAGS.decay_step
DECAY_RATE = FLAGS.decay_rate

LOG_DIR = FLAGS.log_dir
DAT_PATH = FLAGS.dat_path
RUN_ARG = FLAGS.run_arg
DATA_FILE = DAT_PATH + "/db_pointcloud_" + RUN_ARG + "/point_cloud.bin"
MODEL_FILE = FLAGS.model_filename + ".ckpt"

OUTPUT_DIR = FLAGS.output_dir
if not OUTPUT_DIR:
    OUTPUT_DIR = LOG_DIR

DB_POINT_CLOUDS, NUM_OBJ = load_3dor_pc(DATA_FILE, NUM_POINTS)
print("Point clouds loaded of number", NUM_OBJ)

BN_INIT_DECAY = 0.5
BN_DECAY_DECAY_RATE = 0.5
BN_DECAY_DECAY_STEP = float(DECAY_STEP)
BN_DECAY_CLIP = 0.99

def get_bn_decay(batch):
    bn_momentum = tf.train.exponential_decay(
                      BN_INIT_DECAY,
                      batch*BATCH_NUM_QUERIES,
                      BN_DECAY_DECAY_STEP,
                      BN_DECAY_DECAY_RATE,
                      staircase=True)
    bn_decay = tf.minimum(BN_DECAY_CLIP, 1 - bn_momentum)
    return bn_decay     

def evaluate():
    with tf.Graph().as_default():
        with tf.device("/gpu:"+str(GPU_INDEX)):
            print("In Graph")
            query= placeholder_inputs(BATCH_NUM_QUERIES, 1, NUM_POINTS)
            positives= placeholder_inputs(BATCH_NUM_QUERIES, POSITIVES_PER_QUERY, NUM_POINTS)
            negatives= placeholder_inputs(BATCH_NUM_QUERIES, NEGATIVES_PER_QUERY, NUM_POINTS)
            eval_queries= placeholder_inputs(EVAL_BATCH_SIZE, 1, NUM_POINTS)

            is_training_pl = tf.placeholder(tf.bool, shape=())
            print(is_training_pl)

            batch = tf.Variable(0)
            bn_decay = get_bn_decay(batch)

            with tf.variable_scope("query_triplets") as scope:
                vecs= tf.concat([query, positives, negatives],1)
                print(vecs)                
                out_vecs= forward(vecs, is_training_pl, bn_decay=bn_decay)
                print(out_vecs)
                q_vec, pos_vecs, neg_vecs= tf.split(out_vecs, [1,POSITIVES_PER_QUERY,NEGATIVES_PER_QUERY],1)
                print(q_vec)
                print(pos_vecs)
                print(neg_vecs)

            saver = tf.train.Saver()
            
        # Create a session
        gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.95)
        config = tf.ConfigProto(gpu_options=gpu_options)
        config.gpu_options.allow_growth = True
        config.allow_soft_placement = True
        config.log_device_placement = False

        config.intra_op_parallelism_threads = 4
        config.inter_op_parallelism_threads = 1

        sess = tf.Session(config=config)


        saver.restore(sess, os.path.join(LOG_DIR, MODEL_FILE))
        print("Model restored.")

        ops = {"query": query,
               "positives": positives,
               "negatives": negatives,
               "is_training_pl": is_training_pl,
               "eval_queries": eval_queries,
               "q_vec":q_vec,
               "pos_vecs": pos_vecs,
               "neg_vecs": neg_vecs}

        db_vectors = []
        for i in range(NUM_OBJ):
            print("Get latent vector for DB obj #" + str(i), flush=True)
            db_vectors.append(get_latent_vector(sess, ops, i))

        with open(os.path.join(OUTPUT_DIR, "prerun.txt"), "w") as handle:
            handle.write(str(NUM_OBJ) + "\n")
            for i in range(NUM_OBJ):
                handle.write(str(db_vectors[i].tolist()) + "\n")

def get_latent_vector(sess, ops, index, query_vec=False):
    is_training=False

    if (query_vec):
        queries = [rotate_point_cloud(DB_POINT_CLOUDS[index])]
    else:
        queries = [DB_POINT_CLOUDS[index]]
    queries = np.array(queries)
    queries = np.expand_dims(queries, axis=1)
    # print("Shape reshaped:", queries.shape)

    # fake_queries=np.zeros((BATCH_NUM_QUERIES-1,1,NUM_POINTS,3)) # not needed because BATCH_NUM_QUERIES is 1
    fake_pos = np.zeros((BATCH_NUM_QUERIES, POSITIVES_PER_QUERY, NUM_POINTS, 3))
    fake_neg = np.zeros((BATCH_NUM_QUERIES, NEGATIVES_PER_QUERY, NUM_POINTS, 3))
    # q = np.vstack((queries,fake_queries))
    #print(queries.shape)

    feed_dict = { ops["query"]:queries, ops["positives"]:fake_pos, ops["negatives"]:fake_neg, ops["is_training_pl"]:is_training }
    output = sess.run(ops["q_vec"], feed_dict=feed_dict)
    #print(output.shape)

    q_output = output[0]
    q_output=np.squeeze(output)

    # print(q_output.shape)
    return q_output


if __name__ == "__main__":
    evaluate()
