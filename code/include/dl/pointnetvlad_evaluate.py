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

from time import process_time

from plyfile import PlyData, PlyElement
import subprocess

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
parser.add_argument("--query_path", default=".")
parser.add_argument("--num_query", type=int, default=1)
parser.add_argument("--result_fdname", default=".")

parser.add_argument("--num_scenes", type=int, default=0) # for indoor dataset only

parser.add_argument("--prerun_dir", default="")

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

DAT_PATH = FLAGS.dat_path
RUN_ARG = FLAGS.run_arg
QUERY_PATH = FLAGS.query_path
NUM_QUERY = FLAGS.num_query

RESULT_FDNAME = FLAGS.result_fdname

NUM_SCENES = FLAGS.num_scenes

LOG_DIR = FLAGS.log_dir

PRERUN_DIR = FLAGS.prerun_dir
if not PRERUN_DIR:
    PRERUN_DIR = LOG_DIR

DATA_FOLDER = DAT_PATH + "/db_pointcloud_" + RUN_ARG + "/"
DATA_FILE = DATA_FOLDER + "point_cloud.bin"

REMESH_FILE = DAT_PATH + "/remesh_" + str(NUM_POINTS + 1) + ".mlx"
REMESH_FILE_BACKUP = DAT_PATH + "/remesh_" + str(NUM_POINTS + 11) + ".mlx"
MESHLABSERVER_PROG = "/home/meshlab-2020.06/meshlabserver"

QUERY_FILES = []
RESULT_FOLDERS = []
for i in range(NUM_QUERY):
    QUERY_FILES.append(QUERY_PATH + "/" + str(i) + "/query.ply")
    RESULT_FOLDERS.append(QUERY_PATH + "/" + str(i) + "/" + RESULT_FDNAME)

MODEL_FILE = FLAGS.model_filename + ".ckpt"

def load_center_point_ids(num_scenes):
    cp_list = []
    for i in range(num_scenes):
        cp_list.append(np.fromfile(os.path.join(DATA_FOLDER, "center_point_ids." + str(i)), dtype=np.int32))
    return cp_list

CP_LIST = load_center_point_ids(NUM_SCENES)
# print(CP_LIST[0])
# print(CP_LIST[1])
# print(CP_LIST[2])

def get_id_in_scene(cp_id):
    ending_id = 0
    for i in range(len(CP_LIST)):
        if (cp_id < ending_id + len(CP_LIST[i])):
            return i, CP_LIST[i][cp_id - ending_id]
        ending_id += len(CP_LIST[i])
    return -1, -1

def load_db_vectors():
    num_obj = 0
    db_vectors = []
    with open(os.path.join(PRERUN_DIR, "prerun.txt"), "r") as handle:
        num_obj = int(handle.readline())
        for i in range(num_obj):
            line = handle.readline()
            db_vectors.append(np.fromstring(line[1:-1], dtype=float, sep=", "))
    return db_vectors, num_obj

DATABASE_VECTORS, NUM_OBJ = load_db_vectors()

def load_query_pcs(filenames, num_query, num_pts):
    ret_pcs = np.array([], dtype=np.float32)
    for filename in filenames:
        query_plyfile_remesh = filename + ".rem.ply"
        subprocess.run([MESHLABSERVER_PROG, "-i", filename, "-o", query_plyfile_remesh, "-l", "x", "-s", REMESH_FILE])
        while not os.path.exists(query_plyfile_remesh):
            continue

        query_plydata_remesh = PlyData.read(query_plyfile_remesh)

        backup_flag = False
        if (len(query_plydata_remesh["vertex"]) < num_pts):
            query_plyfile_remesh_backup = query_plyfile_remesh + ".backup.ply"
            subprocess.run([MESHLABSERVER_PROG, "-i", query_plyfile_remesh, "-o", query_plyfile_remesh_backup, "-l", "x", "-s", REMESH_FILE_BACKUP])
            while not os.path.exists(query_plyfile_remesh_backup):
                continue
            query_plydata_remesh = PlyData.read(query_plyfile_remesh_backup)
            backup_flag = True

        for j in range(num_pts):
            for k in range(3):
                ret_pcs = np.append(ret_pcs, query_plydata_remesh["vertex"][j][k])

        subprocess.run(["rm", query_plyfile_remesh])
        if backup_flag:
            subprocess.run(["rm", query_plyfile_remesh_backup])

    ret_pcs = np.reshape(ret_pcs, (num_query, num_pts, 3))

    return ret_pcs

QUERY_POINT_CLOUDS = load_query_pcs(QUERY_FILES, NUM_QUERY, NUM_POINTS)

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


        # print(len(DATABASE_VECTORS))
        # for i in range(NUM_OBJ):
        #     print(DATABASE_VECTORS[i])

        database_nbrs = KDTree(DATABASE_VECTORS)

        output_size = NUM_OBJ//10
        if (NUM_SCENES > 0): # for indoor datasets, we output less number, since they are easier to reach the correct results
            output_size = output_size // 2

        query_start = process_time()
        
        QUERY_VECTORS=[]

        for i in range(NUM_QUERY):
            print("Processing query #" + str(i))
            QUERY_VECTORS.append(get_latent_vector_query(sess, ops, QUERY_POINT_CLOUDS[i]))
        distances, indices = database_nbrs.query(np.array(QUERY_VECTORS), k=output_size)

        # print(QUERY_VECTORS[0])
        # print(str(distances))
        # print(str(indices))

        query_end = process_time()
        each_query_time = (query_end - query_start) / NUM_QUERY
        print("Time for processing each query: " + str(each_query_time))

        for i in range(NUM_QUERY):
            if not os.path.exists(RESULT_FOLDERS[i]): os.mkdir(RESULT_FOLDERS[i])

            with open(os.path.join(RESULT_FOLDERS[i], "stat.txt"), "w") as handle:
                handle.write("user_time=" + str(each_query_time))

            if (NUM_SCENES > 0): # for indoor dataset
                with open(os.path.join(RESULT_FOLDERS[i], "output.raw.txt"), "w") as handle:
                    handle.write(str(output_size) + "\n")
                    for j in range(output_size):
                        index_obj = indices[i][j]
                        index_scene, index_in_scene = get_id_in_scene(index_obj)
                        handle.write(str(index_obj) + " " + str(index_scene) + " " + str(index_in_scene) + "\n")
            else: # for non-indoor dataset
                with open(os.path.join(RESULT_FOLDERS[i], "output.txt"), "w") as handle:
                    handle.write(str(output_size) + "\n")
                    for j in range(output_size):
                        handle.write(str(indices[i][j]) + " 1\n") # 1 is to indicate that this output is always valid


def get_latent_vector_query(sess, ops, query_pc):
    # print("Query shape(0):", query_pc.shape)
    queries = np.array([query_pc])
    queries = np.expand_dims(queries, axis=1)
    # print("Query shape(1):", queries.shape)

    # fake_queries=np.zeros((BATCH_NUM_QUERIES-1,1,NUM_POINTS,3)) # not needed because BATCH_NUM_QUERIES is 1
    fake_pos = np.zeros((BATCH_NUM_QUERIES, POSITIVES_PER_QUERY, NUM_POINTS, 3))
    fake_neg = np.zeros((BATCH_NUM_QUERIES, NEGATIVES_PER_QUERY, NUM_POINTS, 3))
    # q = np.vstack((queries,fake_queries))

    feed_dict = { ops["query"]:queries, ops["positives"]:fake_pos, ops["negatives"]:fake_neg, ops["is_training_pl"]:False }
    output = sess.run(ops["q_vec"], feed_dict=feed_dict)
    # print("Sess output shape:", output.shape)

    q_output = output[0]
    q_output = np.squeeze(output)

    # print("Query output shape:", q_output.shape)
    return q_output


if __name__ == "__main__":
    evaluate()
