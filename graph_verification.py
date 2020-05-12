# -*- coding: utf-8 -*-
"""
Created on Mon May 11 19:03:13 2020

@author: zohre
"""

import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import Normalizer
from numpy import *
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import Normalizer
import Graph_Sampling 
from networkx.algorithms.approximation import min_weighted_vertex_cover
import time
from networkx.algorithms import has_path
#from matrix_completion import svt_solve, calc_unobserved_rmse
import scipy.sparse as sps
import itertools
import pandas as pds
import pickle
import math
import random
import numpy as np
import copy
## SVD
import seaborn as sns
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms import community


#g = nx.read_edgelist("./atp/dataset/fb.txt", create_using= nx.Graph(),nodetype=int)
g = nx.read_edgelist('./atp/dataset/government_edges.csv', delimiter=',',nodetype = int, create_using = nx.DiGraph())
num_edges=g.number_of_edges()
num_nodes=g.number_of_nodes()


#path = nx.single_source_shortest_path(g,0)

communities_generator = community.girvan_newman(g)
top_level_communities = next(communities_generator)
next_level_communities = next(communities_generator)
sorted(map(sorted, next_level_communities))
