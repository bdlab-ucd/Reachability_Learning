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
from sklearn import decomposition, datasets, model_selection, preprocessing, metrics

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
import seaborn as sns
import networkx as nx
from sklearn.decomposition import NMF
from sklearn.metrics import mean_squared_error
from datetime import datetime
import psutil
import os




## set input parameters 
query_budget=100000000
rank=64
depth_limit_search=4
source_node=5
num_traverse=2
bfs_list=[]
dfs_list=[]
sources=[5,3406,3,10,16]
backbone=False
jump_prob=.2
burning_prob=.2
test_ratio=.04
from scipy.sparse.linalg import norm
##### input dataset


def get_result (reach):        
    model = NMF(n_components= 64, init='random', random_state=0)
    W1 = model.fit_transform(reach)
    H1 = model.components_
                           
    HT=np.transpose(H1)
    
    test_size=int(num_nodes*test_ratio)
    reach_test=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    #pred=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    
    nzeroo=np.argwhere(reach != 0)
    print(nzeroo)
    print(len(nzeroo))
    count=0
    while count<test_size:
        x=random.randint(0,test_size-1)
        y=random.randint(0,test_size-1)
        if g.has_node(x) and g.has_node(y): 
            if (x,y) not in nzeroo and has_path(g,x,y):
                reach_test[[x],[y]]=1 
                #pred[[x],[y]]=W1[x].dot(HT[y])
                  #print(W1[x].dot(HT[y]))
                count=count+1
                  
    count=0
    while count<test_size:
            x=random.randint(0,test_size-1)
            y=random.randint(0,test_size-1)
            if  g.has_node(x) and g.has_node(y): 
                if [[x],[y]] not in nzeroo:
                    reach_test[[x],[y]]=0
                    #pred[[x],[y]]=W1[x].dot(HT[y])
                    count=count+1              
                   
                  
    
    #pred_arr=pred.toarray()
    #reach_test_dense=reach_test.toarray()
    
    #res=mean_squared_error(reach_test_dense, pred_arr)
    res=model.reconstruction_err_
    fnorm=norm(reach, 'fro')
    print(res)
    print(fnorm)
    return str(res)+"__"+str(res/fnorm)


#'web-Google.txt','soc-LiveJournal1.txt','cit-patent.edges'
for dataa in ['soc-LiveJournal1.txt']:
    sources=[]
    g = nx.read_edgelist("./Reachability_Learning/dataset/"+dataa, create_using= nx.DiGraph(),nodetype=int)
    highdegree=sorted(g.degree, key=lambda x: x[1], reverse=True)
    sources.append(highdegree[0][0])
    sources.append(highdegree[1][0])
    
    f = open("lognew4_time_"+dataa, "a")
    #soc-LiveJournal1.txt
    #cit-patent.edges
    #g = nx.karate_club_graph()
    print(datetime.now())
    num_edges=g.number_of_edges()
    num_nodes=g.number_of_nodes()
    print(dataa)
    print("\n")
    print("number of nodes: {}".format(num_nodes) )

    print("number of edges: {}".format(num_edges) )

    s1=datetime.now()
   
    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)

    min_vertex=min_weighted_vertex_cover(g)
    count=0
    while count<query_budget:
        i=random.choice(list(min_vertex))
        j=random.choice(list(min_vertex))
        if(has_path(g,i,j)) and i!=j and i<num_nodes and j< num_nodes:
                reach[[i],[j]]=1
                count=count+1
    e1=datetime.now()
    f.write("backbone!\n")
    res=get_result(reach)
    f.write(str(res))
    print("backbone\n")
    print(res)
    f.write("time backbone: "+str(e1-s1))
    f.write("\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss))  # in bytes 
    f.write("\n")
    
    s2=datetime.now()
    bfs_list=[] 
    if not backbone:
        for i in range(num_traverse):
            bfs_depth=sorted(list(nx.bfs_tree(g, source=sources[i], depth_limit=3).edges()))
            bfs_list.append(bfs_depth)
    else: 
        for i in range(num_traverse):
            bfs_depth=sorted(list(nx.bfs_tree(g, source=list(min_vertex)[i], depth_limit=3).edges()))
            bfs_list.append(bfs_depth)

    merged = list(itertools.chain(*bfs_list))

    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    e2=datetime.now()
    print("hop\n")
    res=get_result(reach)
    print(res)
    f.write("hop!\n")
    f.write(str(res))
    f.write("\n")
    f.write("time hop: "+str(e2-s2)+"\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss)) 
    f.write("\n")
    
    bfs_list=[] 
    s3=datetime.now()
    for i in range(num_traverse):
        bfs_depth=sorted(list(nx.bfs_tree(g, source=list(min_vertex)[i], depth_limit=3).edges()))
        bfs_list.append(bfs_depth)

    merged = list(itertools.chain(*bfs_list))

    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    e3=datetime.now()
    f.write("\n")
    f.write("time hopbackbone: "+str(e3-s3)+"\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss))  
    f.write("\n")
    print("hopbackbone\n")
    res=get_result(reach)
    print(res)
    f.write("hopbackbone!\n")
    f.write(str(res))
    f.write("\n")
    
    s4=datetime.now()
    dfs_list=[]
    if not backbone:
        for i in range(num_traverse):
            dfs_depth=sorted(list(nx.dfs_tree(g, source=sources[i], depth_limit=3).edges()))
            dfs_list.append(dfs_depth)
    else: 
        for i in range(num_traverse):
            dfs_depth=sorted(list(nx.dfs_tree(g, source=list(min_vertex)[i], depth_limit=3).edges()))
            dfs_list.append(dfs_depth)

    merged = list(itertools.chain(*dfs_list))

    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    e4=datetime.now()
    f.write("\n")
    f.write("time grail: "+str(e4-s4)+"\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss)) 
    f.write("\n")
    print("grail\n")
    res=get_result(reach)
    print(res)
    f.write("grail!\n")
    f.write(str(res))
    f.write("\n")
    
    s5=datetime.now()

    dfs_list=[]
    
    for i in range(num_traverse):
        dfs_depth=sorted(list(nx.dfs_tree(g, source=list(min_vertex)[i], depth_limit=3).edges()))
        dfs_list.append(dfs_depth)

    merged = list(itertools.chain(*dfs_list))

    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    
    e5=datetime.now()
    f.write("\n")
    f.write("time grailbackbone: "+str(e5-s5)+"\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss))  
    f.write("\n")
    print("grailbackbone\n")
    res=get_result(reach)
    print(res)
    f.write("grailbackbone!\n")
    f.write(str(res))
    f.write("\n")
    
    s6=datetime.now()
    object2=Graph_Sampling.SRW_RWF_ISRW()
    sample2= object2.random_walk_sampling_with_fly_back(g,query_budget,jump_prob)
    random_jump=sample2.edges()
    merged=random_jump
    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    
    e6=datetime.now()
    f.write("\n")
    f.write("time random jump: "+str(e6-s6)+"\n")
    print("random jump\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss)) 
    f.write("\n")
    res=get_result(reach)
    f.write("random jump\n")
    f.write(str(res))
    f.write("\n")
    print(res)
    object4=Graph_Sampling.ForestFire()
    sample4 = object4.forestfire(g,query_budget) 
    FF=sample4.edges()
    merged=FF
    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    s7=datetime.now()  
    object4=Graph_Sampling.ForestFire()
    sample4 = object4.forestfire(g,query_budget) 
    FF=sample4.edges()
    merged=FF
    reach=sps.lil_matrix((num_nodes, num_nodes), dtype=np.int8)
    count=0
    for (i,j) in merged:
        if count<query_budget:
            if i< num_nodes and j < num_nodes:
                reach[[i],[j]]=1
                count=count+1
            else:
                break
    e7=datetime.now()
    f.write("\n")
    f.write("time FF: "+str(e7-s7)+"\n")
    process = psutil.Process(os.getpid())
    f.write(str(process.memory_info().rss)) 
    f.write("\n")
    print("FF\n")
    res=get_result(reach)
    f.write("FF\n")
    f.write(str(res))
    print(res)
    f.close()    