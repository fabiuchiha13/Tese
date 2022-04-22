import networkx as nx
import pandas as pd
import numpy as np
from cdlib import algorithms, readwrite, evaluation, NodeClustering
from cdlib.benchmark import LFR
import matplotlib.pyplot as plt
import time
import csv
from infomap import Infomap
from pathlib import Path
import pickle
import os

def averageDegree(graph):
    degrees = [val for (node, val) in graph.degree()]
    sum = 0
    for d in degrees:
        sum += d
    return sum/len(degrees)

def kin(graph):
    nodes = graph.nodes()
    sum = 0
    for n in nodes:
        sum += graph.in_degree(n)
    return sum/len(nodes)

def kout(graph):
    nodes = graph.nodes()
    sum = 0
    for n in nodes:
        sum += graph.out_degree(n)
    return sum/len(nodes)

def APL(graph):
    if (nx.is_directed(graph)):
        largestComponent = 0
        largestAPL = 0
        for C in (graph.subgraph(c) for c in nx.weakly_connected_components(graph)):
            if largestComponent<len(C.nodes):
                largestComponent = len(C.nodes)
                apl = nx.average_shortest_path_length(C)
                largestAPL = apl
        return largestAPL
    else:
        for C in (graph.subgraph(c) for c in nx.connected_components(graph)):
            return nx.average_shortest_path_length(C)

def transcriptionFactors(graph):
    sourceNodes = []
    for e in graph.edges():
        sourceNodes.append(e[0])
    return len(set(sourceNodes))

def targetGenes(graph):
    targetNodes = []
    for e in graph.edges():
        targetNodes.append(e[1])
    return len(set(targetNodes))

def networkInfo(graph):
    print("Transcription Factors:", transcriptionFactors(graph))
    print("Target Genes:", targetGenes(graph))
    print("Average degree:", averageDegree(graph))
    print("Internal average degree:", kin(graph))
    print("External average degree:", kout(graph))
    print("Clustering coefficient:", nx.average_clustering(graph))
    print("Average Path Length (highest value):", APL(graph))

def plotDegreeDistribution(graph):
    fig = plt.figure(figsize=(6*1.61803398875, 6))
    ax = plt.axes((0.2, 0.2, 0.70, 0.70), facecolor='w')
    d = np.array(nx.degree_histogram(graph))
    y = d / len(graph.nodes)
    x = np.arange(len(y))
    ax.plot(x,y,"go")
    ax.set_xlabel("k")
    ax.set_ylabel("Pk")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title("Degree distribution")
    #ax.legend()
    fig.savefig(("Images/DegreeDistribution_%s.png" % (graph.name)))
    plt.close(fig)

def txtToNetworkx(fileName, directed):
    raw = pd.read_csv(fileName, sep="\t", usecols=[0,2], header=None)
    raw.columns = ['Source', 'Target']
    if directed:
        network = nx.from_pandas_edgelist(raw, source='Source', target='Target', edge_attr=None, create_using=nx.DiGraph())
    else:
        network = nx.from_pandas_edgelist(raw, source='Source', target='Target', edge_attr=None)
    network.name = fileName.split(".")[0]
    return network

def dataToJSON(data, filename):
    with open(filename, 'wb') as object_file:
        pickle.dump(data, object_file)

def JSONtoData(fileName):
    with open(fileName, 'rb') as object_file:
        data = pickle.load(object_file)
        return data


def infomap(graph):
    infomapWrapper = Infomap("--two-level --silent")
    for e in graph.edges():
        infomapWrapper.addLink(int(e[0]), int(e[1]))
    infomapWrapper.run();
    communities = [[]] * infomapWrapper.num_top_modules
    for node in infomapWrapper.tree:
        if node.is_leaf:
            if communities[node.module_id-1]:
                communities[node.module_id-1].append(node.node_id)
            else:
                communities[node.module_id-1] = [node.node_id]
    return NodeClustering(communities, graph, method_name="infomap")

def goorflistToDicts():
    with open('Ontology/goorflist.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter =' ')
        next(csvreader)
        
        ca_goorflist = {}
        cg_goorflist = {}
        cp_goorflist = {}
        ct_goorflist = {}
        zb_goorflist = {}
        sc_goorflist = {}
        kl_goorflist = {}
        yl_goorflist = {}
        km_goorflist = {}
        kp_goorflist = {}
        
        for row in csvreader:
            goid = row[0].split("\t")[0]
            goorf = row[0].split("\t")[1]
            
            if 'albicans' in row[1]:
                if goorf in ca_goorflist:
                    auxList = ca_goorflist[goorf]
                    auxList.append(goid)
                    ca_goorflist[goorf] = auxList
                else:
                    ca_goorflist[goorf] = [goid]
            else:
                if 'glabrata' in row[1]:
                    if goorf in cg_goorflist:
                        auxList = cg_goorflist[goorf]
                        auxList.append(goid)
                        cg_goorflist[goorf] = auxList
                    else:
                        cg_goorflist[goorf] = [goid]
                else: 
                    if 'parapsilosis' in row[1]:
                        if goorf in cp_goorflist:
                            auxList = cp_goorflist[goorf]
                            auxList.append(goid)
                            cp_goorflist[goorf] = auxList
                        else:
                            cp_goorflist[goorf] = [goid]
                    else: 
                        if 'tropicalis' in row[1]:
                            if goorf in ct_goorflist:
                                auxList = ct_goorflist[goorf]
                                auxList.append(goid)
                                ct_goorflist[goorf] = auxList
                            else:
                                ct_goorflist[goorf] = [goid]
                        else: 
                            if 'bailii' in row[1]:
                                if goorf in zb_goorflist:
                                    auxList = zb_goorflist[goorf]
                                    auxList.append(goid)
                                    zb_goorflist[goorf] = auxList
                                else:
                                    zb_goorflist[goorf] = [goid]
                            else: 
                                if 'cerevisiae' in row[1]:
                                    if goorf in sc_goorflist:
                                        auxList = sc_goorflist[goorf]
                                        auxList.append(goid)
                                        sc_goorflist[goorf] = auxList
                                    else:
                                        sc_goorflist[goorf] = [goid]
                                else: 
                                    if 'lactis' in row[1]:
                                        if goorf in kl_goorflist:
                                            auxList = kl_goorflist[goorf]
                                            auxList.append(goid)
                                            kl_goorflist[goorf] = auxList
                                        else:
                                            kl_goorflist[goorf] = [goid]
                                    else: 
                                        if 'lipolytica' in row[1]:
                                            if goorf in yl_goorflist:
                                                auxList = yl_goorflist[goorf]
                                                auxList.append(goid)
                                                yl_goorflist[goorf] = auxList
                                            else:
                                                yl_goorflist[goorf] = [goid]
                                        else: 
                                            if 'marxianus' in row[1]:
                                                if goorf in km_goorflist:
                                                    auxList = km_goorflist[goorf]
                                                    auxList.append(goid)
                                                    km_goorflist[goorf] = auxList
                                                else:
                                                    km_goorflist[goorf] = [goid]
                                            else: 
                                                if 'phaffii' in row[1]:
                                                    if goorf in kp_goorflist:
                                                        auxList = kp_goorflist[goorf]
                                                        auxList.append(goid)
                                                        kp_goorflist[goorf] = auxList
                                                    else:
                                                        kp_goorflist[goorf] = [goid]
        dataToJSON(ca_goorflist, "Ontology/goorflist_Ca")
        dataToJSON(cg_goorflist, "Ontology/goorflist_Cg")
        dataToJSON(cp_goorflist, "Ontology/goorflist_Cp")
        dataToJSON(ct_goorflist, "Ontology/goorflist_Ct")
        dataToJSON(zb_goorflist, "Ontology/goorflist_Zb")
        dataToJSON(sc_goorflist, "Ontology/goorflist_Sc")
        dataToJSON(kl_goorflist, "Ontology/goorflist_Kl")
        dataToJSON(yl_goorflist, "Ontology/goorflist_Yl")
        dataToJSON(km_goorflist, "Ontology/goorflist_Km")
        dataToJSON(kp_goorflist, "Ontology/goorflist_Kp")

def gotermsToDicts():
    with open('Ontology/goterms.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter ='\t')
        next(csvreader)
        
        gotermsDict = {}
        for row in csvreader:
            terms = row
            goid = terms[0]
            terms.remove(goid)
            gotermsDict[goid] = terms
        dataToJSON(gotermsDict, "Ontology/gotermsDict")
    
def goparentsToDicts():
    with open('Ontology/goparents.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter ='\t')
        next(csvreader)
        i = 0
        goparentsDict = {}
        gochildsDict = {}
        for row in csvreader:
            i += 1
            parent = row[0]
            child = row[1]
            if parent in goparentsDict:
                auxList = goparentsDict[parent]
                auxList.append(child)
                goparentsDict[parent] = auxList
            else:
                goparentsDict[parent] = [child]
            if child in gochildsDict:
                auxList = gochildsDict[child]
                auxList.append(parent)
                gochildsDict[child] = auxList
            else:
                gochildsDict[child] = [parent]
        dataToJSON(goparentsDict, "Ontology/goparentsDict")
        dataToJSON(gochildsDict, "Ontology/gochildsDict")
        