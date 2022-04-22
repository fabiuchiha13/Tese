{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "4c85d15b",
   "metadata": {},
   "outputs": [],
   "source": [
    "import networkx as nx\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "from cdlib import algorithms, readwrite, evaluation, NodeClustering\n",
    "from cdlib.benchmark import LFR\n",
    "import matplotlib.pyplot as plt\n",
    "import time\n",
    "import csv\n",
    "from infomap import Infomap"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "19b500ee",
   "metadata": {},
   "outputs": [],
   "source": [
    "def averageDegree(graph):\n",
    "    degrees = [val for (node, val) in graph.degree()]\n",
    "    sum = 0\n",
    "    for d in degrees:\n",
    "        sum += d\n",
    "    return sum/len(degrees)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "723d7c7f",
   "metadata": {},
   "outputs": [],
   "source": [
    "def kin(graph):\n",
    "    nodes = graph.nodes()\n",
    "    sum = 0\n",
    "    for n in nodes:\n",
    "        sum += graph.in_degree(n)\n",
    "    return sum/len(nodes)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "4abf07b1",
   "metadata": {},
   "outputs": [],
   "source": [
    "def kout(graph):\n",
    "    nodes = graph.nodes()\n",
    "    sum = 0\n",
    "    for n in nodes:\n",
    "        sum += graph.out_degree(n)\n",
    "    return sum/len(nodes)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "114bf684",
   "metadata": {},
   "outputs": [],
   "source": [
    "def APL(graph):\n",
    "    if (nx.is_directed(graph)):\n",
    "        largestComponent = 0\n",
    "        largestAPL = 0\n",
    "        for C in (graph.subgraph(c) for c in nx.weakly_connected_components(graph)):\n",
    "            if largestComponent<len(C.nodes):\n",
    "                largestComponent = len(C.nodes)\n",
    "                apl = nx.average_shortest_path_length(C)\n",
    "                largestAPL = apl\n",
    "        return largestAPL\n",
    "    else:\n",
    "        for C in (graph.subgraph(c) for c in nx.connected_components(graph)):\n",
    "            return nx.average_shortest_path_length(C)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "068e7519",
   "metadata": {},
   "outputs": [],
   "source": [
    "def transcriptionFactors(graph):\n",
    "    sourceNodes = []\n",
    "    for e in graph.edges():\n",
    "        sourceNodes.append(e[0])\n",
    "    return len(set(sourceNodes))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "d4836caa",
   "metadata": {},
   "outputs": [],
   "source": [
    "def targetGenes(graph):\n",
    "    targetNodes = []\n",
    "    for e in graph.edges():\n",
    "        targetNodes.append(e[1])\n",
    "    return len(set(targetNodes))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "11bb4930",
   "metadata": {},
   "outputs": [],
   "source": [
    "def networkInfo(graph):\n",
    "    print(\"Transcription Factors:\", transcriptionFactors(graph))\n",
    "    print(\"Target Genes:\", targetGenes(graph))\n",
    "    print(\"Average degree:\", averageDegree(graph))\n",
    "    print(\"Internal average degree:\", kin(graph))\n",
    "    print(\"External average degree:\", kout(graph))\n",
    "    print(\"Clustering coefficient:\", nx.average_clustering(graph))\n",
    "    print(\"Average Path Length (highest value):\", APL(graph))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "5232b60a",
   "metadata": {},
   "outputs": [],
   "source": [
    "def plotDegreeDistribution(graph):\n",
    "    fig = plt.figure(figsize=(6*1.61803398875, 6))\n",
    "    ax = plt.axes((0.2, 0.2, 0.70, 0.70), facecolor='w')\n",
    "    d = np.array(nx.degree_histogram(graph))\n",
    "    y = d / len(graph.nodes)\n",
    "    x = np.arange(len(y))\n",
    "    ax.plot(x,y,\"go\")\n",
    "    ax.set_xlabel(\"k\")\n",
    "    ax.set_ylabel(\"Pk\")\n",
    "    ax.set_yscale('log')\n",
    "    ax.set_xscale('log')\n",
    "    ax.set_title(\"Degree distribution\")\n",
    "    #ax.legend()\n",
    "    fig.savefig((\"Images/DegreeDistribution_%s.png\" % (graph.name)))\n",
    "    plt.close(fig)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "21835d04",
   "metadata": {},
   "outputs": [],
   "source": [
    "def txtToNetworkx(fileName, directed):\n",
    "    raw = pd.read_csv(fileName, sep=\"\\t\", usecols=[0,2], header=None)\n",
    "    raw.columns = ['Source', 'Target']\n",
    "    if directed:\n",
    "        network = nx.from_pandas_edgelist(raw, source='Source', target='Target', edge_attr=None, create_using=nx.DiGraph())\n",
    "    else:\n",
    "        network = nx.from_pandas_edgelist(raw, source='Source', target='Target', edge_attr=None)\n",
    "    network.name = fileName.split(\".\")[0]\n",
    "    return network"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c3f98643",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.10"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
