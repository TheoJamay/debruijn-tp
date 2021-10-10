#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Théo Jamay"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Théo Jamay"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Théo Jamay"
__email__ = "jamay.theo@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open (fastq_file, "r") as fq_file : 
        text = fq_file.readlines()
        for i in range(1,len(text),4) :
            yield text[i].strip()


def cut_kmer(read, kmer_size):
    n=len(read)
    for i in range(0,n-(kmer_size-1)):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for read in read_fastq(fastq_file) :
        for kmer in cut_kmer(read, kmer_size) :
            if kmer not in kmer_dict :
                kmer_dict[kmer] = 1
            else :
                kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    G = nx.DiGraph()
    #G.add_nodes_from(kmer_dict)
    for key in kmer_dict :
        G.add_edge(key[:-1], key[1:],weight=kmer_dict[key])
    #G.add_edge(value1, value2, weight=poids)
    #nx.add_edge(value1, value2 [, weight] )
    return(G)


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list :
        graph.remove_nodes_from(path)
    return graph


def std(data):
    return round(statistics.stdev(data),1)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    if std(weight_avg_list) > 0 :
        i_max = weight_avg_list.index(max(weight_avg_list))
        i = 0
        for path in path_list :
            if i != i_max :
                if delete_entry_node == False :
                    path = path[1:]
                if delete_sink_node == False :
                    path = path[:-1]
                graph.remove_nodes_from(path)
            i += 1
    elif std(path_length) > 0 :
        i_max = path_length.index(max(path_length))
        i = 0
        for path in path_list :
            if i != i_max :
                if delete_entry_node == False :
                    path = path[1:]
                if delete_sink_node == False :
                    path = path[:-1]
                graph.remove_nodes_from(path)
            i += 1
    else : 
        i_max = random.randint(0,len(weight_avg_list))
        i = 0
        for path in path_list :
            if i != i_max :
                if delete_entry_node == False :
                    path = path[1:]
                if delete_sink_node == False :
                    path = path[:-1]
                graph.remove_nodes_from(path)
            i += 1
    return graph
    

def path_average_weight(graph, path):
    #graph.subgraph(path).edges(data=True)
    weight = statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
    return weight


def solve_bubble(graph, ancestor_node, descendant_node):
    weight_avg_list = []
    path_length = []
    path_list = nx.all_simple_paths(graph,source = ancestor_node,target = descendant_node)
    for p in path_list :
        weight = path_average_weight(graph, p)
        weight_avg_list.append(weight)
        path_length.append(len(p))
    path_list = nx.all_simple_paths(graph,source = ancestor_node,target = descendant_node)    
    graphii = select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False)
    return graphii
    
    

def simplify_bubbles(graph):
    bubble = False
    for node in graph :
        if node in graph :
            pred_list = []
            pred = graph.predecessors(node)
            for p in pred :
                pred_list.append(p)
            if len(pred_list) > 1 :
                for i in range(0,len(pred_list)) :
                    for j in range(i+1,len(pred_list)) :
                        noeud_ancêtre = nx.lowest_common_ancestor(graph, pred_list[i], pred_list[j])
                        if noeud_ancêtre != None :
                            bubble = True
                            break
                    if bubble:                
                        graph = simplify_bubbles(solve_bubble(graph, noeud_ancêtre, node))


def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    #print(list(graph.nodes))
    in_nodes = []
    for node in graph :
        pred = graph.predecessors(node)
        #print("check : {}".format(node))
        #print(len(pred))
        j = 0
        for i in pred :
            j += 1
            #print(i,end="")
        #print("\n")
        if j == 0 :
            in_nodes.append(node)
        #print(graph.predecessors(node))
    #print(graph.nodes[0])
    #graph.DiGraph.predecessors(n)
    return(in_nodes)

def get_sink_nodes(graph):
    out_nodes = []
    for node in graph :
        suc = graph.successors(node)
        j = 0
        for i in suc :
            j += 1
        if j == 0 :
            out_nodes.append(node)
    return(out_nodes)

def get_contigs(graph, starting_nodes, ending_nodes):
    list_contig = []
    for i in range(0,len(starting_nodes)) :
        for j in range(0,len(ending_nodes)) :
            if nx.has_path(graph, starting_nodes[i], ending_nodes[j]) == True :
                contig = nx.all_simple_paths(graph, starting_nodes[i], ending_nodes[j])
                for c in contig :
                    seq = c[0]
                    for nucl in range(1,len(c)) :
                        seq += c[nucl][-1]
                t = (seq, len(seq))
                list_contig.append(t)
    return(list_contig)


def save_contigs(contigs_list, output_file):
    with open(output_file,"w") as out_file :
        for cont in range(0,len(contigs_list)) :
            out_file.write(">contig_{} len={}\n".format(cont, contigs_list[cont][1]))
            text = fill(contigs_list[cont][0])
            out_file.write("{}\n".format(text))
        

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # for i in read :
    #     print(i)
    #kmer = cut_kmer(read,args.kmer_size)
    #print(kmer)
    # for j in kmer :
    #     print(j)
    
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    #starting_nodes = get_starting_nodes(graph)
    #ending_nodes = get_sink_nodes(graph)
    #contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    #output_file = "contig.fasta"
    #save_contigs(contigs_list, output_file)
    # for u,v in graph.edges() :
    #     allpath = nx.all_simple_paths(graph,u,v)
    #     for p in allpath :
    #         if len(p) > 2 :
    #             print(p)
    # for i in starting_nodes :
    #     for j in ending_nodes :
    #         path_list = nx.all_simple_paths(graph,source = i,target = j)
    #         for p in path_list :
    #             path_average_weight(graph, p)
    #         remove_paths(graph, path_list, True, True)
    simplify_bubbles(graph)
    #nx.all_simple_paths.statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

    #path_average_weight(graph, path)
    # for key in kmer_dict :
    #     print(key,kmer_dict[key])
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
