#!/usr/bin/env python3

import time
from collections import defaultdict
import sys
import heapq

class WeightedDirectedGraph():
    def __init__(self, nr_vs):
        self.nr_vs = nr_vs #number of vertices
        self.es = defaultdict(set) #dictionary of sets that, for every vertex, stores heads of edges leaving it
        self.ws = dict() #dictionary that stores weight of each edge

    def show_nr_vs(self):
        return self.nr_vs

    def add_edge(self, t, h, w):
        (self.es[t]).add(h)
        self.ws[(t,h)] = w

    def get_edges(self, t):
        return self.es[t]
    
    def get_weight(self,t,h):
        return self.ws[(t,h)]

    def removeEdge(self, t, h):
        (self.es[t]).remove(h)
        del self.ws[(t,h)]

    def display(self):
        nr_es_present = 0
        for t in self.es:
            for h in self.es[t]:
                print("edge from %s to %s with weight %s" % (t, h, self.ws[(t,h)]))
                nr_es_present += 1
        print("there are %s edges in total" % nr_es_present)

    def run_Dijkstra(self,s,f):
        #This implementation of Dijkstra will compute the shortest distance from s to all other vertices, and the shortest path to f.
        self.ps = [None for v in range(self.nr_vs)] #list of parent of each vertex
        self.ds = [sys.maxsize for v in range(self.nr_vs)] #list of distances of each vertex
        q = [] #priority queue for Dijkstra
        
        heapq.heappush(q,[0,s]) #distance and vertex are paired, queue will sort by distance
        self.ds[s] = 0

        while q:
            dv = heapq.heappop(q)
            v = dv[1]
        
            for h in self.es[v]:
                new_d = self.ds[v] + self.ws[(v,h)]
                if new_d < self.ds[h]:
                    self.ds[h] = new_d
                    self.ps[h] = v
                    heapq.heappush(q,[new_d,h])
    
        if (f == None): #not interested in shortest path to a final destination
            return None, self.ds
    
        path = []
        end = f
        while end is not None:
            path.append(end)
            end = self.ps[end]

        path.reverse()

        return path, self.ds

    def run_Bellman_Ford(self,s):
        self.A = []
        self.A.append([sys.maxsize for v in range(self.nr_vs)])
        self.A[0][s] = 0

        k = 1
        while k<self.nr_vs:
            self.A.append(self.A[k-1][:])
            for w in range(self.nr_vs):
                for v in self.es[w]:
                    self.A[k][v] = min(self.A[k][v],self.A[k-1][w]+self.ws[(w,v)])
            if self.A[k][:] == self.A[k-1][:]:
                break
            k = k+1
        else:
            return None
    
        return self.A[-1][:]

    def run_Johnson(self):
        #extend graph by an extra (imaginary) vertex that is connected with distance zero to all other vertices
        self.nr_vs += 1
        for v in range(self.nr_vs-1):
            self.es[self.nr_vs-1].add(v)
            self.ws[(self.nr_vs-1,v)] = 0
        
        #calculate vertex-weight for reweighting
        self.vws = self.run_Bellman_Ford(self.nr_vs-1)
        
        if (self.vws == None):
            print('Bellman-Ford algorithm has detected a negative cycle')
            quit()
        
        #remove imaginary vertex
        del self.es[self.nr_vs-1]
        for v in range(self.nr_vs-1):
            del self.ws[(self.nr_vs-1,v)]
        self.nr_vs -= 1
        
        #reweighting edges
        for (u,v) in self.ws:
            self.ws[(u,v)] += self.vws[u] - self.vws[v]
    
        #calculate distances
        self.D = [[None for t in range(self.nr_vs)] for h in range(self.nr_vs)]
        edited = 0
        for u in range(self.nr_vs):
            _, self.D[u][:] = self.run_Dijkstra(u,None)

        #undo reweighting in distances
        for u in range(self.nr_vs):
            for v in range(self.nr_vs):
                self.D[u][v] += self.vws[v] - self.vws[u]

        return self.D

if __name__ == "__main__":
    file_name =  'Johnsontest3.txt'

    start_time = time.time()

    with open(file_name, 'r') as f:
        nr_vs, _ = f.readline().strip().split()
        nr_vs = int(nr_vs)
        graph = WeightedDirectedGraph(nr_vs)

        for line in f:
            t, h, w = line.strip().split()
            t, h, w = int(t)-1, int(h)-1, int(w)
            graph.add_edge(t,h,w)

    allD = graph.run_Johnson()
    print(min(min(allD)))
    
    end_time = time.time()
    print(end_time - start_time)


