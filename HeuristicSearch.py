import bisect
from collections import defaultdict, deque
import re
import typing
import numpy
import matplotlib.pyplot as plt
import pandas
import tkinter
import turtle
import copy
import itertools
import cv2
import abc
import string
import time
import cv2

class DSU:
    def __init__(self):
        self.dsu = []
        self.sizes = []
    def find(self,x):
        if(dsu[x]==x):
            return x
        return dsu[x] = find(dsu[x])
    def make_set(self,n):
        for i in range(n):
            self.dsu[i] = i
            self.sizes[i] = u
    def join(self,a,b):
        a = find(a);
        b = find(b);
    
        if (a != b):
            if (sizes[a] < sizes[b]) :
                dsu[a] = b;
                sizes[b] += sizes[a];
            else:
                dsu[b] = a;
                sizes[a] += sizes[b];

class HeuristicSearch:
    #board will be represented as a 2D matrix of colors
    def __init__(self, board):
        self.board=board



