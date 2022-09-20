import bisect
from collections import defaultdict, deque
import multiprocessing
from concurrent.futures import *
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
import random

'''
VIEW THE SOLVING OF THE FLOW FREE CSP AS BOTH A SAT PROBLEM AND A MAPF(MULTI-AGENT PATH FINDING PROBLEM)
https://www.youtube.com/watch?v=DIcRFQ2xzlA
https://cse442-17f.github.io/Conflict-Driven-Clause-Learning/
https://courses.cs.washington.edu/courses/cse507/17wi/lectures/L02.pdf
'''

class CNFType:
    def __init__(self, node, pos):
        self.node=node
        self.pos=pos
    def __hash__():
        pass
    def __probe__():
        pass

'''
ENCODES THE CONSTRAINTS FOR THE FLOW FREE PROBLEM INTO OUR ALGORITHM
'''
class SAT_Optimizations:
    #implement an incremental heuristic here
    def __init__(self):
        pass


class FlowFreeCNFConditions:
    CLAUSE_ID = 0
    #each variable in our CNF will be represented with a 4-tuple (x_pos,y_pos,color,boolean_value)
    def __init__(self, board, num_cols, col_list):
        self.board = board
        self.col_list = col_list
        (x,y) = board.shape
        self.CNF = defaultdict(list) #represented as a bipartite graph
        self.CNF_Matrix = np.zeroes( (x*y*num_cols) )
    def generate_endpoint_cells(self):
        (i,j) = board.shape
        endpoints = []
        non_endpoints = []
        for idx1 in range(i):
            for idx2 in range(j):
                if(board[idx1][idx2] == (0,0,0) ):
                    endpoints.append( (idx1,idx2) )
                else:
                    non_endpoints.append( (idx1,idx2) )
        self.endpoints = endpoints
        self.non_endpoints = non_endpoints
    def instantiate_case_one(self):
        for (x,y) in self.non_endpoints:
            for j in range(self.num_cols):
                self.CNF[CLAUSE_ID].append( (x,y,j,1) )
        CLAUSE_ID += 1
    def instantiate_case_two(self):
        for (x,y) in self.non_endpoints:
            for j in range(num_cols):
                for k in range(num_cols):
                    if(j==k):
                        continue
                    self.CNF[CLAUSE_ID].
            CLAUSE_ID += 1
    def instantiate_case_three():
        CLAUSE_ID += 1
        pass
    def instantiate_case_four():
        CLAUSE_ID += 1
        pass
    def instantiate_case_five():
        CLAUSE_ID += 1
        pass
    def instantitate_case_six():
        CLAUSE_ID += 1
        pass
    
                    
class CNFSolver(FlowFreeCNFConditions):
    #represent our CSP as a bipartite graph with edges from (literals/variables -> clauses)
    #for small testing, clause parameter is provided
    def __init__(self, board, board_x, board_y, num_cols, col_list):
        super(CNFSolver,self).__init__(board, num_cols)
        self.board_x = board_x
        self.board_y = board_y
        self.num_cols = num_cols
        self.deg = {}
        self.instantiate_case_one()
        self.instantiate_case_two()
        self.instantiate_case_three()
        self.instantiate_case_four()
        self.instantiate_case_five()
        self.instantiate_case_six()
    def __init__(self, ORIGINAL_CLAUSE):
        self.graph = defaultdict(list)
        self.deg = {}
        self.ORIGINAL_CLAUSE = ORIGINAL_CLAUSE
        self.num_vars = 0
        self.implication_graph = defaultdict(list)
        self.satisfied_clauses = []
        self.DPLL_OPS = 0
        self.CDCL_OPS = 0
        self.CDCL_Level = 1
        self.CDCL_Stack = []
    #implication graph function seems to be working properly (stress test more later)
    def construct_implication_graph(self, tmp_assignment):
        print("ORIGINAL CLAUSE: {}".format(self.ORIGINAL_CLAUSE) )
        print("PARTIAL ASSIGNMENT: {}".format(tmp_assignment))

        if(len(tmp_assignment)==0):
            return

        for i in self.ORIGINAL_CLAUSE:
            ACTUAL_CLAUSE = i[0]
            print("ACTUAL CLAUSE: {}".format(ACTUAL_CLAUSE))
            rep = {}
            for (x,y) in tmp_assignment:
                rep[x] = y
            unsatisfied = []
            satisfied = []
            IS_SATISFIED = False
            for j in ACTUAL_CLAUSE:
                (g,h) = j
                if(g not in rep.keys()):
                    satisfied.append( (g,h) )
                else:
                    if(rep[g]==h):
                        satisfied.append( (g,h) )
                        IS_SATISFIED = True
                    elif( (rep[g] != h) ):
                        unsatisfied.append( (g,h) )
            #if len(satisfied)==1 and IS_SATISFIED is false, it means that the current clause is a unit clause under the given partial assignment
            if((not IS_SATISFIED) and (len(satisfied)==1) and (i[1] not in self.satisfied_clauses) ):
                self.satisfied_clauses.append(i[1])
                for (f,s) in unsatisfied:
                    if( (f,s) not in self.implication_graph[unsatisfied[0].__hash__()]):
                        (uf1,uf2) = satisfied[0]
                        print("key: {}, {}, value: {}, {}".format(f,rep[f],uf1,uf2))
                        self.implication_graph[ (f,rep[f]).__hash__()].append( (satisfied[0],self.CDCL_Level) )
    def display_implication_graph(self):
        k = self.implication_graph.keys()
        for i in k:
            print("GRAPH: {}".format(self.implication_graph[i]))


    #find UIP in our implication graph
    def find_unique_implication_point(self):
        pass

    #Given a partial assignment, determine if any contradictions arise (will be called secondhand in our unit prop)
    def check_contradictions(self,assignment):
        d = {}
        for (v,t) in assignment:
            d[v]=t
        CLAUSE_INDEX=-1
        for i in self.ORIGINAL_CLAUSE:
            CONDITION_MET = False
            for (f,s) in i[0]:
                if(f not in d.keys()):
                    CONDITION_MET = True
                    break
                if(d[f] == s):
                    CONDITION_MET=True
                    break
            if(not CONDITION_MET):
                return (True,i[1])
        return (False,-1)

    #unit propagation(or boolean constraint propagation) to reduce our CNF
    #TODO: Implement UNSAT Finding with the unit_propagation() as well (DONE)
    def BOOLEAN_CONSTRAINT_PROPAGATION(self, reduc, CDCL_Mode=False):
        LIM_COUNTER = 5
        CUR_COUNTER = 0
        assignment = []
        print("BOOLEAN CONSTRAINT PROPAGATION")
        print(len(reduc))
        while(True):
            CUR_COUNTER += 1
            n_reduc = []
            added = {}
            truth = {}
            for i in reduc:
                print(i)
                if(len(i[0])==1):
                    (x,y) = i[0][0]
                    added[x] = True
                    assignment.append( (x,y) )
                    truth[x] = y
            for i in reduc:
                tmp_one = i[0][:]
                CLAUSE_INDEX = -1
                if(CDCL_Mode):
                    CLAUSE_INDEX = i[1]
                add = True
                for j in range(len(i[0])):
                    (f,s) = i[0][j]
                    #print("{},{}".format(f,s))
                    #make this O(log N) with a bsearch
                    if(f in added.keys()):
                        if(s == truth[f]):
                            add = False
                            #print("FALSE")
                        else:
                            tmp_one = tmp_one[0:j] + tmp_one[j+1:]
                if(add):
                    if(not CDCL_Mode):
                        n_reduc.append(tmp_one)
                    else:
                        n_reduc.append( [tmp_one,CLAUSE_INDEX] )
                  
            if(len(added)==0):
                break
            reduc = n_reduc[:]
        print("END")
        
        if(CDCL_Mode):
            (truth, conflicting_clause) = self.check_contradictions(assignment)
            if(truth == True):   
                return (-1,conflicting_clause)
            return (reduc, assignment)  
        else:
            return (reduc,assignment)
   #reduction with pure literal elimination built-in.
    '''
    what we'll do is augment our CNF data structure to maintain the original clause which each subclause in the reduced CNF
    came from, and then use that when constructing our implication graph
    '''
    def SAT_REDUCTION(self, reduc, var, truth_value, partial_assignment, CDCL_Mode = False):
        print("CURRENT CLAUSE: {}".format(reduc))
        print("CURRENT ASSIGNMENT: {}".format(partial_assignment))
        print("SAT REDUCTION")
        polarity = defaultdict(set)
        for i in reduc:
            for j in i[0]:
                (variable, boolean) = j
                polarity[variable].add(boolean)
        ret = []
        UNSAT = False
        for i in reduc:
            tmp = []
            satisfied = False
            CLAUSE_INDEX = -1
            if(CDCL_Mode):
                CLAUSE_INDEX = i[1]
            for j in i[0]:
                (x,y) = j
                if( (x==var and y==truth_value) ):
                    satisfied = True
                elif(len(polarity[x])==1):
                    tmp.append( (x,y) )
                    if((x==var) and (len(polarity[x])==1) and (y != truth_value) ):
                      UNSAT = True
                elif(x == var and y != truth_value):
                    continue
                else:
                    tmp.append(j)
            if(satisfied == False):
                if(not CDCL_Mode):
                    ret.append(tmp)
                else:
                    ret.append( [tmp,CLAUSE_INDEX] )
            if(CDCL_Mode and satisfied == False):#creating implication graph DAG
                made_guesses = []
                remain_var = 0
                truth_var = 0
                cnt = 0
                for (v,t) in tmp:
                    if( not ( (v,t) in partial_assignment) ):
                        remain_var = v
                        truth_var = t
                        cnt += 1
                    else:
                        made_guesses.append( (v,t) )
        print("END")      
        if(UNSAT):
            return [-1]
        return ret
    #DPLL Heuristic(Precursor to the Conflict-Based Clause Learning Algorithm) (COMPLETE)
    def DPLL_Heuristic(self, reduc, assignment):
        print("CURRENT CLAUSE: {}".format(reduc))
        print("CURRENT ASSIGNMENT: {}".format(assignment))
        print('\n')
        self.DPLL_OPS += 1
        undetermined_vars = set()
        (rex, acc) = self.BOOLEAN_CONSTRAINT_PROPAGATION(reduc,True)
        reduc = rex[:]
        assignment.extend(acc)
        if(len(reduc)==0):
            return assignment
        #print("{},{}".format(reduc,assignment))
        for i in reduc:
            for j in i[0]:
                (x,y) = j
                undetermined_vars.add(x)
        self.num_vars = max(self.num_vars,len(undetermined_vars))
        ans = []
        a1 = assignment[:]
        a2 = assignment[:]
        print(undetermined_vars)
        for var in undetermined_vars:
            if(len(ans)==0):
                ret1 = self.SAT_REDUCTION(reduc,var,1,assignment,True)
                if(ret1[0] == -1):
                    continue
                #print( (var,1) )
                a1.append( (var, 1) )
                #print(a1)
                ans = self.DPLL_Heuristic(ret1,a1)
                break
            else:
                ret2 = self.SAT_REDUCTION(reduc,var,0,assignment,True)
                if(ret2[0] == -1):
                    continue
                a2.append( (var,0) )
                ans = self.DPLL_Heuristic(ret2,a2)
                break
                
        return ans
    def ADD_ARBITARY_VARIABLES(self, assignment):
        print("TOTAL VARIABLE COUNT: {}".format(self.num_vars))
        d = [-1 for i in range(1,self.num_vars+2,1)]
        for i in assignment:
            (f,s) = i
            d[f] = 1
        for i in range(len(d)):
            if(d[i] == -1):
                assignment.append( (i, (0,1) ) )
    def analyze_conflict(self, CNF):
        pass
    def verify(clause, assignment):
        pass

    #make this method iterative for nonchronological backtracking and not recursive (unlike my DPLL implementation)
    '''
    Conflict Driven Clause Learning uses two new heuristics from the DPLL Implementation:

    1. Clause Learning
    2a. Nonchronological Backtracking
    2b. Implication Graph

    '''

    def CDCL_Heuristic(self, reduc):
        assignment = []
        LIM = 100
        CNT = 0
        while(CNT < LIM):
            CNT += 1
            self.CDCL_OPS += 1
            self.CDCL_Level += 1
            print('\n')
            self.construct_implication_graph(assignment)
            self.display_implication_graph()
            undetermined_vars = set()
            #remove field giving us knowledge of where the clause originated from
            (rex, acc) = self.BOOLEAN_CONSTRAINT_PROPAGATION(reduc,True)
            if(rex==-1):
                print("CONFLICT REACHED!")
                #here, we implement our clause learning, add it to reduc, and use nonchronological backtracking heuristic
                #self.display_implication_graph()
                pass
            else:
                reduc = rex[:]
                assignment.extend(acc)
            if(len(reduc)==0):
                return assignment
            for i in reduc:
                for j in i[0]:
                    (x,y) = j
                    undetermined_vars.add(x)
            self.num_vars = max(self.num_vars,len(undetermined_vars))
            ans = []
            a1 = assignment[:]
            a2 = assignment[:]
            #improve this deciding heuristic
            for var in undetermined_vars:
                sav2 = self.SAT_REDUCTION(reduc,var,1,assignment,True)
                sav1 = self.SAT_REDUCTION(reduc,var,0,assignment,True)
                if(sav2 == [-1]):
                    print("ONE")
                    reduc = sav1[:]
                    #print( (var,1) )
                    #print("REDUCTION ONE: {}".format(sav1))
                    a1.append( (var, 1) )
                    #print(a1)                
                    assignment = a1
                    #ans = self.CDCL_Heuristic(ret1,a1, guess_stack)
                    break
                elif(sav1 == [-1]):
                    print("TWO")
                    #print("REDUCTION TWO: {}".format(sav2) )
                    reduc = sav2[:]
                    a2.append( (var,0) )
                    assignment = a2
                    break
        return assignment
           

def CONSTRUCT_TEST_CASE(MAX_VAR, NUM_CLAUSES):
    bin_rep = []
    case = []
    for i in range( (2**MAX_VAR) + 4):
        bin_rep.append( bin(i)[2:] )
    for i in range(NUM_CLAUSES):
        r = random.randint(0,2**MAX_VAR)
        s = str(bin_rep[r])
        CNF_CLAUSE = []
        for j in range(len(s)):
            if(ord(s[j])==49):
                CNF_CLAUSE.append( (j,1) )
            else:
                CNF_CLAUSE.append( (j,0) )
        case.append( (CNF_CLAUSE, i) )
    return case



    
# 0 = FALSE, 1 = TRUE
def TEST_SET_NONINCREMENTAL():
    PLS_DPLL = [ [(2,0),(3,0),(4,0),(5,1)],[(1,0),(5,0),(6,1)],[(5,0),(7,1)],[(1,0),(6,0),(7,0)],[(1,0),(2,0),(5,1)],[(1,0),(3,0),(5,1)     ],[(1,0),(4,0),(5,1)],[(1,0),(2,1),(3,1),(4,1),(5,1),(6,0)] ]
    
    PLS_CDCL =  [ [ [(2,0),(3,0),(4,0),(5,1)],0],[ [(1,0),(5,0),(6,1)], 1], [ [(5,0),(7,1)], 2],[ [(1,0),(6,0),(7,0)],3],[[(1,0),(2,0),(5,1)],4],[ [(1,0),(3,0),(5,1)],5],[[(1,0),(4,0),(5,1)],6],[[(1,0),(2,1),(3,1),(4,1),(5,1),(6,0)],7] ]
    
    PLS_CDCL_TWO = [ [[(1,0),(2,1),(4,0)], 0], [ [(1,0), (2,0),(3,1)], 1], [ [(3,0),(4,0)],2], [[(4,1),(5,1),(6,1)],3],[[(5,0),(7,1)],4], [[(6,1),(7,1),(8,0)],5] ]



    give_case = CONSTRUCT_TEST_CASE(4,4)
    cop_case = give_case[:]
    print("CASE: {}".format(give_case))
    c = CNFSolver(PLS_CDCL_TWO[:])
    c2 = CNFSolver(PLS_CDCL_TWO[:])
    
    '''
    OPTIMIZED_CNF_ASSIGNMENT = c.DPLL_Heuristic(PLS_CDCL,[])
    c.ADD_ARBITARY_VARIABLES(OPTIMIZED_CNF_ASSIGNMENT)
    print("ANSWER: {}".format(OPTIMIZED_CNF_ASSIGNMENT))
    print("DPLL OPERATIONS: {}".format(c.DPLL_OPS))
    print("ORIGINAL CNF: {}".format(PLS_DPLL))
    '''
    '''
    CDCL_CNF_ASSIGNMENT = c.CDCL_Heuristic(PLS_CDCL_TWO[:],[])
    c.ADD_ARBITARY_VARIABLES(CDCL_CNF_ASSIGNMENT)
    print("ANSWER: {}".format(CDCL_CNF_ASSIGNMENT))
    print("CDCL OPERATIONS: {}".format(c.CDCL_OPS))
    print("ORIGINAL CNF: {}".format(PLS_CDCL_TWO))
    '''

    #TESTING STUFF FOR IMPLICATION GRAPH
    '''
    start_time = time.time()
    print("BEFORE: {}".format(PLS_CDCL_TWO))
    ps = []
    PLS_CDCL_TWO = c2.SAT_REDUCTION(PLS_CDCL_TWO,1,0,ps,True)
    print("AFTER: {}".format(PLS_CDCL_TWO))
    print("TEMPORARY ASSIGNMENT: {}".format(ps) )
    c2.BOOLEAN_CONSTRAINT_PROPAGATION(PLS_CDCL_TWO, True)
    c2.construct_implication_graph(ps)
    c2.display_implication_graph()
    print("ELAPSED TIME: {}".format(time.time() - start_time) )
    '''

    start_time = time.time()
    ans = c.DPLL_Heuristic(PLS_CDCL_TWO,[])
    ans2 = c2.CDCL_Heuristic(PLS_CDCL)
    print("ANS: {}".format(ans))
    print("ELAPSED TIME: {}".format(time.time()-start_time))

TEST_SET_NONINCREMENTAL()



'''
v = VideoCapture(0)
while(True):
    ret, frame = v.read()
'''
