from matplotlib import pyplot as plt
import math
import numpy

class HMM(object):
    ''' Simple Hidden Markov Model implementation.  User provides
        transition, emission and initial probabilities in dictionaries
        mapping 2-character codes onto floating-point probabilities
        for those table entries.  States and emissions are represented
        with single characters.  Emission symbols comes from a finite.  '''
    
    def __init__(self, A, E, I):
        ''' Initialize the HMM given transition, emission and initial
            probability tables. '''
        
        # put state labels to the set self.Q
        self.Q, self.S = set(), set() # states and symbols
        for a, prob in A.items():
            asrc, adst = a[0], a[1]
            self.Q.add(asrc)
            self.Q.add(adst)
        # add all the symbols to the set self.S
        for e, prob in E.items():
            eq, es = e[0], e[1]
            self.Q.add(eq)
            self.S.add(es)
        
        self.Q = sorted(list(self.Q))
        self.S = sorted(list(self.S))
        
        # create maps from state labels / emission symbols to integers
        # that function as unique IDs
        qmap, smap = {}, {}
        for i in range(len(self.Q)): qmap[self.Q[i]] = i
        for i in range(len(self.S)): smap[self.S[i]] = i
        lenq = len(self.Q)
        
        # create and populate transition probability matrix
        self.A = numpy.zeros(shape=(lenq, lenq), dtype=float)
        for a, prob in A.items():
            asrc, adst = a[0], a[1]
            self.A[qmap[asrc], qmap[adst]] = prob
        # make A stochastic (i.e. make rows add to 1)
        self.A /= self.A.sum(axis=1)[:, numpy.newaxis]
        
        # create and populate emission probability matrix
        self.E = numpy.zeros(shape=(lenq, len(self.S)), dtype=float)
        for e, prob in E.items():
            eq, es = e[0], e[1]
            self.E[qmap[eq], smap[es]] = prob
        # make E stochastic (i.e. make rows add to 1)
        self.E /= self.E.sum(axis=1)[:, numpy.newaxis]
        
        # initial probabilities
        self.I = [ 0.0 ] * len(self.Q)
        for a, prob in I.items():
            self.I[qmap[a]] = prob
        # make I stochastic (i.e. adds to 1)
        self.I = numpy.divide(self.I, sum(self.I))
        
        self.qmap, self.smap = qmap, smap
        
        # Make log-base-2 versions for log-space functions
        self.Alog = numpy.log2(self.A)
        self.Elog = numpy.log2(self.E)
        self.Ilog = numpy.log2(self.I)
    
    def jointProb(self, p, x):
        ''' Return joint probability of path p and emission string x '''
        p = list(map(self.qmap.get, p)) # turn state characters into ids
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        tot = self.I[p[0]] # start with initial probability
        for i in range(1, len(p)):
            tot *= self.A[p[i-1], p[i]] # transition probability
        for i in range(0, len(p)):
            tot *= self.E[p[i], x[i]] # emission probability
        return tot
    
    def jointProbL(self, p, x):
        ''' Return log2 of joint probability of path p and emission
            string x.  Just like self.jointProb(...) but log2 domain. '''
        p = list(map(self.qmap.get, p)) # turn state characters into ids
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        tot = self.Ilog[p[0]] # start with initial probability
        for i in range(1, len(p)):
            tot += self.Alog[p[i-1], p[i]] # transition probability
        for i in range(0, len(p)):
            tot += self.Elog[p[i], x[i]] # emission probability
        return tot
    
    def viterbi(self, x):
        ''' Given sequence of emissions, return the most probable path
            along with its probability. '''
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        nrow, ncol = len(self.Q), len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # backtrace
        # Fill in first column
        for i in range(0, nrow):
            mat[i, 0] = self.E[i, x[0]] * self.I[i]
        # Fill in rest of prob and Tb tables
        for j in range(1, ncol):
            for i in range(0, nrow):
                ep = self.E[i, x[j]]
                mx, mxi = mat[0, j-1] * self.A[0, i] * ep, 0
                for i2 in range(1, nrow):
                    pr = mat[i2, j-1] * self.A[i2, i] * ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi
        # Find final state with maximal probability
        omx, omxi = mat[0, ncol-1], 0
        for i in range(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Backtrace
        i, p = omxi, [omxi]
        for j in range(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = ''.join(map(lambda x: self.Q[x], p[::-1]))
        return omx, p # Return probability and path
    
    def viterbiL(self, x):
        ''' Given sequence of emissions, return the most probable path
            along with log2 of its probability.  Just like viterbi(...)
            but in log2 domain. '''
        x = list(map(self.smap.get, x)) # turn emission characters into ids
        nrow, ncol = len(self.Q), len(x)
        mat   = numpy.zeros(shape=(nrow, ncol), dtype=float) # prob
        matTb = numpy.zeros(shape=(nrow, ncol), dtype=int)   # backtrace
        # Fill in first column
        for i in range(0, nrow):
            mat[i, 0] = self.Elog[i, x[0]] + self.Ilog[i]
        # Fill in rest of log prob and Tb tables
        for j in range(1, ncol):
            for i in range(0, nrow):
                ep = self.Elog[i, x[j]]
                mx, mxi = mat[0, j-1] + self.Alog[0, i] + ep, 0
                for i2 in range(1, nrow):
                    pr = mat[i2, j-1] + self.Alog[i2, i] + ep
                    if pr > mx:
                        mx, mxi = pr, i2
                mat[i, j], matTb[i, j] = mx, mxi
        # Find final state with maximal log probability
        omx, omxi = mat[0, ncol-1], 0
        for i in range(1, nrow):
            if mat[i, ncol-1] > omx:
                omx, omxi = mat[i, ncol-1], i
        # Backtrace
        i, p = omxi, [omxi]
        for j in range(ncol-1, 0, -1):
            i = matTb[i, j]
            p.append(i)
        p = ''.join(map(lambda x: self.Q[x], p[::-1]))
        return omx, p # Return log probability and path