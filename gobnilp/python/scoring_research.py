#!/usr/bin/env python3


__docformat__ = 'restructuredtext'

import pandas as pd

from math import lgamma, log
from collections import Counter

from itertools import combinations
import operator

import sys

epsilon = 0.01

def lg(n,x):
    return lgamma(n+x) - lgamma(x)

def h(dist):
    tot = sum(dist)
    res = 0.0
    for n in dist:
        if n > 0:
            res += n*log(n/tot)
    return res

def hsum(distss):
    res = 0.0
    for dists in distss:
        res += sum([h(d) for d in dists])
    return res

def fa(dist,alpha,r):
    res = -lg(sum(dist),alpha)
    alphar = alpha/r
    for n in dist:
        if n > 0:
            res += lg(n,alphar)
    return res

    
def diffa(dist,alpha,r):
    """
    Compute the derivative of local-local BDeu score
    """
    res = 0.0
    for n in dist:
        for i in range(n):
            res += 1.0/(i*r+alpha)
    for i in range(sum(dist)):
        res -= 1.0/(i+alpha)
    return res

def onepositive(dist):
    """
    Is there only one positive count in `dist`?
    """
    npos = 0
    for n in dist:
        if n > 0:
            npos += 1
            if npos > 1:
                return False
    return True

def ub(dists,alpha,r):
    naives = 0.0
    best_diff = 0.0
    for dist in dists:
        naive_ub = h(dist)
        naives += naive_ub
        iffirst_ub = min(len(dist)*-log(r),naive_ub)
        if not onepositive(dist) and diffa(dist,alpha,r) >= 0:
            iffirst_ub = min(iffirst_ub,fa(dist,alpha/2.0,r))
        diff = iffirst_ub - naive_ub
        if diff < best_diff:
            best_diff = diff
    return best_diff + naives
    
#def ub(dist,alpha,r):
#    # start with easy upper bound
#    ub = min(-log(r),h(dist))
#    if onepositive(dist) and diffa(dist,alpha,r) >= 0:
#        ub = min(ub,fa(dist,alpha/2.0,r))
#    return ub

class Data:
    """
    A dataset of complete discrete data
    """
    
    def __init__(self,fname):
        """"
        Parameters:

        - `fname`: Filename containing the data. It is assumed that:
        
        1. All values are separated by whitespace
        2. Comment lines start with a '#'
        3. The first line is a header line stating the names of the variables
        4. The second line states the arities of the variables
        5. All other lines contain the actual data
        """
        
        data = pd.read_table(fname,
                             delim_whitespace=True,
                             comment='#')
        arities = [int(x) for x in data.iloc[0]]
        self._data = data[1:]
        self._arities = dict(zip(list(self._data),arities))
        self._variables = list(self._arities.keys())
        self._varidx = {}
        for i, v in enumerate(self._variables):
            self._varidx[v] = i
        self.get_atoms()

    def upper_bound_james(self,child,parents,alpha=1.0):
        """
        Compute an upper bound on proper supersets of parents
        """
        for pa in parents:
            alpha /= self._arities[pa]
        r = self._arities[child]
        this_ub = 0.0
        for dists in self.atoms_for_parents(child,parents).values():
            this_ub += ub(dists,alpha,r)
        return this_ub

    def upper_bound_cassio(self,child,parents,aq,counts):
        """
        Cassio's upper bound
        """
        m = 0
        for child_counts in self.atoms_for_parents(child,parents).values():
            m += max([len(cc) for cc in child_counts])-1
        return -len(counts)*log(self._arities[child]) - log(2.0/aq + 1)*m
            

    
    def upper_bound_weak(self,child,posfamilyinsts):
        """
        Compute a weak upper bound on supersets of some parent set 
        for some family = child + parents, where
        - `child`: the child of the family
        - `posfamilyinsts`: the number of instantiations of the family 
        which occur at least once in the data
        """
        return -posfamilyinsts * log(self._arities[child])
    
    def bdeu_score(self,child,parents,alpha=1.0):
        """
        Returns a tuple (score, ubs) where
        - score is the BDeu score of a particular child and parent set
        - ubs is a dictionary of mapping the names of upper bounds to 
        upper bounds on the BDeu scores of supersets of the parent set.

        Parameters:

        - `child`: the child of the family
        - `parents`: the parents of the family (an iterable of variables)
        - `alpha`: the effective sample size
        """
        aq = alpha
        for parent in parents:
            aq /= self._arities[parent]
        aqr = aq / self._arities[child]
        counts = self._data.groupby(list(parents)+[child],sort=True).size()
        bdeu_score = 0.0
        if len(parents) == 0:
            nij = 0
            for nijk in counts:
                bdeu_score += lg(nijk,aqr)
                nij += nijk
            bdeu_score -= lg(nij,aq)
        else:
            cnt = Counter()
            for idx, nijk in counts.iteritems():
                cnt[idx[:-1]] += nijk
                bdeu_score += lg(nijk,aqr)
            for nij in cnt.values():
                bdeu_score -= lg(nij,aq)
        return bdeu_score, {'weak':self.upper_bound_weak(child,len(counts)), 'cassio':self.upper_bound_cassio(child,parents,aq,counts), 'james':self.upper_bound_james(child,parents,alpha)}

    def all_bdeu_scores(self,alpha=1.0,palim=3):
        """
        Exhaustively compute all BDeu scores and upper bounds for all families up to `palim`
        Return a dictionary dkt where dkt[child][parents] = bdeu_score
        """
        score_dict = {}
        for i, child in enumerate(self._variables):
            print(child)
            potential_parents = frozenset(self._variables[:i]+self._variables[i+1:])
            child_dkt = {}
            for pasize in range(palim+1):
                for parents in combinations(potential_parents,pasize):
                    child_dkt[frozenset(parents)] = self.bdeu_score(child,parents,alpha)
            score_dict[child] = child_dkt
        return score_dict

    def atoms_for_parents(self,child,parents):
        """
        Return a dictionary whose keys are instantiations of the parents
        with positive count in the data and whose values are lists of lists
        of lists of child counts. 

        If dkt is the returned dictionary, then e.g.
        dkt[0,1,0] = [ [1,2], [0,4] ]
        means that there are 3 parents and there are 2 full parent instantiations for 
        the instantiation (0,1,0) one with child counts [1,2] and one with child counts [0,4]
        
        An exception will be raised if the _atoms attribute has not yet been computed.
        """
        # get indices of parents in vector of all possible parents for child
        child_idx = self._varidx[child]
        parentset = frozenset(parents)
        pa_idxs = []
        for i, parent in enumerate(self._variables[:child_idx]+self._variables[child_idx+1:]):
            if parent in parentset:
                pa_idxs.append(i)
                
        dkt = {}
        for fullinst, childcounts in self._atoms[child].items():
            inst = tuple([fullinst[i] for i in pa_idxs])
            try:
                dkt[inst].append(childcounts)
            except KeyError:
                dkt[inst] = [childcounts]
        return dkt

                
    def get_atoms(self):
        """
        Compute a dictionary whose keys are child variables
        and whose values are dictionaries mapping instantiations
        of all the other parents to a list of counts for the child variable
        for that instantiation. Only parent set instantations with a positive
        count in the data are included.

        The dictionary is stored as the value of self._atoms
        """
        # assume that following is quicker than creating dictionary from scratch
        counts = self._data.groupby(self._variables).size()
        dktall = {}
        for child_idx, child in enumerate(self._variables):
            def keyfn(item): return item[0][:child_idx]+item[0][child_idx+1:]
            child_counts = sorted(counts.iteritems(),key=keyfn)
            last = None
            child_dist = []
            dkt = {}
            for inst, v in child_counts:
                inst = inst[:child_idx]+inst[child_idx+1:]
                if last is not None and inst != last:
                    dkt[last] = child_dist
                    child_dist = []
                last = inst
                child_dist.append(v)
            dkt[last] = child_dist
            dktall[child] = dkt
        self._atoms = dktall


    def pruned_local_scores(self,palim=3,verbose=False):
        """
        Return a dictionary for each child variable where the keys
        are the child variables and the values map parent sets to 
        local scores.

        Not all parent sets are included. Only those parent set of cardinality
        at most `palim` can be included. Also, if it can be established that
        a parent set can not be a parent set for the child in an optimal Bayesian 
        network, then it is not included.

        The data and the particular local score used are determined by `self`

        Args:
            palim(int): maximal size of a parent set
            verbose(bool): whether messages on progress should be printed
                           to standard output

        """
        score_dict = {}
        p = len(self._variables)
        palim = min(palim,p-1)
        for child_idx, child in enumerate(self._variables):
            score, ub = self.local_score(child,[])
            # 'True' means 'score' is the actual local score (of the empty parent set)
            # as opposed to a superior local score of some subset of the parent set
            previous_layer = {():(score,ub,True)}
            for pasize in range(palim):
                new_layer = {}
                for oldpa in previous_layer:
                    last_idx = -1 if oldpa == () else self._varidx[oldpa[-1]]
                    for newpa_idx in range(last_idx+1,p):
                        if newpa_idx == child_idx:
                            continue
                        parents = oldpa + (self._variables[newpa_idx],)

                        # check that all subsets exist in previous layer,
                        # and get best (i.e. highest) score and best (i.e. lowest) upper bound
                        bss = None
                        lub = None
                        for i in range(len(parents)):
                            try:
                                old_score, old_ub, _ = previous_layer[parents[:i]+parents[i+1:]]
                            except KeyError:
                                if verbose:
                                    print(parents[:i]+parents[i+1:],'is missing')
                                bss = None
                                break
                            bss = old_score if bss is None else max(bss,old_score)
                            lub = old_ub if lub is None else min(lub,old_ub)


                        if bss is None or bss > lub:
                            # some subset is exponentially pruned, so don't score
                            # or: best we can hope for for 'parents' (lub) is worse than some existing subset
                            # of parents (bss)
                            if verbose:
                                print('Prune!',child,parents,bss,lub)
                            continue

                        score, ub = self.local_score(child,parents)
                        
                        # shouldn't need the following line since ub should already be as low as possible
                        ub = min(ub,lub)

                        if ub <= bss:
                            # neither this parent set or any of its supersets worth keeping
                            if verbose:
                                print('Exponentially pruning:', child, parents, ub, bss)
                            continue
                      
                        new_layer[parents] = (max(score,bss),ub,(score>bss))

                # only store scores (val[0]) which are indeed the local scores of the relevant parents
                # val[2]==True if this is the case
                child_dkt.update({parents:val[0] for (parents,val) in new_layer.items() if val[2]})        
                previous_layer = new_layer
            score_dict[child] = child_dkt
        return score_dict

        
    def pruned_bdeu_scores(self,alpha=1.0,palim=3,verbose=False):
        """
        Return a dictionary for each child variable where the keys
        are the child variables and the values map parent sets to 
        BDeu scores.

        Not all parent sets are included. Only those parent set of cardinality
        at most `palim` can be included. Also, if it can be established that
        a parent set can not be a parent set for the child in an optimal Bayesian 
        network, then it is not included.

        Parameters:

        - `alpha`: the effective sample size (prior parameter)
        - `palim`: maximal size of a parent set
        - `verbose`: whether messages on progress should be printed
        to standard output
        """
        score_dict = {}
        p = len(self._variables)
        palim = min(palim,p-1)
        for child_idx, child in enumerate(self._variables):
            score, ubdkt = self.bdeu_score(child,[],alpha)
            ub = min(ubdkt.values())
            child_dkt = {():score}
            previous_layer = {():(score,ub,True)}
            for pasize in range(palim):
                new_layer = {}
                for oldpa in previous_layer:
                    last_idx = -1 if oldpa == () else self._varidx[oldpa[-1]]
                    for newpa_idx in range(last_idx+1,p):
                        if newpa_idx == child_idx:
                            continue
                        parents = oldpa + (self._variables[newpa_idx],)

                        # check that all subsets exist in previous layer,
                        # and get best score and upper bound
                        bss = None
                        lub = None
                        for i in range(len(parents)):
                            try:
                                old_score, old_ub, _ = previous_layer[parents[:i]+parents[i+1:]]
                            except KeyError:
                                if verbose:
                                    print(parents[:i]+parents[i+1:],'is missing')
                                bss = None
                                break
                            bss = old_score if bss is None else max(bss,old_score)
                            lub = old_ub if lub is None else min(lub,old_ub)


                        if bss is None or bss > lub:
                            # some subset is exponentially pruned, so don't score
                            # or: best we can hope for for 'parents' is worse than some existing subset
                            # of parents
                            if verbose:
                                print('Prune!',child,parents,bss,lub)
                            continue

                        score, ubdkt = self.bdeu_score(child,parents,alpha)
                        ub = min(ubdkt.values())
                        
                        # shouldn't need the following line since ub should already be as low as possible
                        ub = min(ub,lub)

                        if ub <= bss:
                            # neither this parent set or any of its supersets worth keeping
                            if verbose:
                                print('Exponentially pruning:', child, parents, ub, bss)
                            continue
                      
                        new_layer[parents] = (max(score,bss),ub,(score>bss))
                
                child_dkt.update({parents:val[0] for (parents,val) in new_layer.items() if val[2]})        
                previous_layer = new_layer
            score_dict[child] = child_dkt
        return score_dict

    
if __name__ == '__main__':
    fname = sys.argv[1]
    palim = int(sys.argv[2])
    d = Data(fname)
    all_scores = d.all_bdeu_scores(1.0,palim)
    datapoints = []
    if palim < len(all_scores)-1:
        get_actual_lub = False
    else:
        get_actual_lub = True
    for child, dkt in all_scores.items():
        print('\nChild=',child)
        for parents, (score,ubdkt) in dkt.items():
            highest_sup = None
            if get_actual_lub:
                for sup in dkt.keys():
                    if sup > parents:
                        if highest_sup is None or highest_sup < dkt[sup][0]:
                            highest_sup = dkt[sup][0]
            if highest_sup is not None and highest_sup > min(ubdkt.values()) + epsilon:
                print('ERROR!')
                sys.exit()
            ubdkt['true'] = highest_sup
            print(','.join(parents),'\nscore={0},ubs={1}\n'.format(score,ubdkt))
            datapoints.append(ubdkt)
    forplot = open(sys.argv[3],'w')
    if 'true' in ubdkt:
        def keyfn(x):
            tmp = x['true']
            return 0 if tmp is None else tmp
    else:
        def keyfn(x): return x['james']
    datapoints.sort(key=keyfn)
    ubnames = list(datapoints[0].keys())
    print(' '.join(ubnames),file=forplot)
    for dpt in datapoints:
        ubs = [dpt[name] for name in ubnames]
        print(' '.join([str(x) for x in ubs]),file=forplot)
    forplot.close()
    sys.exit()

    #print(d.bdeu_score('Cancer',['Smoking','Tuberculosis'],1.0))
    score_dict = d.pruned_bdeu_scores(verbose=False,palim=palim)
    print(len(score_dict))
    for child, dkt in score_dict.items():
        print(child,len(dkt))
        for parents, score in sorted(dkt.items(), key=operator.itemgetter(1), reverse=True):
            print(score, len(parents), ' '.join(parents))


