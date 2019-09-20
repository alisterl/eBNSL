#!/usr/bin/env python3


__docformat__ = 'restructuredtext'

from math import lgamma, log

import numpy as np

from itertools import combinations

import sys

from numba import jit, njit

import argparse

from scipy.special import digamma

import pandas as pd
        


# START functions for upper bounds

@jit(nopython=True)
def lg(n,x):
    return lgamma(n+x) - lgamma(x)

def h(counts):
    '''
    log P(counts|theta) where theta are MLE estimates computed from counts
    ''' 
    tot = sum(counts)
    #print(tot)
    res = 0.0
    for n in counts:
        if n > 0:
            res += n*log(n/tot)
    return res

def chisq(counts):
    tot = sum(counts)
    t = len(counts)
    mean = tot/t #Python 3 - this creates a float
    chisq = 0.0
    for n in counts:
        chisq += (n - mean)**2
    return chisq/mean

def hsum(distss):
    res = 0.0
    for dists in distss:
        res += sum([h(d) for d in dists])
    return res

@jit(nopython=True)
def fa(dist,sumdist,alpha,r):
    #res = -lg(sumdist,alpha)
    res = lgamma(alpha) - lgamma(sumdist+alpha)
    alphar = alpha/r
    k = 0
    for n in dist:
        if n > 0:
            #res += lg(n,alphar)
            #res += (lgamma(n+alphar) - lgamma(alphar))
            res += lgamma(n+alphar)
            k += 1
    return res - k*lgamma(alphar)

def diffa(dist,alpha,r):
    """Compute the derivative of local-local BDeu score

    numba does not support the scipy.special functions -- we are using the
    digamma function in defining the entropy of the Chisq
    distribution. For this end I added a python script which contains code
    for a @vectorize-d digamma function. If we need to use anything from
    scipy.special we will have to write it up ourselves.

    """
    args = [n+alpha/r for n in dist] + [alpha,sum(dist)+alpha,alpha/r]
    z = digamma(args)
    return sum(z[:r+1]) - z[r+1] - r*z[r+2] 

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

@njit
def array_equal(a,b,pasize):
    for i in range(pasize):
        if a[i] != b[i]:
            return False
    return True

#@njit
#def get_elems(a,idxs):
#    return np.a

@njit
def upper_bound_james_fun(atoms_ints,atoms_floats,pasize,alpha,r,idxs):
    ub = 0.0
    local_ub = 0.0
    best_diff = 0.0
    pasize_r = pasize+r
    lr = -log(r)
    oldrow = atoms_ints[idxs[0]]
    for i in idxs:
        row = atoms_ints[i]
        
        if not array_equal(oldrow,row,pasize):
            ub += min(local_ub+best_diff,lr)
            best_diff = 0.0
            local_ub = 0.0
            oldrow = row

        mle = atoms_floats[i]
        local_ub += mle
        if row[-1]: # if chi-sq condition met
            diff = fa(row[pasize:pasize_r],row[pasize_r],alpha/2.0,r) - mle
            if diff < best_diff:
                best_diff = diff
    ub += min(local_ub+best_diff,lr)
    return ub


#@jit(nopython=True)
def ub(dists,alpha,r):
    '''
    Args:
        dists (iter): list of (child-value-counts,chisqtest,mles) lists, for some
          particular instantiation ('inst1') of the current parents.
          There is a child-value-counts list for each non-zero instantiation ('inst2') of the biggest
          possible superset of the current set of parents. inst1 and each inst2 have the same
          values for the current parents.
        alpha (float): ESS/(num of current parent insts)
        r (int): arity of the child

    Returns:
       float: An upper bound on the BDeu local score for (implicit) child
        over all possible proper supersets of the current parents
    '''
    naives = 0.0
    best_diff = 0.0
    for (dist,sumdist,ok_first,naive_ub) in dists:
        naives += naive_ub
        #iffirst_ub = naive_ub
        #iffirst_ub = min(len(dist)*-log(r),naive_ub)
        if ok_first:
            diff = fa(dist,sumdist,alpha/2.0,r) - naive_ub
            #iffirst_ub = min(iffirst_ub,fa(dist,sumdist,alpha/2.0,r))
            #diff = iffirst_ub - naive_ub
            if diff < best_diff:
                best_diff = diff
    return best_diff + naives

# START functions for upper bounds


@jit(nopython=True)
def _process_tree_contab_inner(arities, contab_tree, tree_index, att_index, variables, alpha_div_arities):
    """
        Recursively processes a contingency table tree that is stored in an array.
        
        Parameters
        ----------
    
        
        - `arities`: contains the arities of every variable
        - `contab_tree`: the tree contingency table
        - `tree_index`: The current array index that is being examined
        - `att_index`: The index of the current attribute in variables
        - `variables`: The variables for which the contingency table was constructed
        - `alpha_div_arities`: Alpha divided by the product of the arities of the variables
    """
    score = 0.0
    non_zero_count = 0
    for val in range(arities[variables[att_index]]):
        next_tree_index = tree_index + val

        if contab_tree[next_tree_index] != 0:
            if att_index == variables.size - 1:
                # count value
                score -= lgamma(alpha_div_arities+contab_tree[next_tree_index])
                non_zero_count += 1
            else:
                add_to_score, add_to_non_zero_count = _process_tree_contab_inner(arities, contab_tree, contab_tree[next_tree_index], att_index+1, variables, alpha_div_arities)
                score += add_to_score
                non_zero_count += add_to_non_zero_count
    return score, non_zero_count
    
@jit(nopython=True)
def _process_tree_contab(arities, contab, variables, alpha_div_arities):
    """
        Recursively processes a contingency table tree that is stored in an array.
        Calculates the BDeu "score component" for that table.
        
        Parameters:
        
        - `arities`: contains the arities of every variable
        - `contab`: the tree contingency table
        - `variables`: The variables for which the contingency table was constructed
        - `alpha_div_arities`: Alpha divided by the product of the arities of the variables
    """
    score, non_zero_count = _process_tree_contab_inner(arities,contab, 0, 0, variables, alpha_div_arities)
    score += non_zero_count*lgamma(alpha_div_arities)  
    return score, non_zero_count

# This must have the same signature as _process_tree_contab
@jit(nopython=True)    
def _process_array_contab(arities, contab, variables, alpha_div_arities):
    non_zero_count = 0
    score = 0
    for count in contab:
        if count != 0:
            non_zero_count+=1
            score -= lgamma(alpha_div_arities+count) 
    
    score += non_zero_count*lgamma(alpha_div_arities)  
    return score, non_zero_count

# This was moved out of the class as it can easily be accelerated with Numba
@jit(nopython=True)
def _bdeu_score_component(contab_generator, arities, variables, alpha, process_contab_function):
    """
        Parameters:

        - `contab_generator`: Function used to generate the contingency tables for `variables`
        - `arities`: A NumPy array containing the arity of every variable
        - `variables`: A Numpy array of ints which are column indices for the variables for which the score component is being constructed
        - `alpha`: Effective sample size
        - `process_contab_function`: Function called to do main computation

        - Returns (score, non_zero_count) where score is the BDeu component score and non_zero_count
          is the number of cells with positive count in the contingency table for `variables`
    """
    alpha_div_arities = alpha / arities[variables].prod()
    
    contab = contab_generator.make_contab(variables)
    score, non_zero_count = process_contab_function(arities, contab, variables, alpha_div_arities)             
    
    return score, non_zero_count

#@jit(nopython=True)
def get_atoms(data,arities):
    '''
    Args: 
        data(np.array): Discrete data as a 2d array of ints

    Returns:
        list: a list `atoms`, where `atoms[i]` is a dictionary mapping instantations
         of variables other than i to a tuple with 3 elements:

          1. the child counts (n1, n2, .., nr) for that instantations
          2. n * the entropy for the empirical distribution (n1/n, n2/n, .., nr/n)
          3. Whether the chi-squared statistic for (n1, n2, .., nr) exceeds r-1 
        
        Chi-squared test comes from "Compound Multinomial Likelihood Functions are Unimodal: 
        Proof of a Conjecture of I. J. Good". Author(s): Bruce Levin and James Reeds
        Source: The Annals of Statistics, Vol. 5, No. 1 (Jan., 1977), pp. 79-87
    
        At least two of the ni must be positive for the inst-tuple pair to be included
        in the dictionary since only in that case is the inst-tuple useful for getting
        a good upper bound.
       
    '''
    fullinsts = []
    for i in range(data.shape[1]):
        fullinsts.append({})
    for row in data:
        row = tuple(row)
        for i, val in enumerate(row):
            fullinsts[i].setdefault(row[:i]+(0,)+row[i+1:],  # add dummy value for ith val
                                    [0]*arities[i])[val] += 1

    # now, for each child i, delete full insts which are 'deterministic'
    # i.e. where only one child value is non-zero
    newdkts_ints = []
    newdkts_floats = []
    for i, dkt in enumerate(fullinsts):
        #print('old',len(dkt))
        #newdkt = {}
        newdkt_ints = []
        newdkt_floats = []
        th = arities[i] - 1
        for inst, childcounts in dkt.items():
            if len([v for v in childcounts if v > 0]) > 1:
                newdkt_ints.append(inst+tuple(childcounts)+(sum(childcounts),chisq(childcounts) > th))
                #newdkt[inst] = (
                #    np.array(childcounts,np.uint64),sum(childcounts),
                #    chisq(childcounts) > th,
                #    h(childcounts))
                newdkt_floats.append(h(childcounts))
        newdkts_ints.append(np.array(newdkt_ints,np.uint64))
        newdkts_floats.append(np.array(newdkt_floats,np.float64))
        #print('new',len(newdkt))
    #fullinsts = newdkts

    #print(newdkts_ints[0])
    #sys.exit()

    #for x in newdkts_ints:
    #    print(x.shape)
    
    return newdkts_ints, newdkts_floats
    

def save_local_scores(local_scores, filename):
    variables = local_scores.keys()
    with open(filename, "w") as scores_file:
        scores_file.write(str(len(variables)))
        for child, dkt in local_scores.items():
            scores_file.write("\n" + child + " " + str(len(dkt.keys())))
            for parents, score in dkt.items():
                #highest_sup = None
                scores_file.write("\n" + str(score) + " " + str(len(parents)) +" "+ " ".join(parents))

def prune_local_score(this_score, parent_set, child_dkt):
    for other_parent_set, other_parent_set_score in child_dkt.items():
        if other_parent_set_score >= this_score and other_parent_set < parent_set:
            return True
    return False

def fromdataframe(df):
    cols = []
    arities = []
    varnames = []
    for varname, vals in df.items():
        varnames.append(varname)
        cols.append(vals.cat.codes)
        arities.append(len(vals.cat.categories))
    return np.transpose(np.array(cols,dtype=np.uint32)), arities, varnames

class BDeuScoresGenerator:
    """
    Used to calcuate BDeuScores for a dataset of complete discrete data
    """
    
    def __init__(self, data_source, variables = None, arities = None,
                 use_adtree = False, rmin=32, max_palim_size=None):
        '''Initialises a `BDeuScoresGenerator` object.

        If  `data_source` is a filename then it is assumed that:

            #. All values are separated by whitespace
            #. Comment lines start with a '#'
            #. The first line is a header line stating the names of the 
               variables
            #. The second line states the arities of the variables
            #. All other lines contain the actual data

        Args:
          data_source (str/array_like/Pandas.DataFrame) : 
            Either a filename containing the data or an array_like object or
            Pandas data frame containing it.

          variables (iter) : 
           Variable names corresponding to columns in the data.
           Ignored if `data_source` is a filename or Pandas DataFrame (since they 
           will supply the variable names). Otherwise if not supplied (`=None`)
           then variables names will be: X1, X2, ...

          arities (iter) : 
           Arities for the variables corresponding to columns in the data.
           Ignored if `data_source` is a filename or Pandas DataFrame (since they 
           will supply the arities). Otherwise if not supplied (`=None`)
           the arity for each variable will be set to the number of distinct values
           observed for that variable in the data.

          use_adtree (bool) : 
           If True an ADTree will be used to compute the necessary
           counts. Otherwise a simpler approach is taken (it is
           recommended that you only use an ADTree if you have a large
           number of records in your dataset).
            
          rmin (int) :
           The rmin value for the ADTree.
           (Only used when using an ADTree.)
        
          max_palim_size (int/None) : 
           The maximum size of any parent set. 
           Setting this can lead to faster computation.
           (Only used when using an ADTree.)
        '''

        if type(data_source) == str:
            with open(data_source, "r") as file:
                line = file.readline().rstrip()
                while line[0] == '#':
                    line = file.readline().rstrip()
                variables = line.split()
                line = file.readline().rstrip()
                while line[0] == '#':
                    line = file.readline().rstrip()
                arities = np.array([int(x) for x in line.split()],dtype=np.uint32)

                for arity in arities:
                    if arity < 2:
                        raise ValueError("This line: '{0}' is interpreted as giving variable arities but the value {1} is less than 2.".format(line,arity))

                # class whose instances are callable functions 'with memory'
                class Convert:
                    def __init__(self):
                        self._last = 0 
                        self._dkt = {}

                    def __call__(self,s):
                        try:
                            return self._dkt[s]
                        except KeyError:
                            self._dkt[s] = self._last
                            self._last += 1
                            return self._dkt[s]


                converter_dkt = {}
                for i in range(len(variables)):
                    # trick to create a function 'with memory'
                    converter_dkt[i] = Convert()
                data = np.loadtxt(file,
                                  dtype=np.uint32,
                                  converters=converter_dkt,
                                  comments='#')

        elif type(data_source) == pd.DataFrame:
            data, arities, variables = fromdataframe(data_source)
        else:
            data = np.array(data_source,dtype=np.uint32)
        self._data = data
        self._data_length = data.shape[0]
        if arities is None:
            self._arities = np.array([x+1 for x in data.max(axis=1)],dtype=np.uint32)
        else:
            self._arities = np.array(arities,dtype=np.uint32)
        if variables is None:
            self._variables = ['X{0}'.format(i) for i in range(1,len(self._arities)+1)]
        else:
            # order of variables determined by header line in file, if file used
            self._variables = variables
        self._varidx = {}
        for i, v in enumerate(self._variables):
            self._varidx[v] = i

        # When Numba compiled classes are imported
        # they are compiled regardless of whether or not they are used.
        # Therefore we only want to import the class that we will actually be
        # using. (Note: it appears that importing a class multiple times does not 
        # lead to it being compiled multiple times.)
        if use_adtree:
            from ADTree import ADTree
            self.process_contab_function = _process_array_contab
            if max_palim_size == None:
                max_palim_size = self._arities.size - 1
            self._contab_generator = ADTree(data, self._arities, rmin, max_palim_size+1)
        else:
            from contab_simple_tree import ContabGenerator
            self.process_contab_function = _process_tree_contab
            self._contab_generator = ContabGenerator(data, self._arities)

        self._atoms = get_atoms(data,self._arities)
        #print(self._atoms)
        #self._bdeu_global_ubs = self.get_bdeu_global_ubs()
        #self._almost_atoms = self.get_almost_atoms()
        #self._bdeu_almost_global_ubs = self.get_bdeu_almost_global_ubs()
        #print(self._bdeu_global_ubs)
        #print(self._bdeu_almost_global_ubs)
        #for v in self._almost_atoms:
        #    for w in v:
        #        print(w)
        #    print()

    def upper_bound_james(self,child,parents,alpha=1.0):
        """
        Compute an upper bound on proper supersets of parents

        Args:
         child (str) : Child variable.
         parents (iter) : Parent variables
         alpha (float) : ESS value for BDeu score

        Returns:
         float : An upper bound on the local score for parent sets
         for `child` which are proper supersets of `parents`

        """
        child_idx = self._varidx[child]
        pa_idxs = sorted([self._varidx[v] for v in parents])
        for pa_idx in pa_idxs:
            alpha /= self._arities[pa_idx]
        r = self._arities[child_idx]

        # each element of atoms_ints is a tuple of ints:
        # (fullinst,childvalcounts,sum(childvalcounts),ok_first)
        # each element of atoms_floats is:
        # sum_n n*log(n/tot), where sum is over childvalcounts
        # and tot = sum(childvalcounts)
        atoms_ints, atoms_floats = self._atoms[0][child_idx], self._atoms[1][child_idx]

        if len(atoms_floats) == 0:
            return 0.0

        # remove cols corresponding to non-parents and order
        p = len(self._arities)
        end_idxs = list(range(p,p+r+2))
        #print(end_idxs)
        atoms_ints_redux = atoms_ints[:,pa_idxs+end_idxs]
        if len(pa_idxs) > 0:
            idxs = np.lexsort([atoms_ints_redux[:,col] for col in range(len(pa_idxs))])
        else:
            idxs = list(range(len(atoms_floats)))

        #print(atoms_ints)

        #pa_idxs.append(1) # add dummy element
        return upper_bound_james_fun(atoms_ints_redux,atoms_floats,len(pa_idxs),alpha,r,idxs)
        
        # this_ub = 0.0
        # for dists in self._atoms_for_parents(child,parents).values():
        #     this_ub += ub(dists,alpha,r)
        # return this_ub

        
    def _atoms_for_parents(self,child,parents):
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
                
        
        # now atoms for the same parent inst are grouped together
        
        # dkt = {}
        # for fullinst, childcounts in self._atoms[child_idx].items():
        #     inst = tuple([fullinst[i] for i in pa_idxs])
        #     try:
        #         dkt[inst].append(childcounts)
        #     except KeyError:
        #         dkt[inst] = [childcounts]
        # return dkt


    def get_almost_atoms(self):
        #alldkts[i][j] is the counts for child i when j is absent.
        p = len(self._atoms)
        alldkts = []
        for i, dkt in enumerate(self._atoms):
            arity = self._arities[i]
            newdkts = []
            for j in range(p):
                newdkt = {}
                # if jth not in parent set ...
                if j != i:
                    for inst, counts in dkt.items():
                        lessinst = inst[:j]+(0,)+inst[j+1:] # add dummy value again
                        if lessinst in newdkt:
                            # add counts
                            for k in range(arity):
                                newdkt[lessinst][k] += counts[k]
                        else:
                            newdkt[lessinst] = counts[:]
                newdkts.append(newdkt)
            alldkts.append(newdkts)
        return alldkts

    def get_bdeu_almost_global_ubs(self):
        allubs = []
        p = len(self._almost_atoms)
        for i, dkts in enumerate(self._almost_atoms):
            ubs = []
            for j in range(p):
                ub = 0.0
                for counts in dkts[j].values():
                    tot = sum(counts)
                    mle = 0.0
                    for c in counts:
                        if c > 0:
                            mle = c * log(c/tot)
                    ub += mle
                ubs.append(ub)
            allubs.append(ubs)
        return allubs
    
    def get_bdeu_global_ubs(self):
        ubs = []
        for dkt in self._atoms:
            ub = 0.0
            for counts in dkt.values():
                tot = sum(counts)
                mle = 0.0
                for c in counts:
                    if c > 0:
                        mle = c * log(c/tot)
                ub += mle
            ubs.append(ub)
        return ubs
                        
    def data(self):
        '''
        The data with all values converted to unsigned integers.

        Returns:
         pandas.DataFrame: The data
        '''

        df = pd.DataFrame(self._data,columns=self._variables)
        arities = self._arities
        for i, (name, data) in enumerate(df.items()):
            # ensure correct categories are recorded even if not
            # all observed in data
            df[name] = pd.Categorical(data,categories=range(arities[i]))
        return df

    def rawdata(self):
        '''
        The data with all values converted to unsigned integers.
        And without any information about variable names or aritiies.

        Returns:
         numpy.ndarray: The data
        '''
        return self._data
    
    def data_length(self):
        '''
        Returns:
         int: The number of datapoints in the data
        '''
        return self._data_length

    def arities(self):
        '''
        Returns:
         numpy.ndarray: The arities of the variables.
        '''
        return self._arities

    def arity(self,v):
        '''
        Args:
         v (str) : A variable name
        Returns:
         int : The arity of `v`
        '''

        return self._arities[self._varidx[v]]

    def variables(self):
        '''
        Returns:
         list : The variable names
        '''
        return self._variables

    def varidx(self):
        '''
        Returns:
         dict : Maps a variable name to its position in the list of variable names.
        '''
        return self._varidx

#    def ll_score(self,child,parents):
#        if len(parents) == 0:
#            return 0.0
#        variables = np.array([self._varidx[x] for x in list(parents)+[child]], dtype=np.uint32)
#        contab = self._contab_generator.make_contab(
    
    def bdeu_score_component(self,alpha,variables):
        '''Compute the BDeu score component for a set of variables
        (from the current dataset).

        The BDeu score for a child v having parents Pa is 
        the BDeu score component for Pa subtracted from that for v+Pa
        
        Args:
         alpha (float) : The effective sample size parameter for the BDeu score
         variables (iter) : The names of the variables

        Returns:
         float : The BDeu score component.
        '''
        
        if len(variables) == 0:
            return lgamma(alpha) - lgamma(alpha + self._data_length), 1
        else:
            variables = np.array([self._varidx[x] for x in list(variables)], dtype=np.uint32)
            variables.sort() #AD tree requires variables to be in order
            return _bdeu_score_component(self._contab_generator, self._arities, variables, alpha, self.process_contab_function)
    
    def bdeu_scores(self,alpha=1.0,palim=None,pruning=True):
        """
        Exhaustively compute all BDeu scores for all child variables and all parent sets up to size `palim`.
        If `pruning` delete those parent sets which have a subset with a better score.
        Return a dictionary dkt where dkt[child][parents] = bdeu_score
        
        Args:
         alpha (float) : ESS for BDeu score
         palim (int/None) : Limit on parent set size
         pruning (bool) : Whether to prune

        Returns:
         dict : dkt where dkt[child][parents] = bdeu_score
        """
        
        if palim == None:
            palim = self._arities.size - 1

        score_dict = {}
        # Initialisation
        # Need to create dict for every child
        # also its better to do the zero size parent set calc here
        # so that we don't have to do a check for every parent set
        # to make sure it is not of size 0 when calculating score component size
        no_parents_score_component = lgamma(alpha) - lgamma(alpha + self._data_length)
        for c, child in enumerate(self._variables):
            score_dict[child] = {
                frozenset([]):
                no_parents_score_component
            }
        
        for pasize in range(1,palim+1):
            for family in combinations(self._variables,pasize): 
                
                variables = np.array([self._varidx[x] for x in list(family)], dtype=np.uint32)
                score_component = _bdeu_score_component(self._contab_generator, self._arities, variables, alpha, self.process_contab_function)
                family_set = frozenset(family)
                for child in self._variables:
                    if child in family_set:
                        parent_set = family_set.difference([child])
                        score_dict[child][parent_set] -= score_component
                        if pruning and prune_local_score(score_dict[child][parent_set],parent_set,score_dict[child]):
                            del score_dict[child][parent_set]
                    else:
                        score_dict[child][family_set] = score_component 
                
                    
        for vars in combinations(self._variables,palim+1):
            variables = np.array([self._varidx[x] for x in list(vars)], dtype=np.uint32)
            vars_set = frozenset(vars)
            score_component = _bdeu_score_component(self._contab_generator, self._arities, variables, alpha, self.process_contab_function)
            for child in vars:
                parent_set = vars_set.difference([child])
                score_dict[child][parent_set] -= score_component
                if pruning and prune_local_score(score_dict[child][parent_set],parent_set,score_dict[child]):
                    del score_dict[child][parent_set]

            
        return score_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BDeu local score generation (for Bayesian network learning)')
    parser.add_argument('data_file', help="Data file in GOBNILP format")
    parser.add_argument('scores_file', help="Where to store the generated scores")
    parser.add_argument('-p','--palim', type=int, help="Parent set size limit")
    parser.add_argument('-a', '--alpha', type=float, default=1.0, help="The effective sample size")
    parser.add_argument('-t', '--adtree', action="store_true", help="Use and ADTree")
    parser.add_argument('-r', '--rmin', type=int, default=32, help = "Rmin value for the ADTree")

    args = parser.parse_args()
    
    local_scores_generator = BDeuScoresGenerator(args.data_file, use_adtree = args.adtree, rmin=args.rmin, max_palim_size=args.palim)
    local_scores = local_scores_generator.bdeu_scores(args.alpha,args.palim)
    save_local_scores(local_scores, args.scores_file)
    
