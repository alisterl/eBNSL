from numba import jitclass, jit, uint32
import numpy as np

contab_generator_spec = [
    ('_data',uint32[:,:]),
    ('_arities',uint32[:]),
]

@jitclass(contab_generator_spec)
class ContabGenerator(object):
    def __init__(self, data, arities):
        
        self._data = data
        self._arities = arities
    
    def make_contab(self, attributes):
        """
            Returns a contingency table for the attributes
            that only contains the counts (it is just a 1D array
            of numbers)
            Note: attributes must be sorted in increasing order!
        """
        
        return _make_contab(self, attributes)
        
    def make_contab_full(self, attributes):
        """
            Returns a contingency table for the given attribute
            There is one column for each attribute that denotes
            the value that that attribute has taken. The last column
            denotes the number of records that there were with that
            combination of attribute values
            Note: attributes must be sorted in increasing order!
        """
        # Could maybe take this method outside the class so that the
        # parallel compilation option could be used but as its
        # not actually used when computing BDeu scores I didn't
        
        skinny_contab = self.make_contab(attributes)
        sub_contab_size = np.empty(attributes.size, np.uint32)
        sub_contab_size[-1] = 1

        for i in range(attributes.size -2, -1, -1):
            sub_contab_size[i] = self._arities[attributes[i+1]] * sub_contab_size[i+1]

        contab_height = self._arities[attributes[0]] * sub_contab_size[0]
            
        # 1 column per variable + 1 column for the count
        contab_width = attributes.size + 1
        
        contab = np.empty((contab_height,contab_width), dtype=np.uint32)
        repeat_each = contab_height
        contab[0:repeat_each,-1] = 0
        tile = 1
        for attribute_index in np.arange(attributes.size):
            arity = self._arities[attributes[attribute_index]]
            repeat_each //= arity
            for t in np.arange(0,contab_height,contab_height//tile,np.uint32):
                offset = 0
                offset_upper = 0
                for val in np.arange(0,arity,1,np.uint32):
                    offset_upper += repeat_each
                    contab[t+offset:t+offset_upper, attribute_index] = val
                    offset = offset_upper
            # Would do this but its not supported in numba...
            #contab[0:contab_height,attribute_index] = np.tile(np.repeat(np.arange(0, arity, 1, np.uint32), repeat_each), tile)
            
            tile *= arity
        
        fill_contab(self,contab,skinny_contab, 0, 0, 0, attributes, sub_contab_size)
            
        return contab   

@jit(nopython=True)
def fill_contab(contab_generator, contab, tree, tree_index, att_index, contab_row, attributes, sub_contab_size):
    arity = contab_generator._arities[attributes[att_index]]
    for val in range(arity):
        next_tree_index = tree_index + val

        if tree[next_tree_index] != 0:
            if att_index == attributes.size - 1:
                # count value
                contab[contab_row+val,-1] = tree[next_tree_index]
            else:
                fill_contab(contab_generator, contab, tree, tree[next_tree_index], att_index+1, contab_row+(sub_contab_size[att_index]*val), attributes, sub_contab_size)
            
        
@jit(nopython=True)
def _double_array_size(array):
    """
        We could just use array.resize, but that isn't supported by Numba...
        Returns an array that contains the same data as the parameter array
        (these elements are placed at the begining of the new array)
        but is twice as large
    """
    #Array.resize not available with Numba :(
    size = array.size
    new = np.empty(size * 2, dtype=array.dtype)
    new[:size] = array
    return new         

@jit(nopython=True)
def _make_contab(contab_generator, attributes):
    """
        Had to be separate from the class so that the parallel compilation option could be used.
        However it won't compile proprely when the parallel option is used...
    """
    
    real_size = contab_generator._arities[attributes[0]]
    contab_for_counts = np.empty(contab_generator._arities[attributes].sum(), dtype=np.uint32)
    contab_for_counts[0:real_size] = 0
    
    
    last_att = attributes[-1]
    for row in range(contab_generator._data.shape[0]):
        index = 0
        for att in range(attributes.size - 1):
            index += contab_generator._data[row,attributes[att]]
            if contab_for_counts[index] == 0:
               
                contab_for_counts[index] = real_size
                new_real_size = real_size + contab_generator._arities[attributes[att+1]]
                if new_real_size > contab_for_counts.size:
                    contab_for_counts = _double_array_size(contab_for_counts)
                contab_for_counts[real_size:new_real_size] = 0
                real_size = new_real_size
            index = contab_for_counts[index]
        contab_for_counts[index + contab_generator._data[row,last_att]]+=1     
        
    return contab_for_counts