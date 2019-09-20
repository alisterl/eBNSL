'''

@author: Josh

Numba optimised ADtree implemented using dynamic arrays
'''

from numba import jitclass, jit, uint32, from_dtype, void, boolean
import numpy as np

ad_array_item_type = np.dtype([('count',np.uint32),('children',np.uint32)])
vary_array_item_type = np.dtype([('mcv',np.uint32),('children',np.uint32)])

ADTree_spec = [
    ('_ad_nodes',from_dtype(ad_array_item_type)[:]),
    ('_vary_nodes',from_dtype(vary_array_item_type)[:]),
    ('_leaf_pointers',uint32[:]),
    ('_max_attributes_in_contab', uint32),
    ('_data',uint32[:,:]),
    ('_arities',uint32[:]),
    ('_rmin',uint32),
    ('_num_ad_nodes',uint32),
    ('_num_vary_nodes',uint32),
    ('_num_leaf_pointers',uint32)
]

@jitclass(ADTree_spec)
class ADTree(object):
    """
        A Numba accellerated ADtree implementation that can be used to produce
        contingency tables (quickly)
        
        Note: we store the number of AD/vary nodes and leaf pointers
        in this class as they need to be passed by reference to
        every iteration of _build_vary_nodes
    """

    def __init__(self, data, arities, rmin, max_attributes_in_contab):
        """
            Parameters:
            - `data`: data set for which we are building this AD tree
            - `arities`: contains the arity for each variable in data
            - `rmin`: the minimum number of records that a node must refer to for it to be expanded (must be > 0)
            - `max_attributes_in_contab`: the maximum number of attributes that will appear in any contingency table
        """
        if rmin <= 0: # Messes up my code if its zero
            raise ValueError("rmin must be > 0")
      
        self._arities = arities
        self._max_attributes_in_contab = max_attributes_in_contab
        
        recordNums = np.arange(0, data.shape[0], 1, np.uint32)
        self._data = data
        
        # With this initial size we ensure that doubling the array will always mean that
        # there is enough space for a new set of nodes
        # have the root node, then will add at most arities.max()-1 at each step
        self._ad_nodes = np.empty(arities.max(),dtype=ad_array_item_type)
        self._ad_nodes[0]['count'] = recordNums.size
        self._ad_nodes[0]['children'] = 0
        
        # With this initial size we ensure that doubling the array will always mean that
        # there is enough space for a new set of nodes
        # will be wanting to add at most arities.size nodes at each step (in fact we
        # will certainly be adding less than this on every step but the first)
        self._vary_nodes = np.empty(arities.size,dtype=vary_array_item_type)
        
        self._num_vary_nodes = 0
        self._num_ad_nodes = 1
        
        self._rmin = rmin
        if recordNums.size < self._rmin:
            self._leaf_pointers = recordNums
            self._num_leaf_pointers = self._leaf_pointers.size
        else:
            self._leaf_pointers = np.empty(rmin -1, dtype=np.uint32)
            self._num_leaf_pointers = 0
            _build_vary_nodes(self, 0, recordNums, 0)
            self._leaf_pointers = self._leaf_pointers[:self._num_leaf_pointers]
        
        self._ad_nodes = self._ad_nodes[:self._num_ad_nodes]
        self._vary_nodes = self._vary_nodes[:self._num_vary_nodes]
       
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
        contab_height = skinny_contab.size
            
        # 1 column per variable + 1 column for the count
        contabWidth = attributes.size + 1
        
        contab = np.empty((contab_height,contabWidth), dtype=np.uint32)
        repeat_each = contab_height
        contab[0:repeat_each,-1] = 0
        tile = 1
        for attribute_index in np.arange(attributes.size):
            arity = self._arities[attributes[attribute_index]]
            repeat_each //= arity
            for t in np.arange(0,contab_height,contab_height//tile,np.uint32):
                offset = 0
                offsetUpper = 0
                for val in np.arange(0,arity,1,np.uint32):
                    offsetUpper += repeat_each
                    contab[t+offset:t+offsetUpper, attribute_index] = val
                    offset = offsetUpper
            # Would do this but its not supported in numba...
            #contab[0:contab_height,attribute_index] = np.tile(np.repeat(np.arange(0, arity, 1, np.uint32), repeat_each), tile)
            
            tile *= arity
                
              
        contab[:,-1] = skinny_contab
        
        
        return contab   

"""
Library functions - most had to be moved outside of the class due to Numba restrictions
For example if both node building functions were methods in the ADTree class it wouldn't compile
perhaps mutual recursion isn't supported with jitclass methods?
"""        

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
def _build_ad_nodes(ad, attribute_num, MCV, recordNums, num_attributes_in_branch):
    """
        Parameters:
        - `ad`: the ad tree that we are building
        - `attribute_num`: the attribute for which we are building the AD nodes
        - `MCV`: the most common value of this attribute
        - `recordNums`: an array containing "pointers" into the dataset for each record that matches the current "query" for the part of the tree that we are currently building
        - `num_attributes_in_branch`: the number of attributes that have been considered in the current branch from the root node (used to determine when we can stop building the tree)

    """
    # This is where the first new node will be placed
    insert_node_at = ad._num_ad_nodes
    
    ad._num_ad_nodes += ad._arities[attribute_num] - 1 # -1 as don't store MCV
    #self._tree_array.resize(self._tree_array.size + extendArrayBy)
    if ad._ad_nodes.size < ad._num_ad_nodes:
        ad._ad_nodes = _double_array_size(ad._ad_nodes)
    
    for each_attribute_value in np.arange(ad._arities[attribute_num]):
        # If not MCV and count is greater than 0
        #child_nums = [eachRecordNum for eachRecordNum in recordNums if self._data[eachRecordNum,attribute_num] == each_attribute_value]
        
        if each_attribute_value != MCV:
            child_nums = recordNums[ad._data[recordNums,attribute_num] == each_attribute_value]
            count = child_nums.size
            ad._ad_nodes[insert_node_at]['count']=count
                           
            # We may not need to bother building any more :)
            # as we shan't be generating contabs with more attributes than this
            # i.e. if max_attributes_in_contab is 4 and num_attributes_in_branch is 4
            # then this is as deep as we need to build because we won't go any
            # deeper than this when building contabs!
            if num_attributes_in_branch < ad._max_attributes_in_contab:
                #Note: count of 0 is a null node (no children)
                if count >= ad._rmin:
                    # Do not build leaf list
                    ad._ad_nodes[insert_node_at]['children'] = ad._num_vary_nodes
                    _build_vary_nodes(ad,attribute_num+1, child_nums, num_attributes_in_branch)
                elif count > 0:
                    # Build leaf list
                    if ad._leaf_pointers.size - ad._num_leaf_pointers < count:
                        ad._leaf_pointers = _double_array_size(ad._leaf_pointers)
                    
                    ad._leaf_pointers[ad._num_leaf_pointers:ad._num_leaf_pointers+count] = child_nums
                    ad._ad_nodes[insert_node_at]['children'] = ad._num_leaf_pointers
                    ad._num_leaf_pointers += count
            insert_node_at += 1
        
@jit(nopython=True)           
def _build_vary_nodes(ad, start_attribute_num, recordNums, num_attributes_in_branch):
    """
        Parameters:
        - `ad`: the ad tree that we are building
        - `start_attribute_num`: if we are building vary nodes it denotes the attribute from which to start constructing vary nodes (i.e. produce vary nodes for attribute 3 and up as the ones for attributes 0, 1 and 2 were built further up in the tree. 
        - `recordNums`: an array containing "pointers" into the dataset for each record that matches the current "query" for the part of the tree that we are currently building
        - `num_attributes_in_branch`: the number of attributes that have been considered in the current branch from the root node (used to determine when we can stop building the tree)
    """
    # Where the first new vary node will be placed
    insert_node_at = ad._num_vary_nodes
    ad._num_vary_nodes += ad._arities.size - start_attribute_num
    
    if ad._vary_nodes.size < ad._num_vary_nodes:
        ad._vary_nodes = _double_array_size(ad._vary_nodes)
                
    
    for attribute_num in np.arange(start_attribute_num, ad._arities.size,1,np.uint32):
        
        MCV = np.bincount(ad._data[recordNums, attribute_num]).argmax()
       
        ad._vary_nodes[insert_node_at]['mcv'] = MCV
        ad._vary_nodes[insert_node_at]['children'] = ad._num_ad_nodes
        
        _build_ad_nodes(ad,attribute_num, MCV, recordNums, num_attributes_in_branch+1)
        
        insert_node_at += 1    

    
@jit(void(ADTree.class_type.instance_type,uint32[:],uint32,uint32,uint32,uint32[:],uint32,uint32[:]),nopython=True, parallel=True)
def _make_contab_internal(ad, attributes, attribute_index, current_index, parent_vnode_number, contab, contab_row, sub_contab_size):
    """
        Had to be separate from the ADTree class to use the parallel compilation option
        also Numba wouldn't compile it without a signature
    """
    
    if attribute_index == attributes.size:
        contab[contab_row] = ad._ad_nodes[current_index]['count']
        return
    
    # AD node doesn't have children
    if ad._ad_nodes[current_index]['count'] == 0:
        # Nothing to do :)
        return    
    
    
    if ad._ad_nodes[current_index]['count'] < ad._rmin:
        # Leaf list
        count = ad._ad_nodes[current_index]['count']
        
        # Data rows that are pointed to
        # only the attributes that we are currently bothered bytearray
        # (which is those in attributes at index attribute_index or higher)
        useful_data = ad._data[
                ad._leaf_pointers[
                    ad._ad_nodes[current_index]['children']: 
                        ad._ad_nodes[current_index]['children'] + count
                ], :
            ][:,attributes[attribute_index:]]
        
        # For some reason with Numba taking a "normal" slice here breaks things
        # with the multiplication below
        # I think that it is because a normal slice [attribute_index:] creates
        # a view in to the original array and "fancy indexing" like below creates a copy
        # at some point in the future when Numba's been update try replacing subs with a normal
        # slice of sub_contab_size
        subs = sub_contab_size[np.arange(attribute_index, attributes.size,1,np.uint32)]

        rows_to_increment = contab_row + np.sum(useful_data * subs,1)
        
        #Irrittantingly I can't find a parallel way of doing this :( (that will actually compile)
        #as rows_to_increment may contain duplicates so can't just increment the slice...
        for row in rows_to_increment:
            contab[row] += 1

    else:
        num_to_modify = sub_contab_size[attribute_index]
        attribute = attributes[attribute_index]
        # Do the sub-table for the MCV first
        # Then we can calculate the others and subtract the value we get from
        # the appropriate value in the MCV table without worrying about
        # subtracting from a 0 value (we're using unsigned ints) or haveing to
        # calculate some total and subtract that from the MCV (would use extra
        # time/memory)
        vary_node_index = ad._ad_nodes[current_index]['children'] + attribute - parent_vnode_number
        #vary_node = self._tree_array[vary_node_index]
        MCV = ad._vary_nodes[vary_node_index]['mcv']
        MCV_contab_row = contab_row + (MCV*num_to_modify)
        _make_contab_internal(ad,attributes, attribute_index+1, current_index, parent_vnode_number, contab, MCV_contab_row, sub_contab_size)
        
        ad_node_index = ad._vary_nodes[vary_node_index]['children']
        for each_attribute_value in np.arange(ad._arities[attribute]):
            if each_attribute_value != MCV:
                _make_contab_internal(ad,attributes, attribute_index+1, ad_node_index, attribute+1, contab, contab_row, sub_contab_size)
                contab[MCV_contab_row:(MCV_contab_row+num_to_modify)] -= contab[contab_row:(contab_row+num_to_modify)]
                
                ad_node_index +=1 
            
            contab_row += num_to_modify       

@jit(nopython=True, parallel=True)
def _make_contab(ad, attributes):
    """
        Had to be separate from the class so that the parallel compilation option could be used
    """
    if attributes.size > ad._max_attributes_in_contab:
        raise ValueError("Number of attributes greater than maximum for contabs from this ADtree")
    sub_contab_size = np.empty(attributes.size, np.uint32)
    sub_contab_size[-1] = 1

    for i in range(attributes.size -2, -1, -1):
        sub_contab_size[i] = ad._arities[attributes[i+1]] * sub_contab_size[i+1]

    contab_height = ad._arities[attributes[0]] * sub_contab_size[0]
    contab_for_counts = np.zeros(contab_height, dtype=np.uint32)    
    
    _make_contab_internal(ad, attributes, 0, 0, 0, contab_for_counts, 0, sub_contab_size)
    return contab_for_counts