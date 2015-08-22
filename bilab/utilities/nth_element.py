
""" function partition, _select, select 
    for partial order
    Ref http://rosettacode.org/wiki/Quickselect_algorithm#Python
"""
import random
import numpy as np

__all__ = ["nth_element"]

def medianOf3(vector, a, b, c, comp):
    #comp a comparator, i.e. a boolean function accepting two parameters a and b,
    #        and returning true if a < b and false if a >= b.
    # in C or javascript
    # comp (A, B) ?
    #   comp (B, C) ? b : comp (A, C) ? c : a :
    #       comp (A, C) ? a : comp (B, C) ? c : b;
    print vector
    print("a={}, b={}, c={}".format(a,b,c))
    A = vector[a]
    B = vector[b]
    C = vector[c]
    if comp(A, B):
        # A < B
        if comp(B, C):
            # A < B < C
            return b
        else:
            # A<B, C<B
            if comp(A, C):
                # A < C < B
                return c
            else:
                # A<B, A>C, B>C -> B > A > C
                return a
    else:
        # A > B
        if comp(A, C):
            # C > A > B 
            return a
        else:
            # A>B and A> C
            if comp(B, C):
                #A>B B<C A>C -> A>C>B
                return c
            else:
                # A>B A>C B>C -> A>B>C
                return b

def partition(vector, left, right, pivotIndex, comp):
    pivotValue = vector[pivotIndex]
    # Move pivot to end
    vector[pivotIndex], vector[right] = vector[right], vector[pivotIndex]  
    storeIndex = left
    for i in range(left, right):
        #if vector[i] < pivotValue:
        if comp(vector[i], pivotValue):
            vector[storeIndex], vector[i] = vector[i], vector[storeIndex]
            storeIndex += 1
    # Move pivot to its final place
    vector[right], vector[storeIndex] = vector[storeIndex], vector[right]  
    return storeIndex
 
def _select(vector, left, nth, right, comp):
    "Returns the n-th smallest, (nth >= 0), element of vector within vector[left:right+1] inclusive."
    while True:
        # select pivotIndex between left and right
        #pivotIndex = random.randint(left, right)
        # meadian (left+rigth)>>1
        pivotIndex = medianOf3(vector, left, right, (left + right)>>1, comp)
        pivotNewIndex = partition(vector, left, right, pivotIndex, comp)
        #pivotDist = pivotNewIndex - left
        pivotDist = pivotNewIndex - left + 1
        if pivotDist == nth:
            return vector[pivotNewIndex]
        elif nth < pivotDist:
            right = pivotNewIndex - 1
        else:
            #nth -= pivotDist + 1
            #left = pivotNewIndex + 1
            nth -= pivotDist
            left = pivotNewIndex + 1

def nth_element(vector, left, nth, right, comp):
    """ Return the k-th smallest, (k >= 0), element of vector within vector[left:right+1].
    left, right default to (0, len(vector) - 1) if omitted
    """
    if left is None:
        left = 0
    lv1 = len(vector) - 1
    if right is None:
        right = lv1
    assert vector and nth >= 0, "Either null vector or k < 0 "
    assert 0 <= left <= lv1, "left is out of range"
    assert left <= right <= lv1, "right is out of range"
    return _select(vector, left, nth, right, comp)
