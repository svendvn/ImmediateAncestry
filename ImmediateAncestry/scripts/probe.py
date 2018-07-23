from ggrandparents_model import find_smallet_equivalence_class
from math import factorial
from itertools import combinations, combinations_with_replacement
import sys, math
from time import time, localtime, strftime
from itertools import product
import re 
from datetime import datetime



@cython.boundscheck(False)
def backward2(sequence,
             constants,
             np.ndarray[DTYPE_t, ndim=1] initial,
             np.ndarray[DTYPE_t, ndim=2] transitions,
             np.ndarray[DTYPE_t, ndim=2] emissions,
             char_map, label_map, name_map):
    """
    Compute the probability of being in some state and generating the rest of
    the sequence.
    This function implements the scaled backward algorithm.
    :param sequence str: a string over the alphabet specified by the model.
    :param constants np.ndarray: an array of the constants used to normalize
                                 the forward table.
    :rtype: np.ndarray
    :return: the scaled backward table.
    """

    sequence = sequence.upper()

    cdef int no_observations = len(sequence)
    cdef int no_states = len(initial)

    cdef np.ndarray[DTYPE_t, ndim=2] M = \
        np.zeros([no_observations, no_states], dtype=DTYPE)

    cdef unsigned int i, j, k, observation
    cdef double prob, state_sum

    M[no_observations-1] = 1.0 / constants[no_observations - 1]

    for i in range(no_observations-2, -1, -1):
        observation = char_map[sequence[i]]
        for j in range(no_states):
            state_sum = 0.0
            for k in range(no_states):
                state_sum += M[(i + 1), k] * transitions[j, k]
            M[i, j] = state_sum * emissions[j, observation]
        M[i] = M[i] / constants[i]

    return M
