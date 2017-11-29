import numpy as np
from mypylib import *


def rotate_matrix(tei, transform_matrix):

    # Perform all contractions
    new_integral1 = np.einsum('ijkl,ia->ajkl', tei, transform_matrix)
    new_integral2 = np.einsum('ajkl,jb->abkl', new_integral1, transform_matrix)
    new_integral3 = np.einsum('abkl,kc->abcl', new_integral2, transform_matrix)
    new_integral4 = np.einsum('abcl,ld->abcd', new_integral3, transform_matrix)
    return new_integral4


integrals = np.arange(16).reshape(2, 2, 2, 2)
C = np.arange(1, 5).reshape(2, 2)

rotated_integrals = rotate_matrix(integrals, C)

