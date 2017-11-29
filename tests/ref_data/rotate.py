import numpy as np
from mypylib import *


def rotate_matrix(tei, transform_matrix):

    # Perform all contractions
    new_integral = np.einsum('ijkl,ia->ajkl', tei, transform_matrix)
    new_integral = np.einsum('ajkl,jb->abkl', new_integral, transform_matrix)
    new_integral = np.einsum('abkl,kc->abcl', new_integral, transform_matrix)
    new_integral = np.einsum('abcl,ld->abcd', new_integral, transform_matrix)

    return new_integral


integrals = np.arange(16).reshape(2, 2, 2, 2)
C = np.arange(1, 5).reshape(2, 2)
print(C)

rotated_integrals = rotate_matrix(integrals, C)
to_file(rotated_integrals, "rotated1.data")
