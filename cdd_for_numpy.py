# from scipy.spatial import ConvexHull as qhull
import numpy as np
import itertools
import cdd #Install via "pip install pycddlib"



def toHrep(mat, input='Vertices'):
    if input == 'Vertices':
        padding = np.ones((len(mat), 1), dtype=int)
    elif input == 'Rays':
        padding = np.zeros((len(mat), 1), dtype=int)
    newmat = np.hstack((padding, mat))
    cddmat_in = cdd.Matrix(newmat, number_type='fraction')
    cddmat_in.rep_type = cdd.RepType.GENERATOR
    cddmat_out = cdd.Polyhedron(cddmat_in).get_inequalities()

    nparray_out = np.empty((cddmat_out.row_size, cddmat_out.col_size), dtype=float)
    for i in range(cddmat_out.row_size):
        row = cddmat_out.__getitem__(i)
        to_multiply = set()
        for j, element in enumerate(row):
            if type(element)!=int:
                to_multiply.add(element.denominator)
                nparray_out[i,j] = element.numerator/element.denominator
            else:
                nparray_out[i, j] = element
        nparray_out[i] = nparray_out[i] * np.prod(list(to_multiply))

    eq_indices = cddmat_out.lin_set
    ineq_indices = frozenset(set(range(cddmat_out.row_size)).difference(eq_indices))
    nparray_out = nparray_out.astype(int)
    inequalities_array = nparray_out[list(ineq_indices)]
    equalities_array = nparray_out[list(eq_indices)]

    inequalities_b_part = inequalities_array[:, :1]
    inequalities_A_part = inequalities_array[:, 1:]
    equalities_b_part = equalities_array[:, :1]
    equalities_A_part = equalities_array[:, 1:]
    final_b = np.vstack((inequalities_b_part, equalities_b_part, equalities_b_part))
    final_A = np.vstack((inequalities_A_part, equalities_A_part, -equalities_A_part))
    if input=='Vertices':
        return np.hstack((final_b, final_A))
    elif input=='Rays':
        return final_A



if __name__ == '__main__':
    test = np.asarray(
           [[1, 0, 1],
            [0, 1, 1],
            [-1, 0, 1],
            [0, -1, 1]])

    test2 = [[3,7], [10,7], [8,2], [5,3]]
    #
    # test3 = np.asarray(
    #        [[1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
    #         [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0],
    #         [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1],
    #         [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
    #         [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
    #         [0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
    #         [0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
    #         [1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0],
    #         [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
    #         [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
    #         [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
    #         [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
    #         [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1]
    #         ])

    test3 = np.asarray(
           [[1, 1, 0, 0],
            [0, 0, 1, 1],
            [1, 0, 0, 1],
            [0, 1, 1, 0]])
    print(toHrep(test, input='Vertices'))
    print(toHrep(test, input='Rays'))
    print(toHrep(test2, input='Vertices'))
    print(toHrep(test2, input='Rays'))
    print(toHrep(test3, input='Vertices'))
    print(toHrep(test3, input='Rays'))
