import numpy as np
from collections import defaultdict


def physi_constrain(a):
    result = []
    passset = ['x^2','y^2','z^2','xy','xz','yz','x^4','y^4','z^4']
    for i, row in enumerate(a):
        # Combine all other rows except current row
        other_rows = [item for j, other_row in enumerate(a) if j != i for item in other_row]
        # Keep the elements if they appear in other rows
        filtered_row = [elem for elem in row if elem in other_rows or elem in passset]
        result.append(filtered_row)

    return result

# select candidates based on the indicator
def selc_candifunc(candifunc, indicator):
    return [[func for func, ind in zip(row_func, row_ind) if ind] for row_func, row_ind in zip(candifunc, indicator)]
   
# build an indicator of a in candifunc            
def indicate(candifunc, a):
    return [[1 if item in res else 0 for item in row] for row, res in zip(candifunc, a)]
           
# find the indices pair of terms needed to be physics constrained  
def find_index(candifunc, elements2find):
    dic = defaultdict(list)
    
    for ind0, sublist in enumerate(candifunc):
        for ind1, item in enumerate(sublist):
            if item in elements2find:
                dic[item].append((ind0, ind1))
    return dict(dic)

# find the indices pair of terms needed to be physics constrained in a flattened candidate function list
def find_index_flatten(candifunc, index):
    indlist = []
    for ind in index:
        n = 0
        for i in range(ind[0]):
            n += len(candifunc[i])
        n = n + ind[0] + ind[1]
        indlist.append(n)

    return indlist


if __name__ == '__main__':

    candifunc = [['x', 'y', 'z', 'xy', 'yz', 'zx', 'y^2', 'z^2', 'x^3'],
                 ['x', 'y', 'z', 'xy', 'yz', 'zx', 'x^2', 'z^2', 'y^3'],
                 ['x', 'y', 'z', 'xy', 'yz', 'zx', 'x^2', 'y^2', 'z^3']]

    candifunc1 = [['x^2', 'xy', 'xz', 'x^2y', 'xyz', 'x^2z', 'xy^2', 'xz^2', 'x^4'],
                 ['xy', 'y^2', 'yz', 'xy^2', 'y^2z', 'xyz', 'x^2y', 'yz^2', 'y^4'],
                 ['xz', 'yz', 'z^2', 'xyz', 'yz^2', 'xz^2', 'x^2z', 'y^2z', 'z^4']]

    indicator = [[0, 0, 1, 0, 0, 0, 0, 1, 1],
     [0, 0, 1, 0, 1, 0, 1, 1, 1],
     [0, 1, 0, 0, 0, 0, 1, 1, 1]]

    selc_candifunc1 = selc_candifunc(candifunc1, indicator)

    phy_candifunc1 = physi_constrain(selc_candifunc1)

    indicator = indicate(candifunc1, phy_candifunc1)

    final_candifunc = selc_candifunc(candifunc, indicator)

    print(final_candifunc)