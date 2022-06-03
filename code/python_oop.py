# Load dependencies
import pandas as pd
import numpy as np
from scipy.stats import beta, dirichlet
import time

# Define class for utility object

### --------------
### UTILITY CLASS
### --------------

class Utility:

    sixtypes = ["ABC", "ACB", "BAC", "BCA", "CAB", "CBA"]
    threetypes = ["A", "B", "C"]

    def translate_pref(self, v):
        '''Returns vector of preference attributes for every voter'''

        conditions = [
            v['A'].ge(v['B']) & v['B'].ge(v['C']),
            v['A'].ge(v['C']) & v['C'].ge(v['B']),
            v['B'].ge(v['A']) & v['A'].ge(v['C']),
            v['B'].ge(v['C']) & v['C'].ge(v['A']),
            v['C'].ge(v['A']) & v['A'].ge(v['B']),
            v['C'].ge(v['B']) & v['B'].ge(v['A'])
        ]

        return pd.DataFrame(
            np.select(conditions, self.sixtypes)
        )

    def types_prop(self, p, system):
        '''Calculate share of each type depending on system'''

        pvec = []

        # Set up variables conditional on system
        if system == "rcv":
            iter = self.sixtypes
            pref = self.pref_rcv
        if system == "plur":
            iter = self.threetypes
            pref = self.pref_plur

        # Store share for every voter type
        for p in iter:
            share = np.mean(pref == p)
            pvec.append(share)

        return np.array(pvec)

    def __init__(self, U):

        if not isinstance(U, pd.DataFrame):
            raise Error('U is not a data frame')

        if not len(U.columns) == 3:
            raise Error('U does not have exactly 3 cols')

        # Set up initial properties
        U.columns = ['A', 'B', 'C']
        self.U = U,
        self.N = len(U.index)

        # Infer preference types
        self.pref_rcv = self.translate_pref(U)
        self.pref_plur = [x[0] for x in self.pref_rcv[0]]

        # Calculate type shares
        self.prop_rcv = self.types_prop(self.pref_rcv, "rcv")
        self.prop_plur = self.types_prop(self.pref_plur, "plur")


# class Iteration:

#     def create_pmat(U, exvec):

#     def __init__(self, U, system, exvec=None, s=85):

#         # Set up expected vector
#         if exvec == None and system == "rcv":
#             exvec = U.prop_rcv
#         if exvec == None and system == "plur":
#             exvec = U.prop_plur

#         self.exvec = exvec
#         self.s = s

### -------------------------
### PIVOT PROB FUNCTIONS
### -------------------------

# Once this part is done, I want to shove all of this into a P-mat class

def get_increment_midpoints(start, end, increments):
    '''Return the midpoints of increments along line'''
    # TODO: change this into numpy.linspace(start, stop, step)

    return np.linspace(start, end, increments, endpoint = False)


    inc_length = (end - start) / increments
    out = []

    point = start + inc_length / 2
    out.append(point)

    # Calculate points sequentially moving from start
    while point <= end - inc_length / 2:
        point = point + inc_length
        out.append(point)
        if len(out) > (increments):
            break

    return np.array(out)


def en_second_pp(alpha, increments=20):

    if len(alpha) != 6:
        print('Alpha length incorrect!')

    # Set up midpoints and integration limits
    midpoints_y = get_increment_midpoints(1/4, 1/2, increments)


    int_limits = 2 * midpoints_y - 0.5
    int_limits[midpoints_y > 1/3] = midpoints_y[midpoints_y > 1/3] / 2
    width_of_channel = np.sqrt(6)/4
    distance_between_midpoints = np.sqrt(2) * (midpoints_y[1] - midpoints_y[0])

    beta_parts = beta.cdf(2 * int_limits, alpha[5], sum([alpha[2], alpha[3]]))

    # Set up coordinates for integration
    y_mat = [
        midpoints_y,
        1/2 - midpoints_y,
        np.array([1/2] * len(midpoints_y))
    ]

    # Alpha for Dirichlet f'n
    alpha_term = [sum([alpha[0], alpha[1]]), alpha[4],
                  sum([alpha[2], alpha[3], alpha[5]])]

    # Calculate density at each increment
    dir_parts = []
    for z in range(0, len(midpoints_y)):
        d = dirichlet.pdf(
            x=[y_mat[0][z], y_mat[1][z], y_mat[2][z]],
            alpha=alpha_term
        )
        dir_parts.append(d)

    dir_parts = np.array(dir_parts)/np.sqrt(3)

    return width_of_channel * distance_between_midpoints * sum(dir_parts * beta_parts)



# Somehow this doesn't give the right outcome despite the formula being correct


def en_first_pp(alpha, increments=100):
    '''Return first-round pivot probabilities'''

    # Set up first-round vote shares
    fp_alpha = [
        sum([alpha[0], alpha[1]]),
        sum([alpha[2], alpha[3]]),
        sum([alpha[4], alpha[5]])
    ]

    # print(fp_alpha)

    midpoints_y = get_increment_midpoints(1/4, 1/3, increments)
    pr_tie_second_y      = [pd.NA] * len(midpoints_y)
    pr_1_beats_3_given_y = [pd.NA] * len(midpoints_y)
    pr_2_beats_3_given_y = [pd.NA] * len(midpoints_y)
    x1 = [pd.NA] * len(midpoints_y)
    x2 = [pd.NA] * len(midpoints_y)
    x3 = [pd.NA] * len(midpoints_y)

    # Calculate densities for all midpoints
    for z in range(0, len(midpoints_y)):
        y = midpoints_y[z]
        k = (
            dirichlet.pdf(
                [y, y, 1 - 2*y],
                fp_alpha
            ) / np.sqrt(3))
        pr_tie_second_y[z] = k
        pr_1_beats_3_given_y[z] = beta.cdf(2 - 1 / (2 * y), alpha[3], alpha[2])
        pr_2_beats_3_given_y[z] = beta.cdf(2 - 1 / (2 * y), alpha[1], alpha[0])

        x1[z] = pr_tie_second_y[z] * pr_1_beats_3_given_y[z] * pr_2_beats_3_given_y[z]
        x2[z] = pr_tie_second_y[z] * (1 - pr_1_beats_3_given_y[z]) * pr_2_beats_3_given_y[z]
        x3[z] = pr_tie_second_y[z] * pr_1_beats_3_given_y[z] * (1 - pr_2_beats_3_given_y[z])

    distance_between_midpoints = np.sqrt(6)*(midpoints_y[2] - midpoints_y[1])
    width_of_channel = 1 / np.sqrt(2)
    fact = distance_between_midpoints * width_of_channel

    ij_ij = fact * sum(x1)
    ij_kj = fact * sum(x2)
    ij_ik = fact * sum(x3)

    return [ij_ij, ij_kj, ij_ik]

def en_pivprob(alpha, increments = 20):
    '''Return Pr(pivotal event)'''

    first = en_first_pp(alpha, increments = increments)
    second = en_second_pp(alpha, increments = increments)

    return [first, second]


def en_all_pivprob(vvec, s, increments = 20):
    '''Calculate all of the pivot probabilities (12 in total)'''

    # Do this for second round
    a1 = [v * s for v in vvec]
    o2 = [1, 0, 4, 5, 2, 3]
    a2 = [v * s for v in [vvec[i] for i in o2]]
    o3 = [3, 2, 5, 4, 0, 1]
    a3 = [v * s for v in [vvec[i] for i in o3]]

    ab = en_second_pp(a1, increments = increments)
    ac = en_second_pp(a2, increments = increments)
    bc = en_second_pp(a3, increments = increments)

    ab2 = en_first_pp(a1, increments = increments)
    ac2 = en_first_pp(a2, increments = increments)
    bc2 = en_first_pp(a3, increments = increments)

    return [ab, ac, bc, ab2, ac2, bc2]

### NEXT STEPS:
##      - construct Pmat
##      - construct EU frame


### ----------------------------
### UNIT TESTS
### ----------------------------

# Set up toy utility frame
test_data = {
    'A': [4, 7, 8, 9, 2, 4],
    'B': [3, 4, 9, 3, 1, 0],
    'C': [1, 8, 6, 4, 5, 6]
}

test_data = pd.DataFrame(test_data)
test_vec = [0.1, 0.1, 0.3, 0.1, 0.1, 0.3]

# Check Utility class behaviour
test_class = Utility(test_data)
test_class.N

# Check if second-round pivots are correct
en_second_pp(test_vec)

# Check if first-round pivots are correct
en_first_pp(test_vec)

# Check if wrapper function works
en_pivprob(test_vec)

# Check if outer wrapper works

start_time = time.time()
en_all_pivprob(test_vec, 85)
print("--- %s seconds ---" % (time.time() - start_time))

en_all_pivprob(test_vec, 85)
