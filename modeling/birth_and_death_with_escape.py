import math

import numpy as np


def gene_birth_death_with_escape_pdf():
    
    step_size = 0.01
    decay_constant=3  # =(1/mean life expectancy of SSD gene)
    rate_escape=0.1   # % of SSD are maintained
    xs = np.arange(0, 4, step_size)

    exp_decay_ys=[math.exp(-1*decay_constant*x) for x in xs]
    ys_that_have_died=[1.0-math.exp(-1*decay_constant*x) for x in xs]
    ys_that_escape=[y*rate_escape for y in ys_that_have_died]
    ys_remaining=[exp_decay_ys[idx]+ys_that_escape[idx] for idx in range(0,len(exp_decay_ys))]

    normalization_factor=1.0/(step_size*sum(ys_remaining))
    normalized_pdf=[normalization_factor*y for y in ys_remaining]
    return normalized_pdf