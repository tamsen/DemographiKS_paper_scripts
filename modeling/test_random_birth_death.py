import math
import numpy as np
import unittest
import random
import scipy.stats as st
from matplotlib import pyplot as plt

#random_numbers_array = np.random.exponential(scale=2.5, size=5)
#print(random_numbers_array)


#https://scicomp.stackexchange.com/questions/1658/define-custom-probability-density-function-in-python


#mean_gene_birth_rate = 0.001359 #Ya-Long Guo reference
#The distributions of λ and μ have means of 0.00162 and 0.00194, (Tilet 2016)
# which were estimated from land plant genomes based on a model of gene copy number evolution
class MyTestCase(unittest.TestCase):

    def test_exp_something(self):
        step_size=0.01
        xs = np.arange(0, 4, step_size)
        print(xs)
        out_png= "/home/tamsen/Data/DemographiKS_output_from_mesx/bd_exp_process.png"
        decay_constant=3
        rate_escape=0.1
        arctan_ys = [(2/math.pi)*math.atan(4*math.pi*x) for x in xs]
        tanh_ys = [math.tanh(2*x) for x in xs]
        exp_decay_ys=[math.exp(-1*decay_constant*x) for x in xs]
        ys_that_have_died=[1.0-math.exp(-1*decay_constant*x) for x in xs]
        ys_that_escape=[y*rate_escape for y in ys_that_have_died]
        ys_remaining=[exp_decay_ys[idx]+ys_that_escape[idx] for idx in range(0,len(exp_decay_ys))]
        plt.plot(xs, exp_decay_ys, label="SSDs not yet decayed", color='k')
        plt.plot(xs, ys_that_have_died, label="SSDs should be dead", color='b')
        plt.plot(xs, ys_that_escape, label="SSDs that escape", color='c')
        plt.plot(xs, ys_remaining, label="ys_remaining", color='r')
        #plt.plot(xs, arctan_ys, label="arctan_ys", color='r')
        #plt.plot(xs, tanh_ys, label="tanh_ys", color='k')
        #https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution

        CDF=[step_size*sum(ys_remaining[0:i]) for i in range(0,len(xs))]
        print(CDF)
        print("last CDF=" + str(CDF[-1]))
        plt.plot(xs, CDF, label="ys_remaining CDF", color='pink')
        num_SSD_genes_remaining=step_size*sum(ys_remaining)
        print("num_SSD_genes_remaining=" + str(num_SSD_genes_remaining))
        plt.legend()
        plt.savefig(out_png)
        plt.close()
        normalization_factor=1.0/num_SSD_genes_remaining
        normalized_pdf=[normalization_factor*y for y in ys_remaining]
        population = xs
        probabilities = normalized_pdf
        multiple_choices = random.choices(population, weights=probabilities, k=500)

        plt.plot(xs, ys_remaining, label="prediction", color='r')
        hist_ys_real, bins, patches = plt.hist(multiple_choices, bins=40, facecolor='b', alpha=0.25,
                                                    density=True, label="observation")
        out_hist_png = "/home/tamsen/Data/DemographiKS_output_from_mesx/bd_exp_hist.png"
        plt.savefig(out_hist_png)
        plt.close()
        self.assertEqual(True, False)  # add assertion here