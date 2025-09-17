import math
import numpy as np
import unittest
import scipy.stats as st
from matplotlib import pyplot as plt

#random_numbers_array = np.random.exponential(scale=2.5, size=5)
#print(random_numbers_array)


#https://scicomp.stackexchange.com/questions/1658/define-custom-probability-density-function-in-python
class my_toy_pdf(st.rv_continuous):
    def _pdf(self,x):
        return 3*x**2  # Normalized over its range, in this case [0,1]

class my_birth_death_pdf(st.rv_continuous):
    def _pdf(self,x):
        gene_birth_death_rate=0.001359
        escape_rate=0.16#(0.00194-0.00162)/0.00194
        #return 100*math.exp(-1*(gene_birth_death_rate*x))
        #return 100 * math.exp(-1 * (0.2 * x))
        return max(0,0.1-(0.01*x))

def simple_fxn(x):
        return max(0,0.1-(0.01*x))

#mean_gene_birth_rate = 0.001359 #Ya-Long Guo reference
#The distributions of λ and μ have means of 0.00162 and 0.00194, (Tilet 2016)
# which were estimated from land plant genomes based on a model of gene copy number evolution
class MyTestCase(unittest.TestCase):
    def test_something(self):

        #my_cv = my_toy_pdf(a=0, b=100, name='my_pdf') #a & b = upper and lower bd of support
        my_cv = my_birth_death_pdf(a=0, b=10, name='my_pdf')

        print(my_cv)
        random_values = my_cv.rvs(size=100)
        random_values = [min(10,r)  for r in random_values]
        print(random_values)
        #bins = bins
        hist_ys_real, bins, patches = plt.hist(random_values, bins=40, facecolor='b', alpha=0.25,
                                                    density=True)



        width2=(bins[2]-bins[1])

        width=(bins[2]-bins[1])/2
        bar_plot_xs = [b + width for b in bins[0:len(bins) - 1]]  # to match bar-plot axes
        expected_pdf_ys = [simple_fxn(x) for x in (bins + [bins[39]+width2])]
        #plt.plot(bar_plot_xs,hist_ys_real, label="observed pdf")
        #plt.plot(bar_plot_xs,expected_pdf_ys, label="expected pdf")
        #plt.plot(bar_plot_xs,expected_pdf_ys, label="expected pdf")
        plt.plot(bins, expected_pdf_ys, label="expected pdf", color='k')


        out_png= "/home/tamsen/Data/DemographiKS_output_from_mesx/bd_process.png"
        plt.savefig(out_png)

        self.assertEqual(True, False)  # add assertion here


    def test_exp_something(self):
        xs = np.arange(0, 4, 0.01)
        print(xs)
        out_png= "/home/tamsen/Data/DemographiKS_output_from_mesx/bd_exp_process.png"
        decay_constant=3
        rate_escape=0.1
        exp_decay_ys=[math.exp(-1*decay_constant*x) for x in xs]
        ys_that_have_died=[1.0-math.exp(-1*decay_constant*x) for x in xs]
        ys_that_escape=[y*rate_escape for y in ys_that_have_died]
        ys_remaining=[exp_decay_ys[idx]+ys_that_escape[idx] for idx in range(0,len(exp_decay_ys))]
        plt.plot(xs, exp_decay_ys, label="SSDs not dead yet", color='k')
        plt.plot(xs, ys_that_have_died, label="SSDs should be dead", color='b')
        plt.plot(xs, ys_that_escape, label="SSDs that escape", color='c')
        plt.plot(xs, ys_remaining, label="ys_remaining", color='r')

        #https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution
        cdf=

        plt.legend()
        plt.savefig(out_png)

        self.assertEqual(True, False)  # add assertion here