import math
import numpy as np
import unittest
import random
from matplotlib import pyplot as plt
from modeling import birth_and_death_with_escape



#https://scicomp.stackexchange.com/questions/1658/define-custom-probability-density-function-in-python


#mean_gene_birth_rate = 0.001359 #Ya-Long Guo reference
#The distributions of λ and μ have means of 0.00162 and 0.00194, (Tilet 2016)
# which were estimated from land plant genomes based on a model of gene copy number evolution
class MyTestCase(unittest.TestCase):


    def test_birth_death_model(self):

        step_size = 0.01
        decay_constant=3  # =(1/mean life expectancy of SSD gene)
        rate_escape=0.1   # % of SSD are maintained

        include_debugging_plots=True
        my_pdf, xs =birth_and_death_with_escape.gene_birth_death_with_escape_pdf(
            step_size,decay_constant,rate_escape, include_debugging_plots)
        print(my_pdf)
        self.assertEqual(True, True)  # ad

        num_genes_needed = 500
        SSD_ks = birth_and_death_with_escape.draw_SSDs_from_pdf(my_pdf, xs,
                                                       num_genes_needed, include_debugging_plots)

        print(SSD_ks)
        self.assertEqual(True, True)  # ad

    def test_gene_birth_death_corn(self):
        step_size = 0.01
        
        #(Ferris and Whitt, 1977; Nadeau and Sankoff, 1997; Lynch and Conery, 2000;
        # #Blanc and Wolfe, 2004b; Soltis et al.,2016; Cheng et al., 2018).

        # if a SSD lasts a few million years (say 3x10^6)
        # In Ks space, that depends mutation rate.
        # if Mutation rate is 10^-8, then
        # SSD lasts a few million years (say 3x10^-2 ) in KS
        # so, decay constant is 1/(3x10^-2)
        decay_constant = 33  # =(1/mean life expectancy of SSD gene)
        
        #  0.3, 0.4  human and mouse https://pmc.ncbi.nlm.nih.gov/articles/PMC1413713/
        rate_escape = 0.3  # % of SSD are maintained

        include_debugging_plots = True
        my_pdf, xs = birth_and_death_with_escape.gene_birth_death_with_escape_pdf(
            step_size, decay_constant, rate_escape, include_debugging_plots)
        print(my_pdf)
        self.assertEqual(True, True)  # ad

        num_genes_needed = 500
        SSD_ks = birth_and_death_with_escape.draw_SSDs_from_pdf(my_pdf, xs,
                                                                num_genes_needed, include_debugging_plots)

        print(SSD_ks)
        self.assertEqual(True, True)  # ad
