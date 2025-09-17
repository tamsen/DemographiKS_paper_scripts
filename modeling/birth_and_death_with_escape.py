import scipy.stats as st

#https://scicomp.stackexchange.com/questions/1658/define-custom-probability-density-function-in-python
class my_pdf(st.rv_continuous):
    def _pdf(self,x):
        return 3*x**2  # Normalized over its range, in this case [0,1]

my_cv = my_pdf(a=0, b=1, name='my_pdf')

# make an exponential distribution from 1 to KS =4, and sample from it
# + percent escape..