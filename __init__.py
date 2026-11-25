from src.Inverting the linking form import *
from sympy import *

a = Rational(1,5)
b = Rational(1,7)
Q1 = Matrix([[a,a],[a,a**2]])
Q2 = Matrix([[b,b],[b,b**2]])
Q = diag(Q1,Q2)
gen_list = [5,25,7,49]
index_list = prime_factor_decomp(gen_list)[2]
B = almost_diagonalize_form(gen_list , Q, index_list)
c = Rational(1,3)

gen_list = [3,15]
Q = Matrix([[c,c],[c,Rational(1,15)]])
print(full_process(gen_list , Q))
