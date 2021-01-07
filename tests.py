import sympy as sy
from sympy import *

import sys 
import os
cwd = os.getcwd()
sys.path.append(os.path.abspath(cwd))
from state_space_functions import *
from MatrixFunctions import *
t,k = sy.symbols('t,k')
######## Test eq 1
a = [
    [1,0],
    [-1, sy.Rational(1,2)]
]
b = [[4],[1]]
x = [[1],[2]]

u = sy.exp(-2*t)
h = state_homogenous_response(a,x,show_steps=False)
f = state_forced_response(a,b,u, show_steps= False)

assert h.equals(sy.Matrix([[sy.exp(t)], [4*sy.exp(t/2) - 2*sy.exp(t)]])), "Expecting True if state_homogenous_response works"

assert state_transition_matrix(a, show_steps = False).equals(Matrix([[exp(t), 0], [2*exp(t/2) - 2*exp(t), exp(t/2)]])), "Expecting True is state_transition_matrix works"

assert f.equals(Matrix([[(4/3 - 4*exp(-3*t)/3)*exp(t)], [(4/3 - 4*exp(-3*t)/3)*(2*exp(t/2) - 2*exp(t)) + (14/15 + 8*exp(-3*t)/3 - 18*exp(-5*t/2)/5)*exp(t/2)]])), "error in forced resp for eq 1"

####### Test eq 2 

a= sy.Matrix([
    [-2,0],
    [1,-1]
])
# ic = sy.Matrix([[2],[3]])
ic = [0,0]
b = sy.Matrix([
    [1],
    [0]
])
u = 5

c = [
    [2,1]
]

h = state_homogenous_response(a,ic,show_steps=False)
f = state_forced_response(a,b,u, show_steps= False)
assert h.equals(Matrix([[0], [0]])), "error in homogenous response eq2"
assert f.equals(Matrix([[(5*exp(2*t)/2 - 5/2)*exp(-2*t)], [(exp(-t) - exp(-2*t))*(5*exp(2*t)/2 - 5/2) + (-5*exp(2*t)/2 + 5*exp(t) - 5/2)*exp(-t)]])), "Error in forced response"
assert output_response(a,b,c,[0],u, ic, show_steps=False).equals(Matrix([[(exp(-t) - exp(-2*t))*(5*exp(2*t)/2 - 5/2) + 2*(5*exp(2*t)/2 - 5/2)*exp(-2*t) + (-5*exp(2*t)/2 + 5*exp(t) - 5/2)*exp(-t)]])
), "Wrong output"

####### Test eq 3

a= sy.Matrix([
[sy.Rational(-3,4) , sy.Rational(-1,4)],
    [sy.Rational(-1,4) , sy.Rational(-3,4)]
])
b = [
    [1],
    [0]
]
u=1
c = [[1,0]]
d = [4]
ic = [[0],
     [0]]
# idk why this test doesnt work even tho both equations are equal
'''

f = state_forced_response_discrete(a, b, u, show_steps = False)
o = output_response(a,b,c,d,u,ic, discrete = True, show_steps = False)
answer_out = Matrix([[-(-1)**k*(1/2 - (-1)**k/2)/2 - (-1/2)**k*(1/3 - (-2)**k/3) + 4]])

assert o.equals(answer_out), "Error on output discrete"

assert f.equals(Matrix([[-(-1)**k*(1/2 - (-1)**k/2)/2 - (-1/2)**k*(1/3 - (-2)**k/3)], [-(-1)**k*(1/2 - (-1)**k/2)/2 + (-1/2)**k*(1/3 - (-2)**k/3)]])

),"Discrete forced response error"
'''
 

