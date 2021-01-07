#!/usr/bin/env python
# coding: utf-8
# This code was extracted from the neat statespace equaitions

# In[2]:


import sympy as sy
from sympy.integrals import inverse_laplace_transform as ilt
from sympy.integrals import laplace_transform as lt
from sympy import pretty_print as pp

    
def display_steps(expr, string_msg = "", is_ipython = False):
    if string_msg != "":
        pp(string_msg)
    if is_ipython:
        display(expr)
    else:
        pp(expr)


# In[37]:


def state_transition_matrix(a, symbol = sy.symbols('t'), show_steps = True):
    '''
        Computes the state transition matrix for continous time systems.
        state_trans_matrix = e^(A*t)
        
        Params:
            A_matrix (list or sympy Matrix): A list of lists for the rows of the A matrix of the system
            symbol (sympy.Symbol): the symbol used to represent the time variable, defaults to 't'
            show_steps (bool, True): Pretty prints the steps of the computation

        Returns:
            sympy Matrix: the state transition matrix

    '''
    a = sy.Matrix(a)
    state_trans = sy.exp(a*symbol)
    if show_steps:
        display_steps("the state trans matrix is")
        display_steps(state_trans)
    return state_trans


def state_homogenous_response(A_matrix, initial_condition_col,symbol = None, discrete = False, show_steps = True):
    '''
    Computes the homogenous response of the state variables at time t (k for discrete).
    x_homo = state_trans_matrix * initial_cond
    
    Params:
        A_matrix (list or sympy Matrix): A list of lists for the rows of the A matrix of the system
        initial_condition_col (list or sympy Matrix): A list of lists for the inital conditions column
        discrete (bool, optional): is the system discrete time?
        symbol (sympy.Symbol): the symbol used to represent the time variable, defaults to 't' for contionous,
        and k for discrete
        show_steps (bool, True): Pretty prints the steps of the computation
    
    Returns:
        sympy Matrix: the homogenous response
     '''
    A_matrix = sy.Matrix(A_matrix)
    initial_condition_col = sy.Matrix(initial_condition_col)
    if discrete:
        symbol = sy.symbols('k')
        state_trans_mat = state_transition_matrix_discrete(A_matrix, symbol = symbol, show_steps = show_steps)
    else:
        symbol = sy.symbols('t')
        state_trans_mat = state_transition_matrix(A_matrix, symbol=symbol, show_steps=show_steps)
        
    x_homo =  state_trans_mat * initial_condition_col
    if show_steps:
        print("Finding the homogenous response")
        display_steps(A_matrix, "A matrix")
        display_steps(initial_condition_col, "Initial condition")
        display_steps(x_homo, "x homogenous response at time t is ")
    return x_homo


def state_forced_response(a, b, u, show_steps = True):
    '''Computes the forced response of the state variables based on the input u.
    x_homo = state_trans_matrix * initial_cond
    
    Params:
        a (list or sympy Matrix): A list of lists for the rows of the A matrix of the system
        b (list or sympy Matrix): A list of lists for the b matrix
        u (list,sympy Matrix, or a constant): the input to the system
        show_steps (bool, True): Pretty prints the steps of the computation
    
    Returns:
        sympy Matrix: the forced response
    '''
    a = sy.Matrix(a)
    b = sy.Matrix(b)
    # T is for tau
    t, T = sy.symbols('t T')
    state_trans_neg = state_transition_matrix(-1*a, T, show_steps = False)
    state_trans = state_transition_matrix(a, t, show_steps = False)
    # This is in case the u param is dependant on t, if it is we need to make it
    # into a T for the integration below
    try:
        u = u.replace(t,T)
    except:
        pass
    unevaluated_conv_integral = sy.Integral(state_trans_neg * b*u ,(T,0,t))
    conv_integral = unevaluated_conv_integral.doit() 
    forced_resp = state_trans * conv_integral
    if show_steps:
        display_steps("Finding the forced response")
        display_steps("A matrix")
        display_steps(a)
        display_steps("B matrix")
        display_steps(b)
        display_steps(state_trans,"State transition matrix: ")
        display_steps(state_trans_neg, "Negative state trans matrix")
        display_steps(unevaluated_conv_integral, "The convolution integral")
        display_steps(conv_integral, "Evaluated convolution integral")
        display_steps(forced_resp, "Foreced response is ")
        
    return forced_resp

def output_response(a,b,c,d,u, initial_condition_col, discrete = False, show_steps=True):
    '''Computes the output response of the state space model.
    y = C x homogenous_response + C x forced_response + Du
    
    Params:
        a (list or sympy Matrix): A list of lists for the rows of the A matrix of the system
        b (list or sympy Matrix): A list of lists for the b matrix
        c (list or sympy Matrix): A list of lists for the b matrix
        d (list or sympy Matrix): A list of lists for the b matrix
        u (list,sympy Matrix, or a constant): the input to the system
        initial_condition_col (list or sympy Matrix): A list of lists for the inital conditions column
        discrete (bool, optional): is the system discrete time?
        show_steps (bool, optional): Pretty prints the steps of the computation
    
    Returns:
        sympy Matrix: the output response of state space model
    '''
        
    a = sy.Matrix(a)
    b = sy.Matrix(b)
    c = sy.Matrix(c)
    d = sy.Matrix(d)
    initial_condition_col = sy.Matrix(initial_condition_col)
    # set show steps to false because it will be redundant
    homo_resp = state_homogenous_response(a,initial_condition_col,discrete = discrete, show_steps = False)
    forced_resp = []
    
    if discrete:
        forced_resp = state_forced_response_discrete(a, b,u, show_steps = False)
    else:
        forced_resp = state_forced_response(a, b,u, show_steps = False)

    full_resp_times_c_col = c * forced_resp + c * homo_resp
    
    y =  full_resp_times_c_col + d*u
    if show_steps:
        print("finding the output response")
        display_steps(a,"A matrix")
        display_steps(b,"B matrix")
        display_steps(c,"C matrix")
        display_steps(d,"D matrix")
        display_steps(initial_condition_col, "initial condition")
        display_steps(homo_resp, "homogenous resp")
        display_steps(forced_resp,"forced resp")
        display_steps(homo_resp + forced_resp, "Full x response ie x homogenous + x forced")
        display_steps(full_resp_times_c_col, "Full response times c col at time t (or k if discrete)")
        display_steps(y, "the output y is ")
    
    return y


# In[23]:


def state_transition_matrix_discrete(a, symbol = sy.symbols('k'), show_steps = True):
    '''
    Computes the state transition matrix for continous time systems.
    state_trans_matrix = e^(A*t)

    Params:
        A_matrix (list or sympy Matrix): A list of lists for the rows of the A matrix of the system
        symbol (sympy.Symbol): the symbol used to represent the time variable, defaults to 't'
        show_steps (bool, True): Pretty prints the steps of the computation

    Returns:
        sympy Matrix: the state transition matrix

    '''
    a = sy.Matrix(a)
    state_trans = a ** symbol
    if show_steps:
        display_steps(state_trans, "the state trans matrix is")
    return state_trans

def state_forced_response_discrete(a, b, u, show_steps = True):
    '''Computes the forced response of the state variables based on the input u.
    x_homo = state_trans_matrix * initial_cond
    
    Params:
        a (list or sympy Matrix): A list of lists for the rows of the A matrix of the system
        b (list or sympy Matrix): A list of lists for the b matrix
        u (list,sympy Matrix, or a constant): the input to the system
        show_steps (bool, True): Pretty prints the steps of the computation
    
    Returns:
        sympy Matrix: the forced response
    '''
    a = sy.Matrix(a)
    b = sy.Matrix(b)
    # T is for tau
    k, i = sy.symbols('k i')
    state_trans_neg = state_transition_matrix_discrete(a, k-1-i, show_steps = False)
    state_trans = state_transition_matrix_discrete(a, k, show_steps = False)
    # This is in case the u param is dependant on k, if it is we need to make it
    # into an i for the summation below
    try:
        u = u.replace(k,i)
    except:
        pass
    unevaluated_conv_sum = sy.Sum(state_trans_neg * b * u ,(i, 0, k-1))
    forced_resp = unevaluated_conv_sum.doit()
    # sympy returns None if the answer is zero
    if not forced_resp:
        forced_resp = sy.zeros(state_trans.shape[0])
    if show_steps:
        print("Finding the forced response")
        display_steps(a, "A matrix")
        display_steps(b, "B matrix")
        display_steps(state_trans,"State transition matrix: ")
        display_steps(state_trans_neg,"Negative state trans matrix")
        display_steps(unevaluated_conv_sum,"The convolution summation")
        display_steps(forced_resp ,"Forced response is ")
        
    return forced_resp


# In[ ]:





# In[29]:


## Testing

# a = sy.Matrix([
#     [0,1],[-2,-3]
# ])
# # state_trans_matrix_inverse_laplace_method(a)
# print("with cayley hamilton")
# state_trans_matrix(a)
# # display("We can also do this for an equivalent answer", sy.exp(a*t))
# # state_var_homo_response(a,[0,1])

# k = sy.symbols('k')
# a = sy.Matrix([
#     [1,0],[k+1, 0]
# ])
# b = [[0],[0]]
# c = [[1, 1/(k+1)]]
# d = [1]
# x = [[1],[1]]
# u = 1
# state_homogenous_response(a,x, discrete=True)


# In[32]:


# k = sy.symbols('k')
# a = sy.Matrix([
#     [1,0],[k+1, 0]
# ])
# b = [[0],[0]]
# c = [[1, 1/(k+1)]]
# d = [1]
# x = [[1],[1]]
# u = 1
# # output_response(a,b,c,d,u,x, discrete=True)
# # state_homogenous_response(a,x,discrete=True)
# sy.Pow(a,k)


# In[17]:


# a= sy.Matrix([
# [sy.Rational(-3,4) , sy.Rational(-1,4)],
#     [sy.Rational(-1,4) , sy.Rational(-3,4)]
# ])
# b = [
#     [1],
#     [0]
# ]
# u=1
# c = [[1,0]]
# d = [4]
# ic = [[0],
#      [0]]
# # matrix_info(a)
# # state_trans = state_trans_matrix(a)
# # state_trans_cayley_hamilton_method(a)
# # state_var_forced_response_discrete(a, b, u)
# a = sy.Matrix([
#     [0,1],
#     [sy.Rational(-1,4), -1]
# ])

# state_transition_matrix_discrete(a)


# In[ ]:




