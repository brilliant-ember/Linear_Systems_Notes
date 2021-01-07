#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sympy as sy


# In[1]:


# Matrix Functions

def round_expr(expr, num_digits):
    # https://www.thetopsites.net/article/54462201.shtml
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(sy.Number)})

def matrix_info(matrix):
    ''' displays the determenant, inverse and charactersitic funtion of a matrix'''
    sMatrix = sy.Matrix(matrix)
    # char func is det(A- (Lambda x I)) = 0
    L = sy.symbols('L')
    LI = L * sy.eye(sMatrix.shape[0])
    charMatrix = sMatrix - LI
    print("A - (Lamda x I): ")
    display(charMatrix)
    inverse = get_inverse(sMatrix)
    display({
        "det": sMatrix.det(),
        "inverse":inverse,
        "charFunction": sy.latex(charMatrix.det())
    })
    
def get_eign(matrix):
    '''displays eignvalues and vectors of a matrix'''
    sMatrix = sy.Matrix(matrix)
    sym_eignvals = list(sMatrix.eigenvals().keys())
    sym_eignvects = []
    for tup in sMatrix.eigenvects():
        for v in tup[2]:
            sym_eignvects.append(list(v))
    display("Sympy eigenvalues: ", sym_eignvals)
    display("Sympy eigenvectors: ", sym_eignvects)

def get_inverse(m):
    '''
    Takes in a matrix t and returns the inverse assuming it exists, the inverse will have it's poles
    seperated by a gaussian dividtion by the determinent
    Params:
        m (list or sympy.Matrix)
    Returns:
        (sympy.Matrix) Inverted matrix or None if matrix doesnt exist
        
    '''
    m = sy.Matrix(m)
    if sy.det(m)!=0:
        m = m.inv()
        m = m.applyfunc(lambda e: sy.factor(e, gaussian=True))
    else:
        print("Not invertable, det=0")
        return None
    return m
A = sy.Matrix([
    [4, 0, 1],
    [2, 3, 2],
    [1, 0, 4]
])
def get_poles_from_matrix(a, show_steps = False):
    '''Takes the A matrix of a system and returns the poles, format is {pole: multiplicity} can be 	 	used to determine stability'''
    ss=sy.symbols('s')
    char_mtx = ss*sy.eye(a.rows) - a
    char_eq = char_mtx.det()
    if show_steps:
        display(char_mtx)
        display(char_eq)
    return(sy.roots(char_eq))


# In[11]:


# No unit testing done here, will do later.

