A set of helper functions and code that will help solve linear systems problems

_note_ the ipynb files are more complete than the .py files.

This repo will help you do the following and more:
1. State transition matrix (both continuous time and discrete)
2. Homogenous response (aka free response or zero input response) of a state space model both continuous and discrete
3. Forced response (aka zero state response) of a state space model
4. Output response of a state space model
5. Determine and Contralabiltiy and Observability of a system, and find the contralabilty  an
6. Get a neat matrix inverse
7. Misc matrix and other functions

The documented code is in the jupyter notebook file `NeatStateSpaceFunctions.ipynb`
which I extracted that code from the notebook and put it in a Python file called `state_space_functions.py` for people who want to import it as a module.

The file `MatrixFunctions.ipynb` which is the same code as `MatrixFunctions` is the matrix helper function such as get_inverse() and others.

The `tests.py` attempts to test `state_space_functions.py` by making it solve state space model homogenous, forced, and output responses and checking the solution. It is incomplete but has few good tests.

The file `StateSpaceModel.ipynb`  has more functions and examples but it is an undocumented, unorganized version unlike `NeatStateSpaceFunctions.ipynb` I didn't remove it because it has examples which you can see to know how to use this code. You can see the example if you scroll down to the **Example usage** section in the `StateSpaceModel.ipynb`file


*Overall, this code helped me understand and solve Linear Systems problems so I wanted to share it, by no means is this a complete package or library.*
