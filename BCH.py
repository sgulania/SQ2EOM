from sympy import *                                                             
from sympy.physics.secondquant import F, Fd,wicks,Commutator,AntiSymmetricTensor,NO,evaluate_deltas                  
from sympy.physics.secondquant import substitute_dummies
from sympy.physics.quantum import Commutator, Dagger, Operator
from fractions import Fraction
from IPython.display import display, Markdown


def get_CC_operators():
    """
    Returns a tuple (T1,T2) of unique operators.
    """
    i = symbols('i', below_fermi=True, cls=Dummy)
    a = symbols('a', above_fermi=True, cls=Dummy)
    t_ai = AntiSymmetricTensor('t', (a,), (i,))
    ai = NO(Fd(a)*F(i))
    i, j = symbols('i,j', below_fermi=True, cls=Dummy)
    a, b = symbols('a,b', above_fermi=True, cls=Dummy)
    t_abij = AntiSymmetricTensor('t', (a, b), (i, j))
    abji = NO(Fd(a)*Fd(b)*F(j)*F(i))

    T1 = t_ai*ai
    T2 = Rational(1, 4)*t_abij*abji
    return (T1, T2)


def level(H,expr):

    pretty_dummies_dict = {
    'above': 'defg',
    'below': 'lmno',
    'general': 'pqrst'
    }

    #display(Markdown
    #    (rf"""Calculating 4 nested commutators"""))
    C = Commutator

    T1, T2 = get_CC_operators()
    T = T1 + T2
    comm1 = wicks(C(H, T))
    comm1 = evaluate_deltas(comm1)
    comm1 = substitute_dummies(comm1)

    T1, T2 = get_CC_operators()
    T = T1 + T2
    comm2 = wicks(C(comm1, T))
    comm2 = evaluate_deltas(comm2)
    comm2 = substitute_dummies(comm2)

    T1, T2 = get_CC_operators()
    T = T1 + T2
    comm3 = wicks(C(comm2, T))
    comm3 = evaluate_deltas(comm3)
    comm3 = substitute_dummies(comm3)

    T1, T2 = get_CC_operators()
    T = T1 + T2
    comm4 = wicks(C(comm3, T))
    comm4 = evaluate_deltas(comm4)
    comm4 = substitute_dummies(comm4)

    eq = H + comm1 + comm2/2 + comm3/6 + comm4/24
    eq = eq.expand()
    eq = evaluate_deltas(eq)
    eq = substitute_dummies(eq, new_indices=True,
            pretty_indices=pretty_dummies_dict)

    return eq

