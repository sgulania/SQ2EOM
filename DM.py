from sympy import *                                                             
from sympy.physics.secondquant import F, Fd, wicks, Commutator,evaluate_deltas,AntiSymmetricTensor                  
from sympy.physics.quantum import Commutator, Dagger, Operator
from fractions import Fraction
from sympy.physics.secondquant import simplify_index_permutations
from sympy.physics.secondquant import PermutationOperator
from sympy.physics.secondquant import substitute_dummies
from sympy.printing.latex import LatexPrinter, print_latex
from IPython.display import display, Markdown
import BCH

def OPDM(L,R,flavor):                                                           

    display(Markdown
        (rf""" Computing OPDM for {flavor} (skipping summation for dummy variables)"""))

    i,j = symbols('i,j' , below_fermi=True)
    a,b = symbols('a,b' , above_fermi=True)
 
    PermutList = [PermutationOperator(i,j),PermutationOperator(a,b)]

    oo = Fd(i)*F(j)
    cc = BCH.level(oo,"SD")
    g_oo = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    g_oo = simplify_index_permutations(g_oo,PermutList)
    index_rule = {'below':  'klmno','above':  'abcde'}
    g_oo = substitute_dummies(g_oo,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ij}')
    final_eq = Eq(gamma, g_oo)
    display(final_eq)

    ov = Fd(i)*F(a)
    cc = BCH.level(ov,"SD")
    g_ov = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    g_ov = simplify_index_permutations(g_ov,PermutList)
    index_rule = {'below':  'jklmn','above':  'bcdef'}
    g_ov = substitute_dummies(g_ov,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ia}')
    final_eq = Eq(gamma, g_ov)
    display(final_eq)
    
    vo = Fd(a)*F(i)
    cc = BCH.level(vo,"SD")
    g_vo = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    g_vo = simplify_index_permutations(g_vo,PermutList)
    index_rule = {'below':  'jklmn','above':  'bcdef'}
    g_vo = substitute_dummies(g_vo,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ai}')
    final_eq = Eq(gamma, g_vo)
    display(final_eq)
    
    vv = Fd(a)*F(b)
    cc = BCH.level(vv,"SD")
    g_vv = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    g_vv = simplify_index_permutations(g_vv,PermutList)
    index_rule = {'below':  'ijklm','above':  'cdefg'}
    g_vv = substitute_dummies(g_vv,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ab}')
    final_eq = Eq(gamma, g_vv)
    display(final_eq)


def TPDM(L,R,flavor):                                                           

    display(Markdown
        (rf""" Computing TPDM for {flavor} (skipping summation for dummy variables)"""))
    
    i,j,k,l = symbols('i,j,k,l' , below_fermi=True)
    a,b,c,d = symbols('a,b,c,d' , above_fermi=True)

    oooo = Fd(i)*Fd(j)*F(l)*F(k)
    cc = BCH.level(oooo,"SD")
    g_oooo = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    PermutList = [PermutationOperator(i,j),PermutationOperator(i,k), \
               PermutationOperator(i,l),PermutationOperator(j,k), \
               PermutationOperator(j,l),PermutationOperator(k,l)]
    g_oooo = simplify_index_permutations(g_oooo,PermutList)
    index_rule = {'below':  'mnop','above':  'abcde'}
    g_oooo = substitute_dummies(g_oooo,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ijkl}')
    final_eq = Eq(gamma, g_oooo)
    display(final_eq)
    
    
    ooov = Fd(i)*Fd(j)*F(a)*F(k)
    cc = BCH.level(ooov,"SD")
    g_ooov = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    PermutList = [PermutationOperator(i,j),PermutationOperator(i,k), \
               PermutationOperator(j,k)]
    g_ooov = simplify_index_permutations(g_ooov,PermutList)
    index_rule = {'below':  'lmnop','above':  'bcdef'}
    g_ooov = substitute_dummies(g_ooov,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ijka}')
    final_eq = Eq(gamma, g_oo)
    display(final_eq)
    
    ooov = Fd(i)*Fd(a)*F(k)*F(j)
    cc = BCH.level(ooov,"SD")
    g_ovoo = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    PermutList = [PermutationOperator(i,j),PermutationOperator(i,k), \
               PermutationOperator(j,k)]
    g_ovoo = simplify_index_permutations(g_ovoo,PermutList)
    index_rule = {'below':  'lmnop','above':  'bcdef'}
    g_ovoo = substitute_dummies(g_ovoo,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{iajk}')
    final_eq = Eq(gamma, g_ovoo)
    display(final_eq)
    
    ovov = Fd(i)*Fd(a)*F(b)*F(j)
    cc = BCH.level(ovov,"SD")
    g_ovov = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    PermutList = [PermutationOperator(i,j),PermutationOperator(a,b)]
    g_ovov = simplify_index_permutations(g_ovov,PermutList)
    index_rule = {'below':  'klmno','above':  'cdef'}
    g_ovov = substitute_dummies(g_ovov,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{iajb}')
    final_eq = Eq(gamma, g_ovov)
    display(final_eq)
    
    ovvv = Fd(i)*Fd(a)*F(c)*F(b)
    cc = BCH.level(ovvv,"SD")
    g_ovvv = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    PermutList = [PermutationOperator(a,b),PermutationOperator(a,c), \
               PermutationOperator(b,c)]
    g_ovvv = simplify_index_permutations(g_ovvv,PermutList)
    index_rule = {'below':  'jklmn','above':  'defg'}
    g_ovvv = substitute_dummies(g_ovvv,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{iabc}')
    final_eq = Eq(gamma, g_ovvv)
    display(final_eq)
    
    oovv = Fd(i)*Fd(j)*F(b)*F(a)
    cc = BCH.level(oovv,"SD")
    g_oovv = evaluate_deltas(wicks(L*cc*R , keep_only_fully_contracted=True))
    PermutList = [PermutationOperator(i,j),PermutationOperator(a,b)]
    g_oovv = simplify_index_permutations(g_oovv,PermutList)
    index_rule = {'below':  'klmn','above':  'cdefg'}
    g_oovv = substitute_dummies(g_oovv,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_{ijab}')
    final_eq = Eq(gamma, g_oovv)
    display(final_eq)
