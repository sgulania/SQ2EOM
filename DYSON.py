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

def OPTDM(Lf1,Rf1,Lf2,Rf2,flavor1,flavor2):                                                           

    display(Markdown
        (rf""" Computing Dyson OPTDM between {flavor1} $\rightarrow$ {flavor2} (skipping summation for dummy variables)"""))


    i = symbols('i' , below_fermi=True)
    a = symbols('a' , above_fermi=True)
 
    index_rule = {'below':  'jklmn','above':  'bcde'}

    oo = Fd(i)
    cc = BCH.level(oo,"SD")
    g_oo = evaluate_deltas(wicks(Lf2*cc*Rf1 , keep_only_fully_contracted=True))
    g_oo = substitute_dummies(g_oo,new_indices=True, pretty_indices=index_rule)
    gamma = Symbol('\gamma_i^{R}')
    final_eq = Eq(gamma, g_oo)
    display(final_eq)
                                                  
    ov = Fd(a)                                                                     
    cc = BCH.level(ov,"SD")                                                        
    g_ov = evaluate_deltas(wicks(Lf2*cc*Rf1 , keep_only_fully_contracted=True))      
    index_rule = {'below':  'jklmn','above':  'bcdef'}                             
    g_ov = substitute_dummies(g_ov,new_indices=True, pretty_indices=index_rule)    
    gamma = Symbol('\gamma_a^{R}')
    final_eq = Eq(gamma, g_ov)
    display(final_eq)
                                                      
    vo = F(i)                                                                      
    cc = BCH.level(vo,"SD")                                                        
    g_vo = evaluate_deltas(wicks(Lf1*cc*Rf2 , keep_only_fully_contracted=True))      
    index_rule = {'below':  'jklmn','above':  'bcdef'}                             
    g_vo = substitute_dummies(g_vo,new_indices=True, pretty_indices=index_rule)    
    gamma = Symbol('\gamma_i^{L}')
    final_eq = Eq(gamma, g_vo)
    display(final_eq)
                                                     
    vv = F(a)                                                                      
    cc = BCH.level(vv,"SD")                                                        
    g_vv = evaluate_deltas(wicks(Lf1*cc*Rf2 , keep_only_fully_contracted=True))      
    index_rule = {'below':  'ijklm','above':  'cdefg'}                             
    g_vv = substitute_dummies(g_vv,new_indices=True, pretty_indices=index_rule)    
    gamma = Symbol('\gamma_a^{L}')
    final_eq = Eq(gamma, g_vv)
    display(final_eq)
