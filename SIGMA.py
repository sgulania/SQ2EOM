from sympy import *                                                             
from sympy.physics.secondquant import F, Fd, wicks, Commutator,evaluate_deltas,AntiSymmetricTensor,NO                  
from sympy.physics.quantum import Commutator, Dagger, Operator
from fractions import Fraction
from sympy.physics.secondquant import simplify_index_permutations
from sympy.physics.secondquant import PermutationOperator
from sympy.physics.secondquant import substitute_dummies
from IPython.display import display, Markdown
import BCH

def RVECTORS(R0,R1,R2,flavor):                                       
 display(Markdown
        (rf""" Computing right sigma amplitudes for {flavor} (skipping summation for dummy variables)"""))

 p, q, r, s = symbols('p,q,r,s', cls=Dummy)
 f = AntiSymmetricTensor('f', (p,), (q,))
 pr = NO((Fd(p)*F(q)))
 v = AntiSymmetricTensor('v', (p, q), (r, s))
 pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
 ham = f*pr + Rational(1, 4)*v*pqsr

 cc  = BCH.level(ham,"SD")
 
 i,j,k = symbols('i,j,k', below_fermi=True)
 a,b,c = symbols('a,b,c', above_fermi=True)

 if flavor == "IP": 
  sig11 = evaluate_deltas(wicks(Fd(i)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'abcdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
 if flavor == "EA": 
  sig11 = evaluate_deltas(wicks(F(a)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'bcdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j)]
  sig11 = evaluate_deltas(wicks(Fd(i)*Fd(j)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'klmno','above':  'abcdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
  sig11 = simplify_index_permutations(sig11,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b)]
  sig11 = evaluate_deltas(wicks(F(b)*F(a)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'cdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
  sig11 = simplify_index_permutations(sig11,PermutList)

 sigma_11 = Symbol('((\overline{H}_{SS}-E_{cc})R_{1})')
 final_eq = Eq(sigma_11, sig11)
 display(final_eq)


 if flavor == "IP": 
  PermutList = [PermutationOperator(i,j)]
  sig12 = evaluate_deltas(wicks(Fd(i)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'abcdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)
 if flavor == "EA": 
  PermutList = [PermutationOperator(a,b)]
  sig12 = evaluate_deltas(wicks(F(a)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'bcdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j),PermutationOperator(j,k),PermutationOperator(i,k)]
  sig12 = evaluate_deltas(wicks(Fd(i)*Fd(j)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'klmno','above':  'abcdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b),PermutationOperator(b,c),PermutationOperator(a,c)]
  sig12 = evaluate_deltas(wicks(F(b)*F(a)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'cdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)

 sigma_12 = Symbol('(\overline{H}_{SD}R_{2})')
 final_eq = Eq(sigma_12, sig12)
 display(final_eq)

 if flavor == "IP": 
  sig21 = evaluate_deltas(wicks(Fd(i)*Fd(j)*F(a)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'klmno','above':  'bcdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
 if flavor == "EA": 
  sig21 = evaluate_deltas(wicks(Fd(i)*F(b)*F(a)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'cdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j)]
  sig21 = evaluate_deltas(wicks(Fd(i)*Fd(j)*Fd(k)*F(a)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'lmno','above':  'bcdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
  sig21 = simplify_index_permutations(sig21,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b)]
  sig21 = evaluate_deltas(wicks(Fd(i)*F(c)*F(b)*F(a)*(cc*R1-R1*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'cdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
  sig21 = simplify_index_permutations(sig21,PermutList)

 sigma_21 = Symbol('(\overline{H}_{DS}R_{1})')
 final_eq = Eq(sigma_21, sig21)
 display(final_eq)


 if flavor == "IP": 
  PermutList = [PermutationOperator(i,j)]
  sig22 = evaluate_deltas(wicks(Fd(i)*Fd(j)*F(a)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'klmno','above':  'bcdefgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)
 if flavor == "EA": 
  PermutList = [PermutationOperator(a,b)]
  sig22 = evaluate_deltas(wicks(Fd(i)*F(b)*F(a)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'cdefgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j),PermutationOperator(j,k),PermutationOperator(i,k)]
  sig22 = evaluate_deltas(wicks(Fd(i)*Fd(j)*Fd(k)*F(a)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'lmno','above':  'bcdefgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b),PermutationOperator(b,c),PermutationOperator(a,c)]
  sig22 = evaluate_deltas(wicks(Fd(i)*F(c)*F(b)*F(a)*(cc*R2-R2*cc) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'defgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)

 sigma_22 = Symbol('((\overline{H}_{DD}-E_{cc})R_{2})')
 final_eq = Eq(sigma_22, sig22)
 display(final_eq)

def LVECTORS(L0,L1,L2,flavor):                                       
 display(Markdown
        (rf""" Computing left sigma amplitudes for {flavor} (skipping summation for dummy variables)"""))

 p, q, r, s = symbols('p,q,r,s', cls=Dummy)
 f = AntiSymmetricTensor('f', (p,), (q,))
 pr = NO((Fd(p)*F(q)))
 v = AntiSymmetricTensor('v', (p, q), (r, s))
 pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
 ham = f*pr + Rational(1, 4)*v*pqsr

 cc  = BCH.level(ham,"SD")
 
 i,j,k = symbols('i,j,k', below_fermi=True)
 a,b,c = symbols('a,b,c', above_fermi=True)

 if flavor == "IP": 
  sig11 = evaluate_deltas(wicks((cc*L1-L1*cc)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'abcdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
 if flavor == "EA": 
  sig11 = evaluate_deltas(wicks((cc*L1-L1*cc)*Fd(a) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'bcdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j)]
  sig11 = wicks(cc*L1*F(j)*F(i)-L1*cc*F(j)*F(i), simplify_kronecker_deltas=True,keep_only_fully_contracted=True)
  index_rule = {'below':  'klmno','above':  'abcdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
  sig11 = simplify_index_permutations(sig11,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b)]
  sig11 = evaluate_deltas(wicks((cc*L1-L1*cc)*Fd(a)*Fd(b) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'cdefg'}
  sig11 = substitute_dummies(sig11,new_indices=True, pretty_indices=index_rule)
  sig11 = simplify_index_permutations(sig11,PermutList)

 sigma_11 = Symbol('(L_{1}(\overline{H}_{SS}-E_{cc}))')
 final_eq = Eq(sigma_11, sig11)
 display(final_eq)


 if flavor == "IP": 
  PermutList = [PermutationOperator(i,j)]
  sig12 = evaluate_deltas(wicks((cc*L2-L2*cc)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'abcdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)
 if flavor == "EA": 
  PermutList = [PermutationOperator(a,b)]
  sig12 = evaluate_deltas(wicks((cc*L2-L2*cc)*Fd(a) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'bcdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j),PermutationOperator(j,k),PermutationOperator(i,k)]
  sig12 = wicks(cc*L2*F(j)*F(i)-L2*cc*F(j)*F(i) ,simplify_kronecker_deltas=True,keep_only_fully_contracted=True)
  index_rule = {'below':  'klmno','above':  'abcdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b),PermutationOperator(b,c),PermutationOperator(a,c)]
  sig12 = evaluate_deltas(wicks((cc*L2-L2*cc)*Fd(a)*Fd(b) , keep_only_fully_contracted=True))
  index_rule = {'below':  'ijklmno','above':  'cdefg'}
  sig12 = substitute_dummies(sig12,new_indices=True, pretty_indices=index_rule)
  sig12 = simplify_index_permutations(sig12,PermutList)

 sigma_12 = Symbol('(L_{2}\overline{H}_{DS})')
 final_eq = Eq(sigma_12, sig12)
 display(final_eq)

 if flavor == "IP": 
  sig21 = evaluate_deltas(wicks((cc*L1-L1*cc)*Fd(a)*F(j)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'klmno','above':  'bcdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
 if flavor == "EA": 
  sig21 = evaluate_deltas(wicks((cc*L1-L1*cc)*Fd(a)*Fd(b)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'cdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j)]
  sig21 = wicks(cc*L1*Fd(a)*F(k)*F(j)*F(i)-L1*cc*Fd(a)*F(k)*F(j)*F(i) ,simplify_kronecker_deltas=True,keep_only_fully_contracted=True)
  index_rule = {'below':  'lmno','above':  'bcdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
  sig21 = simplify_index_permutations(sig21,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b)]
  sig21 = evaluate_deltas(wicks((cc*L1-L1*cc)*Fd(a)*Fd(b)*Fd(c)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'cdefgh'}
  sig21 = substitute_dummies(sig21,new_indices=True, pretty_indices=index_rule)
  sig21 = simplify_index_permutations(sig21,PermutList)

 sigma_21 = Symbol('(L_{1}\overline{H}_{SD})')
 final_eq = Eq(sigma_21, sig21)
 display(final_eq)


 if flavor == "IP": 
  PermutList = [PermutationOperator(i,j)]
  sig22 = evaluate_deltas(wicks((cc*L2-L2*cc)*Fd(a)*F(j)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'klmno','above':  'bcdefgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)
 if flavor == "EA": 
  PermutList = [PermutationOperator(a,b)]
  sig22 = evaluate_deltas(wicks((cc*L2-L2*cc)*Fd(a)*Fd(b)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'cdefgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)
 if flavor == "DIP": 
  PermutList = [PermutationOperator(i,j),PermutationOperator(j,k),PermutationOperator(i,k)]
  sig22 = wicks(cc*L2*Fd(a)*F(k)*F(j)*F(i)-L2*cc*Fd(a)*F(k)*F(j)*F(i) , simplify_kronecker_deltas=True,keep_only_fully_contracted=True)
  index_rule = {'below':  'lmno','above':  'bcdefgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)
 if flavor == "DEA": 
  PermutList = [PermutationOperator(a,b),PermutationOperator(b,c),PermutationOperator(a,c)]
  sig22 = evaluate_deltas(wicks((cc*L2-L2*cc)*Fd(a)*Fd(b)*Fd(c)*F(i) , keep_only_fully_contracted=True))
  index_rule = {'below':  'jklmno','above':  'defgh'}
  sig22 = substitute_dummies(sig22,new_indices=True, pretty_indices=index_rule)
  sig22 = simplify_index_permutations(sig22,PermutList)

 sigma_22 = Symbol('(L_{2}(\overline{H}_{DD}-E_{cc}))')
 final_eq = Eq(sigma_22, sig22)
 display(final_eq)
