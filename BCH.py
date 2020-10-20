from sympy import *                                                             
from sympy.physics.secondquant import F, Fd,Commutator,AntiSymmetricTensor      
from sympy.physics.quantum import Commutator, Dagger, Operator                  
from fractions import Fraction                                                  
                                                                                
def level(main,expr):                                                           
 #setup T                                                                        
 k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12 = symbols('k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12' , below_fermi=True, cls=Dummy)
 d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 = symbols('d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12' , above_fermi=True, cls=Dummy)
                                                                                
 TA1 = AntiSymmetricTensor('T',(d1,),(k1,))*Fd(d1)*F(k1)                       
 TA2 = Fraction(1, 4)*AntiSymmetricTensor('T',(d2,d3),(k2,k3))*Fd(d2)*Fd(d3)*F(k3)*F(k2)
 TA  = TA1+TA2                                                                  
                                                                                
 TB1 = AntiSymmetricTensor('T',(d4,),(k4,))*Fd(d4)*F(k4)                       
 TB2 = Fraction(1, 4)*AntiSymmetricTensor('T',(d5,d6),(k5,k6))*Fd(d5)*Fd(d6)*F(k6)*F(k5)
 TB  = TB1+TB2                                                                  
                                                                                
 TC1 = AntiSymmetricTensor('T',(d7,),(k7,))*Fd(d7)*F(k7)                       
 TC2 = Fraction(1, 4)*AntiSymmetricTensor('T',(d8,d9),(k8,k9))*Fd(d8)*Fd(d9)*F(k9)*F(k8)
 TC  = TC1+TC2                                                                  
                                                                                
 TD1 = AntiSymmetricTensor('T',(d10,),(k10,))*Fd(d10)*F(k10)                   
 TD2 = Fraction(1, 4)*AntiSymmetricTensor('T',(d11,d12),(k11,k12))*Fd(d11)*Fd(d12)*F(k11)*F(k12)
 TD = TD1+TD2                                                                   
                                                                                
 # BCH expansion                                                                
 if expr == "SD":                                                               
  BCH_expansion = main + Commutator(main,TA) + Fraction(1, 2)*Commutator(Commutator(main,TA),TB) + \
                                Fraction(1, 6)*Commutator(Commutator(Commutator(main,TA),TB),TC) + \
                                Fraction(1, 24)*Commutator(Commutator(Commutator(Commutator(main,TA),TB),TC),TD)
  return BCH_expansion
