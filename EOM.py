from sympy import *                                                             
from sympy.physics.secondquant import F, Fd,AntiSymmetricTensor,NO
from fractions import Fraction                                                  

def R0(expr):

 if   expr == "IP":
  R0 = 0
  return R0
 elif expr == "DIP":
  R0 = 0
  return R0
 elif expr == "EA":
  R0 = 0
  return R0
 elif expr == "DEA":
  R0 = 0
  return R0
 elif expr == "EE":
  R0 = Fraction(1, 1)*AntiSymmetricTensor('R0',(),())
  return R0
 elif expr == "CCSD":
  R0 = 1
  return R0


def R1(expr):                                                           

 i1,i2,i3 = symbols('i1,i2,i3' ,below_fermi=True, cls=Dummy)                          
 a1,a2,a3 = symbols('a1,a2,i3' ,above_fermi=True, cls=Dummy)

 if   expr == "IP":                                                             
  R1 = Fraction(1, 1)*AntiSymmetricTensor('r',(),(i1,))*(F(i1))                                                            
  return R1                                                                     
 elif expr == "DIP":                                                            
  R1 = Fraction(1, 2)*AntiSymmetricTensor('r',(),(i1,i2))*F(i2)*F(i1)          
  return R1
 elif expr == "EA":
  R1 = Fraction(1, 1)*AntiSymmetricTensor('r',(a1,),())*Fd(a1)
  return R1
 elif expr == "DEA":
  R1 = Fraction(1, 2)*AntiSymmetricTensor('r',(a1,a2),())*NO(Fd(a1)*Fd(a2))
  return R1                                                                   
 elif expr == "EE":
  R1 = Fraction(1, 1)*AntiSymmetricTensor('r',(a1,),(i1,))*Fd(a1)*F(i1)
  return R1
 elif expr == "CCSD":
  R1 = 0
  return R1
                                                                                
def R2(expr):                                                           
 i1,i2,i3,i4,i5 = symbols('i1,i2,i3,i4,i5' ,below_fermi=True, cls=Dummy)
 a1,a2,a3,a4,a5 = symbols('a1,a2,a3,a4,a5' ,above_fermi=True, cls=Dummy)                          
 if   expr == "IP":                                                             
  R2 = Fraction(1, 2)*AntiSymmetricTensor('r',(a1,),(i2,i3))*Fd(a1)*F(i3)*F(i2)                  
  return R2                                                                     
 elif expr == "DIP":                                                            
  R2 = Fraction(1, 6)*AntiSymmetricTensor('r',(a1,),(i3,i4,i5))*Fd(a1)*F(i5)*F(i4)*F(i3)          
  return R2
 elif expr == "EA":
  R2 = Fraction(1, 2)*AntiSymmetricTensor('r',(a2,a3),(i1,))*Fd(a2)*Fd(a3)*F(i1)
  return R2
 elif expr == "DEA":
  R2 = Fraction(1, 6)*AntiSymmetricTensor('r',(a3,a4,a5),(i1,))*NO(Fd(a3)*Fd(a4)*Fd(a5)*F(i1))
  return R2                                                                   
 elif expr == "EE":
  R2 = Fraction(1, 4)*AntiSymmetricTensor('r',(a2,a3),(i2,i3))*Fd(a2)*Fd(a3)*F(i3)*F(i2)
  return R2
 elif expr == "CCSD":
  R2 = 0
  return R2


def L0(expr):

 if   expr == "IP":
  L0 = 0
  return L0
 elif expr == "DIP":
  L0 = 0
  return L0
 elif expr == "EA":
  L0 = 0
  return L0
 elif expr == "DEA":
  L0 = 0
  return L0
 elif expr == "EE":
  L0 = 0
  return L0
 elif expr == "CCSD":
  L0 = 1
  return L0


def L1(expr):                                                           
 j1,j2,j3 = symbols('j1,j2,j3' ,below_fermi=True, cls=Dummy)                          
 b1,b2,b3 = symbols('b1,b2,b3' ,above_fermi=True, cls=Dummy)

 if   expr == "IP":                                                             
  L1 = Fraction(1, 1)*AntiSymmetricTensor('l',(j1,),())*Fd(j1)
  return L1                                                                     
 elif expr == "DIP":                                                            
  L1 = Fraction(1, 2)*AntiSymmetricTensor('l',(j1,j2),())*Fd(j1)*Fd(j2)          
  return L1
 elif expr == "EA":
  L1 = Fraction(1, 1)*AntiSymmetricTensor('l',(),(b1,))*F(b1)
  return L1
 elif expr == "DEA":
  L1 = Fraction(1, 2)*AntiSymmetricTensor('l',(),(b1,b2))*F(b2)*F(b1)
  return L1                                                                   
 elif expr == "EE":
  L1 = Fraction(1, 1)*AntiSymmetricTensor('l',(j1,),(b1,))*Fd(j1)*F(b1)
  return L1
 elif expr == "CCSD":
  L1 = Fraction(1, 1)*AntiSymmetricTensor('l',(j1,),(b1,))*Fd(j1)*F(b1)
  return L1
                                                                                
def L2(expr):                                                           
 j1,j2,j3,j4,j5 = symbols('j1,j2,j3,j4,j5' ,below_fermi=True, cls=Dummy)
 b1,b2,b3,b4,b5 = symbols('b1,b2,b3,b4,b5' ,above_fermi=True, cls=Dummy)                          
 if   expr == "IP":                                                             
  L2 = Fraction(1, 2)*AntiSymmetricTensor('l',(j2,j3),(b1,))*Fd(j2)*Fd(j3)*F(b1)                  
  return L2                                                                     
 elif expr == "DIP":                                                            
  L2 = Fraction(1, 6)*AntiSymmetricTensor('l',(j3,j4,j5),(b1,))*Fd(j3)*Fd(j4)*Fd(j5)*F(b1)          
  return L2
 elif expr == "EA":
  L2 = Fraction(1, 2)*AntiSymmetricTensor('l',(j1,),(b2,b3))*Fd(j1)*F(b3)*F(b2)
  return L2
 elif expr == "DEA":
  L2 = Fraction(1, 6)*AntiSymmetricTensor('l',(j1,),(b3,b4,b5))*Fd(j1)*F(b5)*F(b4)*F(b3)
  return L2                                                                   
 elif expr == "EE":
  L2 = Fraction(1, 4)*AntiSymmetricTensor('l',(j2,j3),(b2,b3))*Fd(j2)*Fd(j3)*F(b3)*F(b2)
  return L2
 elif expr == "CCSD":
  L2 = Fraction(1, 4)*AntiSymmetricTensor('l',(j2,j3),(b2,b3))*Fd(j2)*Fd(j3)*F(b3)*F(b2)
  return L2

