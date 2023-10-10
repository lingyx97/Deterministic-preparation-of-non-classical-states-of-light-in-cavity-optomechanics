from __future__ import division
from pyparsing import (Literal, CaselessLiteral, Word, Combine, Group, Optional,
                       ZeroOrMore, Forward, nums, alphas, oneOf, delimitedList)
import math
import operator
from numpy import *
from qutip import *

from generateOps import geOps

__author__ = 'Paul McGuire'
__version__ = '$Revision: 0.0 $'
__date__ = '$Date: 2009-03-20 $'
__source__ = '''http://pyparsing.wikispaces.com/file/view/fourFn.py
http://pyparsing.wikispaces.com/message/view/home/15549426
'''
__note__ = '''
All I've done is rewrap Paul McGuire's fourFn.py as a class, so I can use it
more easily in other places.
'''
def nonCommutativeTimes(*argList):
    result=1
    for i in argList[::-1]:
        # print(i)
        result=i*result
    return result

class NumericStringParser(object):
    '''
    Most of this code comes from the fourFn.py pyparsing example

    '''

    def pushFirst(self, strg, loc, toks):
        self.exprStack.append(toks[0])

    def pushUMinus(self, strg, loc, toks):
        if toks and toks[0] == '-':
            self.exprStack.append('unary -')

    def __init__(self,dimHS,args):
        """
        expop   :: '^'
        multop  :: '*' | '/'
        addop   :: '+' | '-'
        integer :: ['+' | '-'] '0'..'9'+
        atom    :: a | ad | b | bd | PI | id |I | E | real | fn '(' expr ')' | '(' expr ')'
        factor  :: atom [ expop factor ]*
        term    :: factor [ multop factor ]*
        expr    :: term [ addop term ]*
        """
        point = Literal(".")
        e = Literal("E")
        fnumber = Combine(Word("+-" + nums, nums) +
                          Optional(point + Optional(Word(nums))) +
                          Optional(e + Word("+-" + nums, nums)))
        ident = Word(alphas, alphas + nums + "_$")
        plus = Literal("+")
        minus = Literal("-")
        mult = Literal("*")
        div = Literal("/")
        lpar = Literal("(").suppress()
        rpar = Literal(")").suppress()
        addop = plus | minus
        multop = mult | div
        expop = Literal("^")
        pi = Literal("Pi")
        
        # comma=Literal(",").suppress()
        
        a=Literal("a")
        ad=Literal("ad")
        b=Literal("b")
        bd=Literal("bd")
        idM=Literal("id")
        imagUnit = Literal("I")
        eta = Literal("eta")
        k = Literal("k")
        psi = Literal("psi")
        # t=Literal("t")
        
        
        expr = Forward()
        atom = ((Optional(oneOf("- +")) +
                 (ident + lpar + Group(delimitedList(expr, delim=r',')) + rpar | pi | e | fnumber | ad | a |
                    bd | b | idM | imagUnit | eta | k | psi).setParseAction(self.pushFirst))
                | Optional(oneOf("- +")) + Group(lpar + expr + rpar)
                ).setParseAction(self.pushUMinus)
        # by defining exponentiation as "atom [ ^ factor ]..." instead of
        # "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-right
        # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = Forward()
        factor << atom + \
            ZeroOrMore((expop + factor).setParseAction(self.pushFirst))
        term = factor + \
            ZeroOrMore((multop + factor).setParseAction(self.pushFirst))
        expr << term + \
            ZeroOrMore((addop + term).setParseAction(self.pushFirst))
        # addop_term = ( addop + term ).setParseAction( self.pushFirst )
        # general_term = term + ZeroOrMore( addop_term ) | OneOrMore( addop_term)
        # expr <<  general_term
        self.bnf = expr
        # map operator symbols to corresponding arithmetic operations
        epsilon = 1e-12
        self.opn = {"+": operator.add,
                    "-": operator.sub,
                    "*": operator.mul,
                    "/": operator.truediv,
                    "^": operator.pow}
        self.fn = {"sin": math.sin,
                   "cos": math.cos,
                   "tan": math.tan,
                   "exp": math.exp,
                   "abs": abs,
                #    "trunc": lambda a: int(a),
                #    "round": round,
                #    "sgn": lambda a: abs(a) > epsilon and cmp(a, 0) or 0
                   "Nct":nonCommutativeTimes
                   }
        self.args={
            "eta":args[0],
            "k":args[1],
            "psi":args[2],
            # "t":args[3]
        }
        
        adN,aN,bdN,bN,idmN=geOps(dimHS)
        
        self.quantumOP={
            "ad":adN,
            "a":aN,
            "bd":bdN,
            "b":bN,
            "id":idmN
        }

    def evaluateStack(self, s):
        op = s.pop()
        if op == 'unary -':
            return -self.evaluateStack(s)
        if op in "+-*/^":
            op2 = self.evaluateStack(s)
            op1 = self.evaluateStack(s)
            if op=="^" and int(real(op2))-op2==0:
                op2=int(real(op2))
            return self.opn[op](op1, op2)
        elif op == "Pi":
            return math.pi  # 3.1415926535
        elif op == "E":
            return math.e  # 2.718281828
        elif op=="I":
            return 1j
        elif op in self.quantumOP:
            return self.quantumOP[op]
        elif op in self.args:
            return self.args[op]
        elif op in self.fn:
            tempArg=[]
            while s[-1] in "adbd^" and len(s)>=1:
                # print(op)
                # input(s[-5:])
                tempArg.append(self.evaluateStack(s))
                if len(s)==0:
                    break
            return self.fn[op](*tempArg)
            # return self.fn[op](x)
        elif op[0].isalpha():
            return 0
        else:
            return float(op)

    def eval(self, num_string, parseAll=True):
        self.exprStack = []
        results = self.bnf.parseString(num_string, parseAll)
        val = self.evaluateStack(self.exprStack[:])
        return val
    
def Nct(*ops):
    res=1
    for op in ops:
        res=op*res
    return res
    
if __name__ == '__main__':
    dimHS=10
    ad,a,bd,b,idm=geOps(dimHS)
    # print(sum(abs((Nct(ad,a**2)-a**2*ad).full())))
    # exit()
    eta=4
    k=1/20
    psi=0.1
    data='((4*eta*k^4*Pi*Nct(b^2, a))/E^(I*psi) + 4*E^(I*psi)*eta*k^4*Pi*Nct(ad, b^2) + (4*eta*k^4*Pi*Nct(bd^2, a))/E^(I*psi) + 4*E^(I*psi)*eta*k^4*Pi*Nct(bd^2, ad))/2 + (4*eta*k^4*Pi*Nct(b, bd^3, a))/(3*E^(I*psi)) + (4*eta*k^4*Pi*Nct(b^3, bd, a))/(3*E^(I*psi)) - (4*E^(I*psi)*eta*k^4*Pi*Nct(bd, ad, b^3))/3 - (4*E^(I*psi)*eta*k^4*Pi*Nct(bd^3, ad, b))/3 '
    # data='Nct(bd, b^2, a)'
    nsp=NumericStringParser(dimHS,[eta,k,psi])
    test=nsp.eval(data)
    testCompare=((4*eta*k**4*pi*a*b**2)/exp(1j*psi) + 4*exp(1j*psi)*eta*k**4*pi*Nct(ad, b**2) + (4*eta*k**4*pi*Nct(bd**2, a))/exp(1j*psi) + 4*exp(1j*psi)*eta*k**4*pi*Nct(bd**2, ad))/2 + (4*eta*k**4*pi*Nct(b, bd**3, a))/(3*exp(1j*psi)) + (4*eta*k**4*pi*Nct(b**3, bd, a))/(3*exp(1j*psi)) - (4*exp(1j*psi)*eta*k**4*pi*Nct(bd, ad, b**3))/3 - (4*exp(1j*psi)*eta*k**4*pi*Nct(bd**3, ad, b))/3
    # testCompare=bd*b**2*a
    # testCompare=a*b**2*bd
    print(test)
    print(testCompare-test)
    print(sum(abs((testCompare-test).full())))
    print(sum(abs((testCompare+testCompare.dag()).full())))
    
    