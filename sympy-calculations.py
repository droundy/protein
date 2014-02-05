from __future__ import division
import matplotlib
matplotlib.use('Agg')
from sympy import *
import numpy as np
import os

rAD,rAT,rE,sDE,sD = symbols(r'\rho_{AD} \rho_{AT} \rho_{E} \sigma_{DE} \sigma_{D}')
kde,kAD,kD,kdD,kE,L = symbols(r'k_{de} k_{AD} k_{D} k_{dD} k_{E} L')
gD,gE,Dd,De,k2 = symbols(r'\gamma_D \gamma_E Dd De k2')
x = symbols(r'x')

eq1 =  Eq(rAD,2*kde*sDE/(L*kAD))

eq2 = Eq(rAT, L*kAD*rAD/(2*(kD + kdD*(sD + sDE))))
eq3 = Eq(rE,kde*sDE/(kE*sD))
eq4 = Eq(gD,rAD*L + rAT*L + 2*sD + 2*sDE)
eq4 = rAD*L + rAT*L + 2*sD + 2*sDE - gD
eq5 = Eq(gE,rE*L + 2*sDE)

eq4_modA = eq4.subs(rAT, solve(eq2, rAT)[0]).subs(rAD,solve(eq1, rAD)[0])
eq5_mod = eq5.subs(rE,solve(eq3, rE)[0])
eq5_solved = Eq(sDE, solve(eq5_mod, sDE)[0])

def put_in_numbers(f):
    return f.subs(kAD,1).subs(kD,.025).subs(kde,.7).subs(kE,.093).subs(kdD,.0015).subs(gD, 1000*L/(pi*0.5*0.5)).subs(gE, 350*L/(pi*0.5*0.5)).subs(Dd,2.5).subs(De,2.5)#.subs(L, 1)

Larr = np.arange(0.05,2.05,0.05)
first_elements = np.array([0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,])
Larr = np.append(first_elements,Larr)

eq4_modA_with_eq5_solved_input = eq4_modA.subs(sDE, solve(eq5_solved, sDE)[0])

mysD = 1.0
sDarr = np.zeros_like(Larr)
for i in xrange(len(Larr)):
    mysD = nsolve(put_in_numbers(eq4_modA_with_eq5_solved_input).subs(L, Larr[i]).subs(sD, x), 2*mysD)
    sDarr[i] = mysD
    print Larr[i], mysD

k2arr = np.zeros_like(Larr)
karr = np.zeros_like(Larr)
lambda_over_2 = np.zeros_like(Larr)

for i in range(len(Larr)):
    sDE_ans = solve(put_in_numbers(eq5_mod).subs(L,Larr[i]).subs(sD,sDarr[i]),sDE)[0]
    rE_ans = solve(put_in_numbers(eq3).subs(sD,sDarr[i]).subs(sDE,sDE_ans),rE)[0]
    rAD_ans = solve(put_in_numbers(eq1).subs(sD,sDarr[i]).subs(sDE,sDE_ans).subs(L,Larr[i]),rAD)[0]
    rAT_ans = solve(put_in_numbers(eq2).subs(sD,sDarr[i]).subs(sDE,sDE_ans).subs(L,Larr[i]).subs(rAD,rAD_ans),rAT)[0]
    print sDE_ans,rE_ans,rAD_ans,rAT_ans

    M = Matrix([[-Dd*k2-kAD, 0                             , 0               , 0             , 2*kde/L    ],
                [ kAD      , -Dd*k2-2*kD/L-2*kdD*(sD+sDE)/L, 0               , -2*kdD*rAT/L  ,-2*kdD*rAT/L],
                [ 0        , 0                             , -De*k2-2*kE*sD/L, -2*kE*rE/L    , 2*kde/L    ],
                [ 0        , kD+kdD*(sD+sDE)               , -kE*sD          , -kE*rE+kdD*rAT, kdD*rAT    ],
                [ 0        , 0                             , kE*sD           , kE*rE         , -kde       ]])
    Mnumerical = put_in_numbers(M).subs(rE,rE_ans).subs(rAD,rAD_ans).subs(rAT,rAT_ans).subs(sD,sDarr[i]).subs(sDE,sDE_ans).subs(L,Larr[i])
    def find_greatest_eigenvalue_real_part(k2val):
        xxx = np.array(Mnumerical.subs(k2, k2val).tolist(), dtype=np.complex)
        return np.real(np.linalg.eigvals(xxx)).max()
    def find_k2_of_instability_between(k2lo, k2hi):
        k2try = 0.5*(k2lo + k2hi)
        if k2hi - k2lo < 1e-7:
            return k2try
        grt = find_greatest_eigenvalue_real_part(k2try)
        #print '   ', k2try, grt
        if grt > 0:
            return find_k2_of_instability_between(k2try, k2hi)
        else:
            return find_k2_of_instability_between(k2lo, k2try)
    tryk2 = 0.1**2
    print Mnumerical.subs(k2, tryk2).eigenvals(), Mnumerical.subs(k2, tryk2).eigenvects(), Mnumerical.subs(k2, tryk2).det()
    print Mnumerical.subs(k2, tryk2)

    print 'numpy eigvals', find_greatest_eigenvalue_real_part(tryk2)
    print 'k2instability', find_k2_of_instability_between(0, 100**2)
    print 'k instability', np.sqrt(find_k2_of_instability_between(0, 100**2))
    print '*********** lambda/2 instability', Larr[i], np.pi/np.sqrt(find_k2_of_instability_between(0, 100**2))

    print 'numpy eigvals', find_greatest_eigenvalue_real_part(tryk2)
    k2arr[i] =  find_k2_of_instability_between(0, 100**2)
    karr[i] =  np.sqrt(find_k2_of_instability_between(0, 100**2))
    lambda_over_2[i] =  np.pi/np.sqrt(find_k2_of_instability_between(0, 100**2))


f = open('sympy-calculations-out.txt', 'w')
for i in range(len(Larr)):
    f.write('%f\t%f\t%f\n' % (Larr[i],karr[i],k2arr[i]))
f.close()

# f = open('sympy-calc-for-mathematica.tex','w')
# f.write(r"""\documentclass[letterpaper,onecolumn,amsmath,amssymb,pre]{revtex4-1}

# \usepackage{breqn}
# \usepackage{graphicx}

# \begin{document}
# """)

# f.write(r'''These are the equations to solve for the homogeneoues situation (solve for the density values, that is):

# \begin{equation}
# 0 = %s
# \end{equation}
# \begin{equation}
# 0 = %s
# \end{equation}
# \end{document}
# ''' % (latex(eq4_modA_with_eq5_solved_input), latex(put_in_numbers(eq4_modA_with_eq5_solved_input).subs(L, 0.5))))
# f.close()
# os.system('pdflatex sympy-calc-for-mathematica.tex > /dev/null')
