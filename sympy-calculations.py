from __future__ import division
import matplotlib
matplotlib.use('Agg')
from sympy import *
import numpy as np

rAD,rAT,rE,sDE,sD = symbols(r'\rho_{AD} \rho_{AT} \rho_{E} \sigma_{DE} \sigma_{D}')
kde,kAD,kD,kdD,kE,L = symbols(r'k_{de} k_{AD} k_{D} k_{dD} k_{E} L')
gD,gE,Dd,De,k2 = symbols(r'\gamma_D \gamma_E Dd De k2')

eq1 =  Eq(rAD,2*kde*sDE/(L*kAD))
eq2 = Eq(rAT, L*kAD*rAD/(2*(kD + kdD*(sD + sDE))))
eq3 = Eq(rE,kde*sDE/(kE*sD))
eq4 = Eq(gD,rAD*L + rAT*L + sD + sDE)
eq4 = rAD*L + rAT*L + sD + sDE - gD
eq5 = Eq(gE,rE*L + sDE)

eq4_modA = eq4.subs(rAT, L*kAD*rAD/(2*(kD + kdD*(sD + sDE)))).subs(rAD,2*kde*sDE/(L*kAD))
eq5_mod = eq5.subs(rE,kde*sDE/(kE*sD))
eq5_solved = Eq(sDE, gE/(L*kde/(sD*kE)+1))

def put_in_numbers(f):
    return f.subs(kAD,1).subs(kD,.025).subs(kde,.7).subs(kE,.093).subs(kdD,.0015).subs(gD, 1000*L/(pi*0.5*0.5)).subs(gE, 350*L/(pi*0.5*0.5)).subs(Dd,2.5).subs(De,2.5)#.subs(L, 1)

Larr = np.arange(0.05,2.05,0.05)
first_elements = np.array([0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,])
Larr = np.append(first_elements,Larr)

sDarr = np.array([0.966079,1.6284,2.10444,2.46,2.73616,2.95971,3.14835,
                  3.31381,3.46384,3.60356,4.83007,6.01588,7.20555,8.39814,9.59224,
                  10.7871,11.9824,13.1779,14.3735,15.5692,16.7648,17.9605,19.1562,
                  20.3518,21.5475,22.7431,23.9387,25.1343,26.3299,27.5254,
                  28.721,29.9165,31.112,32.3075,33.503,34.6985,35.894,37.0894,
                  38.2849,39.4804,40.6758,41.8712,43.0667,44.2621,45.4575,
                  46.6529,47.8483,49.0437,50.2391])

# Larr = Larr[0:3]
# sDarr = sDarr[0:3]

k2arr = np.array([])
karr = np.array([])
for i in range(len(Larr)):
    sDE_ans = solve(put_in_numbers(eq5_mod).subs(L,Larr[i]).subs(sD,sDarr[i]),sDE)[0]
    eq3 = put_in_numbers(Eq(rE,kde*sDE/(kE*sD))).subs(sD,sDarr[i]).subs(sDE,sDE_ans)
    rE_ans = solve(eq3)[0]
    eq1 = put_in_numbers(Eq(rAD,2*kde*sDE/(L*kAD))).subs(sD,sDarr[i]).subs(sDE,sDE_ans).subs(L,Larr[i])
    rAD_ans = solve(eq1)[0]
    eq2 = put_in_numbers(Eq(rAT, L*kAD*rAD/(2*(kD + kdD*(sD + sDE))))).subs(sD,sDarr[i]).subs(sDE,sDE_ans).subs(L,Larr[i]).subs(rAD,rAD_ans)
    rAT_ans = solve(eq2)[0]

    M = Matrix([[-Dd*k2-kAD,0,0,0,2*kde/L],[kAD, -Dd*k2-2*kD/L-2*kdD*(sD+sDE)/L, 0, -2*kdD*rAT/L, 2*kdD*rAT/L],
                [0,0,-De*k2-2*kE*sD/L,-2*kE*rE/L,2*kde/L],[0,kD+kdD*(sD+sDE),-kE*sD,-kE*rE+kdD*rAT,kdD*rAT],
                [0,0,kE*sD,kE*rE,-kde]])

    M_det = collect(put_in_numbers(M.det()).subs(rE,rE_ans).subs(rAD,rAD_ans).subs(rAT,rAT_ans).subs(sD,sDarr[i]).subs(sDE,sDE_ans).subs(L,Larr[i]),k2)
    k2_from_det = solve(M_det)
    pos = False
    for j in k2_from_det:
        if j>0:
            k2_from_det = j
            if pos == True:
                print "There were two positive solutions to the det equation."
                print exit(1)
                pos = True
    k2arr = np.append(k2arr,k2_from_det)
    k = sqrt(k2_from_det)
    karr = np.append(karr,k)
    print ""
    print Larr[i]
    print k2_from_det
    print ""

f = open('sympy-calculations-out.txt', 'w')
for i in range(len(Larr)):
    f.write('%f\t%f\t%f\n' % (Larr[i],karr[i],k2arr[i]))
f.close()
