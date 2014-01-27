from __future__ import division
import matplotlib
matplotlib.use('Agg')
from sympy import *
import numpy as np
import difflib
import matplotlib.pyplot as plt


rAD,rAT,rE,sDE,sD = symbols(r'\rho_{AD} \rho_{AT} \rho_{E} \sigma_{DE} \sigma_{D}')
kde,kAD,kD,kdD,kE,L = symbols(r'k_{de} k_{AD} k_{D} k_{dD} k_{E} L')
gD,gE,Dd,De,k2 = symbols(r'\gamma_D \gamma_E Dd De k2')
u = 2

eq1 =  Eq(rAD,2*kde*sDE/(L*kAD))
eq2 = Eq(rAT, L*kAD*rAD/(2*(kD + kdD*(sD + sDE))))
eq3 = Eq(rE,kde*sDE/(kE*sD))
eq4 = Eq(gD,rAD*L + rAT*L + sD + sDE)
eq4 = rAD*L + rAT*L + sD + sDE - gD
eq5 = Eq(gE,rE*L + sDE)

eq4_modA = eq4.subs(rAT, L*kAD*rAD/(2*(kD + kdD*(sD + sDE)))).subs(rAD,2*kde*sDE/(L*kAD))
eq5_mod = eq5.subs(rE,kde*sDE/(kE*sD))
eq5_solved = Eq(sDE, gE/(L*kde/(sD*kE)+1))

f = open('learn-sympy.tex','w')
f.write(r"""\documentclass[letterpaper,onecolumn,amsmath,amssymb,pre]{revtex4-1}

\usepackage{breqn}
\usepackage{graphicx}

\begin{document}
""")

f.write(r'''These are the equations to solve for the homogeneoues situation (solve for the density values, that is):

\begin{equation}
%s
\end{equation}
\begin{equation}
%s
\end{equation}
\begin{equation}
%s
\end{equation}
\begin{equation}
%s
\end{equation}
\begin{equation}
%s
\end{equation}
''' % (latex(eq1),latex(eq2),latex(eq3),latex(eq4),latex(eq5)))

f.write(r'''This is Equation 4, but with the first two ($\rho_{AD}$ and $\rho_{AT}$) subbed in:
\begin{equation}
%s
\end{equation}

This is Equation 5, but with the third ($\rho_E$) subbed in:
\begin{equation}
%s
\end{equation}
or
\begin{equation}
%s
\end{equation}
''' % (latex(eq4_modA),latex(eq5_mod),latex(eq5_solved)))

eq4_modB = eq4_modA.subs(sDE, gE/(L*kde/(sD*kE)+1))

f.write('and plugging this answer into Equation 4:')
f.write('\n\\begin{dmath}\n')
f.write('%s\n' % latex(eq4_modB))
f.write('\end{dmath}\n')

f.write('''The following can be skipped but I left in in the latex as a check.
There is a hard coded bottom of this fraction term in the python so we want to make sure its
the same thing:\n''')
f.write('\n\\begin{dmath}\n')
f.write('%s\n' % latex(simplify(collect(ratsimp(eq4_modB), sD))))
f.write('\end{dmath}\n')

bottom = kAD*(L**2*sD*kdD*kde**2 + L**2*kD*kde**2 + L*gE*sD*kE*kdD*kde +
              2*L*sD**2*kE*kdD*kde + 2*L*sD*kD*kE*kde + gE*sD**2*kE**2*kdD +
              sD**3*kE**2*kdD + sD**2*kD*kE**2)


f.write('''\nAt this point solving for $\\sigma_D$ in this
equation takes too long so Im going to substitute in proper values for the variables first.
These values are right in the paper in the Model and Methods section.  The gammas, which
are equal to the amount of proteins in each unit area column, come out to:
$\\gamma_D = 1000*L/(\\pi*r^2)$ and $\\gamma_D = 350*L/(\\pi*r^2)$ where $r = .5\mu m$.
Also it's Rat simplified:
''')


pi = 3.14159265

def put_in_numbers(f):
    return f.subs(kAD,1).subs(kD,.025).subs(kde,.7).subs(kE,.093).subs(kdD,.0015).subs(gD, 1000*L/(pi*0.5*0.5)).subs(gE, 350*L/(pi*0.5*0.5)).subs(Dd,2.5).subs(De,2.5)#.subs(L, 1)

#print eq4_modB*bottom
f.write('\n\\begin{dmath}\n')
f.write('0 = %s\n' % latex(simplify(collect(put_in_numbers(ratsimp(eq4_modB*bottom)), sD))))
f.write('\end{dmath}\n')

f.write('''\nSympy refuses to solve this equation, complaining about it being a quartic, so I solved it in mathematica instead.  The solution that is dependent upon L is too big to put into sympy, so what Ive done instead is solved at a series of L values.  These and their corresponding values Ive plotted.''')


Larr = np.arange(0.05,2.05,0.05)
first_elements = np.array([0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,])
Larr = np.append(first_elements,Larr)
#print Larr

sDarr = np.array([0.966079,1.6284,2.10444,2.46,2.73616,2.95971,3.14835,
                  3.31381,3.46384,3.60356,4.83007,6.01588,7.20555,8.39814,9.59224,
                  10.7871,11.9824,13.1779,14.3735,15.5692,16.7648,17.9605,19.1562,
                  20.3518,21.5475,22.7431,23.9387,25.1343,26.3299,27.5254,
                  28.721,29.9165,31.112,32.3075,33.503,34.6985,35.894,37.0894,
                  38.2849,39.4804,40.6758,41.8712,43.0667,44.2621,45.4575,
                  46.6529,47.8483,49.0437,50.2391])

def sDcalc(L,c):
    return c[0]*L**3 + c[1]*L**2 + c[2]*L**1 + c[3]

polynomial_constants = np.polyfit(Larr, sDarr, 3)
sDarr_calc = np.array([])
for Li in Larr:
    sDarr_calc = np.append(sDarr_calc,sDcalc(Li,polynomial_constants))

#print sDarr_calc

plt.plot(Larr,sDarr)
plt.plot(Larr,sDarr_calc)
plt.xlabel('L')
plt.ylabel('sigmaD')
plt.xlim([0,.3])
plt.ylim([0,15])
plt.savefig('learn-sympy-plot.pdf')

f.write(r'''
\begin{figure}
\includegraphics[width=5in]{learn-sympy-plot}
\label{fig:10-01}
\end{figure}
''')

f.write('''\nI tried to fit the first density $\sigma_D$ to the
solution curve so I could get an easier function of L to put into
the matrix but it didn't work out (at least not yet) so now I'll
write a function to ultimately get k^2 of L\n''' )

f.write('''\nAt this point the calculations for k^2 as a function
of L happen in the sympy-calculations file. These are output and picked
up by this file, which will plot them.\n''')

data = np.loadtxt('sympy-calculations-out.txt')
print "len = ",len(data)

Larr = np.array([])
karr = np.array([])
k2arr = np.array([])
wavelength = np.array([])
for i in range(len(data)):
    Larr = np.append(Larr,data[i][0])
    wavelength = np.append(wavelength,2*pi/data[i][1])
    karr = np.append(karr,data[i][1])
    k2arr = np.append(k2arr,data[i][2])

for i in range(len(Larr)):
    if Larr[i] == 0.3 or Larr[i] == 0.25:
        print "L = %.2f and half wavelength = %f" % (Larr[i],wavelength[i]/2)

plt.figure()
plt.plot(Larr,wavelength/2)
plt.xlabel('Slab thickness ($\mu m$)')
plt.ylabel('Half wavelength ($\mu m$)')
plt.xlim([.1,.7])
plt.ylim([0,5])
plt.savefig('sympy-plot-wavelength.pdf')

f.write(r'''
\begin{figure}
\includegraphics[width=5in]{sympy-plot-k}
\label{fig:10-01}
\end{figure}
''')



f.write(r'\end{document}')




# def cut_string_for_latex(ans):
#     size_ch = 100
#     num_chunks = floor(len(str(ans))/size_ch)
#     ans_list = []
#     mod = 0
#     old_mod = 0
#     for i in range(num_chunks):
#         mod = 0
#         while ans[(i+1)*size_ch+mod] != ' ' and (i+1)*size_ch+mod < len(ans)-1:
#             mod += 1
#         char_list = ans[i*size_ch+old_mod : (i+1)*size_ch+mod]
#         ans_list.append(''.join(char_list))
#         old_mod = mod
#     ans_list.append(''.join(ans[(num_chunks)*size_ch+old_mod:]))
#     return ans_list



