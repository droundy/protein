from __future__ import division
import sys
import math

print "Just to make sure, the command line should read:"
print "python quadrupole-for-stadium.py <Volume> <radius>"
print "or"
print "python quadrupole-for-stadium.py <length_of_rectangle> <radius> L"

Vol = 0
L = 0
L_given = False

for arg in sys.argv:
    if arg == "L":
        L_given = True

pi = math.pi
T = 0.25 # thickness of pancake cell
R = float(sys.argv[2])

if L_given:
    L = float(sys.argv[1])
    Vol = L*(2.0*R*T + pi*T*T) + pi*R*R*T + pi*pi*T*T*R
else:
    Vol = float(sys.argv[1])
    L = ( (Vol - pi*R*R*T - pi*pi*T*T*R) / (2.0*R*T + pi*T*T) )

Qzz = L*L*L*R/3.0 + (L*L*pi/2.0)*(R*R + R*R*R*R/4.0)

Qyy = 4.0*L*R*R*R/3.0 + pi*R*R*R*R/2.0

print ""
print "Volume = ",Vol
print "Radius of the stadium shape R = ",R
print "Full length of rectangle part of stadium shape L = ",L
print "Ratio Qzz/Qyy = ",Qzz/Qyy
print ""

# From running the quad-for-diag script, I get answers:
# For 18_65-28_65 94
# Total cell volume is = 4.123500
# Qzz/Qyy = 3.681

# For 18_55-18_55 95
# total cell volume is 4.195
# Qzz/Qyy = 3.858

# Results of this script:
# With Vol = 4.123 and R = 1.18 and L = 2.93 the Qzz/Qyy = 3.72
# which is close enough for 18_65-28_65 94

# With Vol = 4.195 and R = 1.18 and L = 3.02 the Qzz/Qyy = 3.90
# which is close enough for 18_55-18_55 95

# These are so close though that I'll just do R = 1.18 and L = 3.00
# and call it good.



# But all of the above is for the quadrupole incorporating the 3 and 2 factors,
# So the below is more accurate (doesn't have those factors)
# For 18_55-18_55 95
# total cell volume is 4.195
# Qzz/Qyy = 2.240

# For 18_65-28_65 94
# Total cell volume is = 4.123500
# Qzz/Qyy = 3.611

# Results of this script:
# With Vol = 4.195 and R = 1.32 and L = 2.35 the Qzz/Qyy = 2.29
# which is close enough for 18_55-18_55 95

# With Vol = 4.118 and R = 1.18 and L = 2.92 the Qzz/Qyy = 3.70
# which is close enough for 18_65-28_65 94

# These are so close though that I'll just do R = 1.18 and L = 3.00
# and call it good.
