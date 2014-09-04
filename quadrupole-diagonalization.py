from __future__ import division
import sys
import math

print "Just to make sure, the command line should read:"
print "python quadrupole-diagonalization.py <Qzz> <Qyy> <Qyz>"

Qzz = float(sys.argv[1])
Qyy = float(sys.argv[2])
Qzy = float(sys.argv[3])
Qyz = float(sys.argv[3])

square_root = math.sqrt( Qzz*Qzz + 2*Qzz*Qyy + Qyy*Qyy - 4*(Qzz*Qyy-Qzy*Qyz)  )

eigenvalue_1 = (Qzz + Qyy + square_root)/2.0

eigenvalue_2 = (Qzz + Qyy - square_root)/2.0

print ""
print "eigenvalue one is =",eigenvalue_1
print "eigenvalue two is =",eigenvalue_2
print ""
print "These will be the diagonalized values Qzz and Qyy"
print "I will assume that Qzz is larger (is this ok to do?)"
print "It seems that you can rotate it so that the longer Q is along the z anyhow, so ok."

new_Qzz = 0
new_Qyy = 0

if eigenvalue_1 > eigenvalue_2:
    new_Qzz = eigenvalue_1
    new_Qyy = eigenvalue_2
else:
    new_Qzz = eigenvalue_2
    new_Qyy = eigenvalue_1

if new_Qzz == 0 or new_Qyy == 0:
    print ""
    print "Qzz or Qyy is zero! Cant give ratio"
    exit(0)

print ""
print "And the ratio Qzz/Qyy =",new_Qzz/new_Qyy
print ""


# For 18_55-18_55 95
# total cell volume is 4.195 and cm_z = 2.760 and cm_y = 2.213
# Qzz = 38.114  Qyy = 35.470  Qzy = 21.031
# These values incorporate the factor 2 for the Qzz and Qyy and
# the factor 3 for the Qzy, and this is wrong, so without factors I'll use:
# Qzz = 19.057  Qyy = 17.735  Qzy = 7.0103

# For 18_65-28_65 94
# Total cell volume is = 4.123500 and cm_z = 3.284 and cm_y = 1.854
# Qzz = 66.198  Qyy = 19.337  Qzy = 9.173
# These values incorporate the factor 2 for the Qzz and Qyy and
# the factor 3 for the Qzy, and this is wrong, so without factors I'll use:
# Qzz = 33.099  Qyy = 9.668  Qzy = 3.057


# For 18_55-18_55 95
# Qzz/Qyy = 3.858
# but this is with the factors 2 and 3. Answers without those factors:
# Qzz/Qyy = 2.240

# This script gives me answers:
# For 18_65-28_65 94
# Qzz/Qyy = 3.681
# but this is with the factors 2 and 3. Answers without those factors:
# Qzz/Qyy = 3.611

