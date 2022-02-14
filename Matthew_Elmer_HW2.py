################################################################################
# Solar array beta angle optimization
# Matthew Elmer
# 2/13/2022
# Part of an AERO 651 Assignment for Dr. Chamitoff at Texas A&M University
################################################################################


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Note: working in inches
RAD_TO_DEG = 180/np.pi
DEG_TO_RAD = 1/RAD_TO_DEG

NUM_POINTS = 500  # WARNING: this will nuke your computer. try 50 points if you want to manipulate the plot with plt.show()
TOL = 1e-4

A = 175 * 120  # in^2
phi = 1366 / 1550  # W/in^2
eta = 0.11  # between 0.10 and 0.12 for ISS solar arrays

def x(betaB):
    try:  # Try to create an iterator object. If it fails, you do not have an iterable input
        iter(betaB)
    except:  # Since you do NOT have an iterable input, return the results for the single input
        return 225 * np.sin(betaB)
    return np.array([x(betaB[i])] for i in range(len(betaB)))  # Since you DO have an iterable object, compute for all values

def y(betaB):
    try:
        iter(betaB)
    except:
        return 600 - 225 * np.cos(betaB)
    return np.array([y(betaB[i]) for i in range(len(betaB))])

def z(betaB):
    try:
        iter(betaB)
    except:
        return np.sqrt(x(betaB)**2 + y(betaB)**2)
    return np.array([z(betaB[i]) for i in range(len(betaB))])

def s(betaS, betaB):
    try:
        iter(betaB)
    except:
        return (225 - z(betaB) * (np.sin(np.pi/2 - betaS - np.tan(x(betaB)/y(betaB))))
                           / (np.sin(np.pi/2 + betaS - betaB)))
    return np.array([s(betaS[i], betaB[i]) for i in range(len(betaS))])

def S(betaS, betaB):
    try:
        iter(betaB)
    except:
        result = s(betaS, betaB) / 175
        if result > 1:
            return 1
        elif result < 0:
            return 0
        return result
    return np.array([S(betaS[i], betaB[i]) for i in range(len(betaS))])

PN = A * phi * eta

def P4A(betaS, betaB):
    try:
        iter(betaB)
    except:
        if S(betaS, betaB) >= 0.2:
            return 0
        return (1 - S(betaS, betaB)) * PN
    return np.array([P4A(betaS[i], betaB[i]) for i in range(len(betaS))])

def P(betaS, betaB):
    try:
        iter(betaB)
    except:
        return (3 * PN + P4A(betaS, betaB)) * np.cos(betaS - betaB)
        # return P4A(betaS, betaB) * np.cos(betaS - betaB)  # comment the line above and uncomment this line to see
        # what the power output results look like with just the shadowed solar array
    return np.array([P(betaS[i], betaB[i]) for i in range(len(betaS))])

betaSolar = np.linspace(0, 70 * DEG_TO_RAD, NUM_POINTS)  # functions work in radians...
betaBGA = np.linspace(0, 70 * DEG_TO_RAD, NUM_POINTS)
betaSolarmesh, betaBGAmesh = np.meshgrid(betaSolar, betaBGA)
power = P(betaSolarmesh, betaBGAmesh)
betaSolarmesh *= RAD_TO_DEG  # ... but I want to see degrees on the axes
betaBGAmesh *= RAD_TO_DEG

betaBGAOptimal_indices = np.argmax(power, axis=0)  # find indices of max power output
betaBGAOptimal = np.array([betaBGA[idx] for idx in betaBGAOptimal_indices])
fig2D, p2D = plt.subplots(constrained_layout=True)
p2D.plot(betaSolar * RAD_TO_DEG, betaBGAOptimal * RAD_TO_DEG)
p2D.set(
    xlabel="Solar Beta Angle (deg)",
    ylabel="Optimal BGA Angle (deg)",
    title="Optimal BGA Angle vs. Solar Beta Angle"
)
fig2D.savefig("Optimal BGA Angle", dpi=400)

fig3D, ax3D = plt.subplots(subplot_kw={"projection": "3d"})
ax3D.set_xlabel("Solar Beta Angle (deg)")
ax3D.set_ylabel("BGA Angle (deg)")
ax3D.set_zlabel("Power (W)")
surface3D = ax3D.plot_surface(betaSolarmesh, betaBGAmesh, power,
                              cmap=cm.Spectral_r, antialiased=False,
                              rcount=NUM_POINTS, ccount=NUM_POINTS)  # these "-count" params allow plot to be prettier

# It's much easier to set the view and let your computer chug out a three-view,
# and with 500 points it's necessary
ax3D.azim = -110
ax3D.elev = 45
fig3D.savefig("Power Surface Front", dpi=900)

ax3D.azim = 0
ax3D.elev = 30
fig3D.savefig("Power Surface Rear", dpi=900)

ax3D.azim = -90
ax3D.elev = 89.99
fig3D.savefig("Power Surface Top", dpi=900)
