# ? Comments that begin with a question mark are hints to students using
# ? this code as a starting point for calculating something different
import math
import numpy as np
from scipy.optimize import root
#from hL_gelatine import hL_gelatine #? Delete this line when working with CT2
from H_steam import H_steam
from hL_black_liquor import hL_black_liquor
from BPE import BPE
from k_black_liquor import k_black_liquor


#def Ts_from_Ps_water(Ps):
#    # ? Note: In task 5.7 we were given the steam pressure
#    # ? You wont be needing Antoines equation in the compulsory task
#    # Calculation of vapor temperature (saturated steam)
#    Ts=1730.630/(-math.log10(Ps)+10.19625)-233.426 # Antoine equation, see handbook
#    return Ts


def evaporator(X, Ts, Tc, F, Tf, xF, xL1, k_constant,use_BPE,vary_k):
    #Calculates residuals of mass and energy balances for an evaporator
    #
    #    Input: X contains guesses of operating conditions in this order
    #       condition      unit     explanation
    #       S,             kg/s     Steam flux
    #       V,             kg/s     Vapor flux (pure water)
    #       L,             kg/s     Liquid flux (pure water)
    #       A               m2      Heat exchanger area
    #
    #    Other input (known parameters)
    #       F              kg/s     Feed flow
    #       xF             kg/kg    Solids fraction in feed
    #       Tf             Celsius  Feed temperature
    #       Ts             Celcius  Fresh steam temperature
    #       TL             Celsius  Product stream temperature
    #       xL             kg/kg    Solids fraction in product stream
    #       k_constant     J/m2,K   Overall heat transfer coefficient
    #       use_BPE        boolean  Take boiling point elevation into account?
    #       vary_k         boolean  Let k depend on e.g. T, x?
    #
    #    Output: residuals that should be zero if all data in X are
    #    internally consistent
    #
    #    Assumptions:
    #       * xV=0, i.e. vapor flow from evaporator is pure water
    #       * heat capacity for water is constant within the temperature range
    #       * involatile substance is a liquid with a constant heat capacity
    #       * no mixing enthalpy, i.e. enthalpy of flows containing water and
    #         the involatile substance can be calculated using heat capacities
    #         for pure water and pure involatile substance



    # ================ Initialisation ========================= 
    #
    # Copy X-values to variables that are easier to understand
    [L1,L2,L3,V1,V2,V3,S,xL2,xL3,A,T1,T2]=X

    # =============== Calculations start here ==================
    #
    if use_BPE:
        b1=BPE(xL2)
        b2=BPE(xL3)
        b3=BPE(xF)
    else:
        b1=0
        b2=0
        b3=0

    T3=Tc+b3

    # Calculation of water vapor enthalpies 
    # ? NOTE: To calculate H for overheated steam, 2nd argument cannot be -1
    [Hs,dummy1,dummy2]=H_steam(Ts,-1) 
    [HV1,dummy1,dummy2]=H_steam(T1,-1)
    [HV2,dummy1,dummy2]=H_steam(T2,-1)
    [HV3,dummy1,dummy2]=H_steam(T3,-1)

    hk1=hL_black_liquor(0,Ts)
    hL1=hL_black_liquor(xL1,T1)
    hk2=hL_black_liquor(0,T1)
    hL2=hL_black_liquor(xL2,T2)
    hk3=hL_black_liquor(0,T2)
    hL3=hL_black_liquor(xL3,T3)
    hF =hL_black_liquor(xF,Tf)
    
        
    if vary_k:
        [k1,dummy]=k_black_liquor(xL2, xL1, Ts, L2, L1)
        [k2,dummy]=k_black_liquor(xL3, xL2, T1, L3, L2)
        [k3,dummy]=k_black_liquor(xF,  xL3, T2, F,  L3)
        # ? Delete line above and write some code here to implement
        # ? handling of boiling point elevation
    else:
        k1=k_constant
        k2=k_constant
        k3=k_constant

    # ================ Calculating residuals: ==================
    #
    # 
    Y=X*0
    # evaporator 1
    Y[0]=L2-L1-V1                           # MB total
    Y[1]=xL2*L2-xL1*L1                      # MB solids
    Y[2]=S*(Hs-hk1)+L2*hL2-L1*hL1-V1*HV1    # EB over evaporator
    Y[3]=S*(Hs-hk1)-k1*A*(Ts-T1)             # EB over heat exchanger


    # evaporator 2
    Y[4]=L3-L2-V2                           # MB total
    Y[5]=xL3*L3-xL2*L2                      # MB solids
    Y[6]=V1*(HV1-hk2)+L3*hL3-L2*hL2-V2*HV2  # EB over evaporator
    Y[7]=V1*(HV1-hk2)-k2*A*(T1-b1-T2)           # EB over heat exchanger


    # evaporator 3
    Y[8]=F-L3-V3                            # MB total
    Y[9]=xF*F-xL3*L3                        # MB solids
    Y[10]=V2*(HV2-hk3)+F*hF-L3*hL3-V3*HV3   # EB over evaporator
    Y[11]=V2*(HV2-hk3)-k3*A*(T2-b2-T3)          # EB over heat exchanger

    return Y


Ts = 158.83         # degrees C               Temperature fresh steam
Tc = 45.81          # degrees C               Temperature Condenser
F  = 15             # kg/s                    Feed flux
Tf = 40             # degrees C               Feed temperature
xF = 0.25           # kg dry matter/kg total  Feed dry matter content
xL1= 0.8            # kg dry matter/kg total  Product dry matter content
k_constant = 0.855  # kW/m2/K                 Overall heat transfer coeff
use_BPE=True
vary_k=True

known=(Ts, Tc, F, Tf, xF, xL1, k_constant,use_BPE,vary_k)

T_guess=90 #degrees C
flow_guess=7 #kg/s
feed_fraction_guess = (xF+xL1)/2
# Unknowns: L1,L2,L3,V1,V2,V3,S,xL2,xL3,A,T1,T2
guess=np.array([*[flow_guess]*7,feed_fraction_guess,feed_fraction_guess,200,T_guess,T_guess]) #Set initial guesses

sol = root(evaporator, guess, args=known, method='hybr')

if not sol.success:
    print('Iteration not successful:',sol.message)
else:
    print('Iteration successful:',sol.message)
# Translate results to variable names easier to understand
[L1,L2,L3,V1,V2,V3,S,xL2,xL3,A,T1,T2]=sol.x[:]
print(f"L1: {L1}")
print(f"L2: {L2}")
print(f"L3: {L3}")
print(f"V1: {V1}")
print(f"V2: {V2}")
print(f"V3: {V3}")
print(f"S: {S}")
print(f"xL2: {xL2}")
print(f"xL3: {xL3}")
print(f"A: {A}")
print(f"T1: {T1}")
print(f"T2: {T2}")
print(f"k1: {k_black_liquor(xL2, xL1, Ts, L2, L1)[0]}")
print(f"k2: {k_black_liquor(xL3, xL2, T1, L3, L2)[0]}")
print(f"k3: {k_black_liquor(xF,  xL3, T2, F,  L3)[0]}")
print(f"Steam eco: {V3/S}")




# You can use the command "print" to print results
# print('Steam flux',S,'kg/s' )
# print('Vapor flux',V,'kg/s' )
# print('Liquid flux',L,'kg/s' )
# print('Area',A,'m2' )
# print('S/Vtot',S/V)
# print('Vtot/S',V/S)

# Recalculate other results
# ? Note: A common beginners method to retrieve these kind of values is
# ? to print out values inside the routine called by root
# ? There are two problems why that should be avoided:
# ? 1. There is in general no guarante that the solver, in its very
# ?    last call of the function, will use the correct solution
# ? 2. Printing out things to the screen every time a function is called
# ?    by the solver slows down the calculation considerably
#Ts=Ts_from_Ps_water(Ps)
#[Hs,P,dummy]=H_steam(Ts,-1)
#[HV,P,dummy]=H_steam(T,-1)
#print('Ts',Ts,'degree C')
#print('Steam enthalpy Hs',Hs,'kJ/kg')
#
##? Variables in Python can easily be printed to files
##? Note: the "w" argument WRITEs over the file if it already exists
##? If you instead want to APPEND your new results to end of file, change to "a"
#with open("output.txt", "w") as f:
#  print(S,V,L,A,Ts,Hs,HV, file=f)
