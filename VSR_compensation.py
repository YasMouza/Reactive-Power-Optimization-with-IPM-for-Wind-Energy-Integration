### VSR Kompensation

### This code does show the Compensation in alpha ventus grid with the VSR as a single compensating element
### Here, there is one value for P on busbar 5 that is definded in advance



# import
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import cmath
from pyomo.environ import *
import matplotlib.pyplot as plt
from pyomo.environ import cos, sin
from pyomo.environ import SolverFactory
solver = SolverFactory('ipopt', executable='C:\\Users\\User\\AppData\\Local\\Ipopt-3.11.1-win64-intel13.1\\bin\\ipopt.exe')



results_list = []

Vpcc = 1.05
model = ConcreteModel()   
n=5
# Variablen
model.V = Var(range(n), within=Reals, bounds=(0.98, 1.1))  # Spannungsbegrenzungen
model.theta = Var(range(n), within=Reals, bounds = (-np.pi, np.pi))
model.delta = Var(range(n), within = Reals, bounds = (-math.pi, math.pi))
model.P = Var(range(n), within = Reals, bounds=(-1.5, 1.5))
model.P[4].setlb(None) 
model.P[4].setub(None) 
model.Q = Var(range(n), within = Reals, bounds =(-1,1))

# operational range of VSR (10 steps)
model.k = Var(bounds = (0, 1))
# auxiliary variable for initializing permissible offset
model.deltaV = Var()

# permissible offset determines the "area" around 1.05 p.u. 
model.permissible_offset_plus = Var(within=Reals, bounds=(0, 0.5))
model.permissible_offset_minus = Var(within=Reals, bounds=(0, 0.5))

# P5 = infeed busbar
model.P[4].fix(1) ## adjust x in fix(x) in order to observe optmization behaviour

# VSR is only switching equipment: Onshore and Offshore Transformer are set as constant
ton = 0
toff = 1

# constants
VINF = 220e3
VON = 110e3
VWP = 30e3
PI = np.pi
F = 50
OMEGA = 2*PI*F

# Base quantities
VB = VINF
SB = 100e6
ZB = VB**2 / SB
YB = 1 / ZB


# Infitine Bus
Skpp = 2000e6   #strong grid equivalent
#Skpp = 400e6   #weak grid equivalent

# SI
Xinf = VINF**2 / Skpp
Rinf = 0.1 * Xinf

# pu
xinf = Xinf / ZB
rinf = Rinf/ ZB
zinf = rinf + 1j*xinf

# Onshore Transformer
xtrf_on_pu = 0.15
Strf_on = 75e6
# pu
xtrf_on = xtrf_on_pu * (SB/Strf_on)
ztrf_on = 1j * xtrf_on
# PI circuit
ytrf_on = 1 / ztrf_on
dv_on = 1.25 / 100
theta_on = (1 + ton*dv_on)
#theta_on = (1 + model.ton*dv_on)
theta_on

# Offshore Transformer
xtrf_off_pu = 0.13
Strf_off = 75e6
# pu
xtrf_off = xtrf_off_pu * (SB/Strf_on)
ztrf_off = 1j * xtrf_off
# PI circuit
ytrf_off = 1 / ztrf_off
dv_off = 1.25 / 100
theta_off = (1 + toff*dv_off)
#theta_off = (1+model.toff*dv_off)

# VSR
Qvsr = 35e6

# SI
bvsr = (Qvsr/SB) * (VB/VON)
#yvsr = (-1) * 1j * scale*bvsr
yvsr = (-1) * model.k*bvsr

# FSR
Qfsr = 35e6

# SI
bfsr = (Qfsr/SB) * (VB/VON)
yfsr = (-1)*1j*bfsr

# Cable
Ccable_km = 280e-9
Rcable_km = 0.018
Xcable_km = 0.123

# SI
length = 70
Ccable = Ccable_km * length
Rcable = Rcable_km * length
Xcable = Xcable_km * length
Bcable = OMEGA * Ccable

# pu
rcable = Rcable * (VB/VON)**2 / ZB
xcable = Xcable * (VB/VON)**2 / ZB
bcable = Bcable / (VB/VON)**2 * ZB

zcable = rcable + 1j * xcable
ycable = 0.5 * 1j * bcable
zcable, ycable

# model.Q[4].fix(0)

model.delta[0].fix(0)

model.V[1].setlb(1.01)  # Untere Grenze für die PCC-Spannung
model.V[1].setub(1.1)   # Obere Grenze für die PCC-Spannung

b01 = b10 = xinf/(xinf**2+rinf**2)
ytrf_on = 1/ xtrf_on
b12 = b21 = ytrf_on/theta_on
b23 = b32 = xcable/(xcable**2+rcable**2)
ytrf_off = 1/xtrf_off
b34 = b43 = ytrf_off/theta_off

yvsr = (-1) * model.k*bvsr
ycable = 0.5 * bcable
yfsr = (-1)*bfsr


# Shunt Admittances
y1g = (1 - theta_on)/(theta_on**2) * ytrf_on
y2g = (theta_on - 1)/theta_on * ytrf_on + yvsr + ycable
y3g = (1 - theta_off)/(theta_off**2) * ytrf_off + ycable + yfsr
y4g = (theta_off - 1)/theta_off * ytrf_off

y1g, y2g, y3g, y4g

b00 = -b01
b11 = -(b01+b12+y1g)
b22 = -(b12+b23+y2g)
b33 = -(b23+b34+y3g)
b44 = -(b34+y4g)

g01 = g10 = -rinf/(rinf**2+xinf**2)
g12 = g21 = 0
g23 = g32 = -rcable/(rcable**2+xcable**2)
g34 = g43 = 0

g00 = -(g01)
g11 = -(g01+g12)
g22 = -(g12+g23)
g33 = -(g23+g34)
g44 = -(g34)


# Constraint für Knotenpunkt 1
model.active_power_constraint_1 = Constraint(expr=
    model.V[0] * model.V[1] * (g01 * cos(model.delta[0] - model.delta[1]) + b01 * sin(model.delta[0] - model.delta[1])) +
    model.V[0] * model.V[0] * (g00 * cos(model.delta[0] - model.delta[0]) + b00 * sin(model.delta[0] - model.delta[0])) == model.P[0]
)

# Constraint für Knotenpunkt 2
model.active_power_constraint_2 = Constraint(expr=
    model.V[1] * model.V[0] * (g10 * cos(model.delta[1] - model.delta[0]) + b10 * sin(model.delta[1] - model.delta[0])) +
    model.V[1] * model.V[2] * (g12 * cos(model.delta[1] - model.delta[2]) + b12 * sin(model.delta[1] - model.delta[2])) +
    model.V[1] * model.V[1] * (g11 * cos(model.delta[1] - model.delta[1]) + b11 * sin(model.delta[1] - model.delta[1])) == model.P[1]
)

# Constraint für Knotenpunkt 3
model.active_power_constraint_3 = Constraint(expr=
    model.V[2] * model.V[1] * (g21 * cos(model.delta[2] - model.delta[1]) + b21 * sin(model.delta[2] - model.delta[1])) +
    model.V[2] * model.V[3] * (g23 * cos(model.delta[2] - model.delta[3]) + b23 * sin(model.delta[2] - model.delta[3])) + 
    model.V[2] * model.V[2] * (g22 * cos(model.delta[2] - model.delta[2]) + b22 * sin(model.delta[2] - model.delta[2]))== model.P[2]
)

# Constraint für Knotenpunkt 4
model.active_power_constraint_4 = Constraint(expr=
    model.V[3] * model.V[2] * (g32 * cos(model.delta[3] - model.delta[2]) + b32 * sin(model.delta[3] - model.delta[2])) +
    model.V[3] * model.V[4] * (g34 * cos(model.delta[3] - model.delta[4]) + b34 * sin(model.delta[3] - model.delta[4])) +
    model.V[3] * model.V[3] * (g33 * cos(model.delta[3] - model.delta[3]) + b33 * sin(model.delta[3] - model.delta[3])) == model.P[3]
)

# Constraint für Knotenpunkt 5
model.active_power_constraint_5 = Constraint(expr=
    model.V[4] * model.V[3] * (g43 * cos(model.delta[4] - model.delta[3]) + b43 * sin(model.delta[4] - model.delta[3])) +
    model.V[4] * model.V[4] * (g44 * cos(model.delta[4] - model.delta[4]) + b44 * sin(model.delta[4] - model.delta[4])) == model.P[4]
)


# Constraint für Knotenpunkt 1
model.reactive_power_constraint_1 = Constraint(expr=
    model.V[0] * model.V[1] * (g01 * sin(model.delta[0] - model.delta[1]) - b01 * cos(model.delta[0] - model.delta[1])) +
    model.V[0] * model.V[0] * (g00 * sin(model.delta[0] - model.delta[0]) - b00 * cos(model.delta[0] - model.delta[0])) == model.Q[0]
)

# Constraint für Knotenpunkt 2
model.reactive_power_constraint_2 = Constraint(expr=
    model.V[1] * model.V[0] * (g10 * sin(model.delta[1] - model.delta[0]) - b10 * cos(model.delta[1] - model.delta[0])) +
    model.V[1] * model.V[2] * (g12 * sin(model.delta[1] - model.delta[2]) - b12 * cos(model.delta[1] - model.delta[2])) +
    model.V[1] * model.V[1] * (g11 * sin(model.delta[1] - model.delta[1]) - b11 * cos(model.delta[1] - model.delta[1])) == model.Q[1]
)

# Constraint für Knotenpunkt 3
model.reactive_power_constraint_3 = Constraint(expr=
    model.V[2] * model.V[1] * (g21 * sin(model.delta[2] - model.delta[1]) - b21 * cos(model.delta[2] - model.delta[1])) +
    model.V[2] * model.V[3] * (g23 * sin(model.delta[2] - model.delta[3]) - b23 * cos(model.delta[2] - model.delta[3])) + 
    model.V[2] * model.V[2] * (g22 * sin(model.delta[2] - model.delta[2]) - b22 * cos(model.delta[2] - model.delta[2])) == model.Q[2]
)

# Constraint für Knotenpunkt 4
model.reactive_power_constraint_4 = Constraint(expr=
    model.V[3] * model.V[2] * (g32 * sin(model.delta[3] - model.delta[2]) - b32 * cos(model.delta[3] - model.delta[2])) +
    model.V[3] * model.V[4] * (g34 * sin(model.delta[3] - model.delta[4]) - b34 * cos(model.delta[3] - model.delta[4])) +
    model.V[3] * model.V[3] * (g33 * sin(model.delta[3] - model.delta[3]) - b33 * cos(model.delta[3] - model.delta[3])) == model.Q[3]
)

# Constraint für Knotenpunkt 5
model.reactive_power_constraint_5 = Constraint(expr=
    model.V[4] * model.V[3] * (g43 * sin(model.delta[4] - model.delta[3]) - b43 * cos(model.delta[4] - model.delta[3])) +
    model.V[4] * model.V[4] * (g44 * sin(model.delta[4] - model.delta[4]) - b44 * cos(model.delta[4] - model.delta[4])) == model.Q[4]
)

# formulate objective function as function in dependency of 1.05 and permissible offset
model.abs_constraint1 = Constraint(expr=model.deltaV <= model.V[1] - 1.05 + model.permissible_offset_plus - model.permissible_offset_minus)
model.abs_constraint2 = Constraint(expr=model.deltaV >= -(model.V[1] - 1.05 + model.permissible_offset_plus - model.permissible_offset_minus))

# expr=1 => find operation point that satisfies permissible offset
model.objective = Objective(expr=1)

# access fileshare wheree ipopt is located 
solver = SolverFactory('ipopt', executable='C:\\Users\\User\\AppData\\Local\\Ipopt-3.11.1-win64-intel13.1\\bin\\ipopt.exe')
solver.solve(model)

P_values = [model.P[i].value for i in range(5)]
Q_values = [model.Q[i].value for i in range(5)]
V_values = [model.V[i].value for i in range(5)]
theta_values = [model.delta[i].value for i in range(5)]

# Creating a DataFrame
results_df = pd.DataFrame({
    'Busbar': list(range(1, 6)),
    'P': P_values,
    'Q': Q_values,
    'V': V_values,
    'Theta': theta_values
})

# Displaying the DataFrame
print(results_df)

# Displaying model.k.value separately
print("Model k value:", model.k.value)

# shall not be "within reals" but rounded to 0.i
print("Model k rounded:", round(model.k.value,1))    
