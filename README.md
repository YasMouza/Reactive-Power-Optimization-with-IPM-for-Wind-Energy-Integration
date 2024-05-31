# Reactive-Power-Optimization-with-IPM-for-Wind-Energy-Integration
This repository contains code that optimizes the voltage at a target busbar with a given set of 
parameters by a test grid. The optimization is done via python.ipopt

3 different compensation operations are being carried out: 
    - VSR_compensation.py
    - ton_toff_compensation.py
    - VSR_ton_toff_compensation.py


1. VSR_compensation.py
A Variable Shunt Reactor (VSR) is the only switching equipment in the test grid. 

2. ton_toff_compensation.py
Onshore and Offshore Transformer are compensating equipment in test grid. VSR is deactivated

3. VSR_ton_toff_compensation.py
Onshore, Offshore Transformer and VSR are compensating equipment in test grid.
