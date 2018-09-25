# 2d-ising-hexagonal

Notes:
How to choose ntrans and nmcs?
At low temperatures, it takes a long time to reach equilibrium state.
The system may trap in local minimum states for a long time.
For example, if I set: J1=5.0, J2=J3=0, kb=1.0, size=25, and T=1.0, it takes about 70000 steps to reach global minimum state
Maybe implement cluster-flip dynamics.
