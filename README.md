# Non-Linear Mechanical System Control: Linearization & Loop Shaping

This project focuses on the modeling and control of a non-linear motorized mechanism. It involves the transition from non-linear differential equations to a robust linear controller capable of meeting strict dynamic and static specifications.

##  Project Overview
The system is characterized by a position-dependent moment of inertia $J(\theta)$, creating a non-linear dynamic. The goal was to design a controller that ensures stability, noise rejection, and precise tracking even in the presence of disturbances.

##  Engineering Pipeline
1. **State-Space Modeling:** Defined the system using state variables ($\theta, \omega$) and identified the equilibrium point $(x_e, u_e)$.
2. **Linearization:** Derived the linear approximation via Jacobian matrices ($A, B, C, D$) to obtain the Transfer Function $G(s)$.
3. **Loop Shaping Design:** - **Static Regulator:** Designed to nullify steady-state error and reject low-frequency disturbances ($50\text{ dB}$ abatement).
   - **Dynamic Regulator:** Implemented a **Lead Network (Rete Anticipatrice)** to guarantee a Phase Margin $M_f \geq 45^\circ$ and meet overshoot constraints ($S\% \leq 8\%$).
4. **Stability Analysis:** Verified robust stability using Bode diagrams and the **Roofline** of the loop function $L(s)$.



##  Performance & Validation
- **Linear Tests:** Successfully met settling time ($T_{a, \epsilon} < 0.1\text{ s}$) and overshoot specs.
- **Non-Linear Validation:** Simulated the controller on the original non-linear model in **Simulink**, testing for asymptotic stability and basin of attraction across various initial conditions.
- **Disturbance Rejection:** High-frequency noise ($n(t)$) was attenuated by over $60\text{ dB}$ through high-frequency pole insertion.

##  Tech Stack
- **Tools:** MATLAB, Simulink.
- **Theory:** State-space representation, Frequency response, Lead-lag compensation, Non-linear dynamics.
