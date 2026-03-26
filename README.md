This is a simple model for the gas-cooled reactor. It is used to have a working system.

The model is taken from: 
[[Application of the Method of Manufactured Solutions to a Close-Coupled Gas-Cooled Reactor and Brayton Cycle Power System.pdf]]
and
[[Transient simulation of the USNC Pylon using a rapid multiphysics model.pdf]]

An extra battery model is added. 

![[Pasted image 20260325142648.png]]

And the turbine side is simplified. 

![[Pasted image 20260325160529.png]]
# Reactor

## Reactor Kinetics

The power of the reactor is modeled with the PKE:

$$
\frac{d P(t)}{dt} = \frac{\rho-\beta}{\Lambda} P(t) + \sum_i \lambda_i C_i 
$$

With the corresponding DNP Groups:

$$
\frac{d C_i(t)}{dt} = \frac{\beta}{\Lambda} P(t) -  \lambda_i C_i 
$$

The reactivity is calculated from temperature feedback, Xenon feedback, and the drum configuration.
so that 

$$\rho = \rho_{cr}+\rho_{T}+\rho_{Xe} $$

with

$$
\rho_{T} = - \alpha (T(t)_{reactor}-T_{ref})
$$

With the reactor temperature $T(t)_{reactor} = (T(t)_{Rin}+T(t)_{Rout})/2$ and the reference temperature $T_{ref}$.

$$
\rho_{Xe} = - \alpha \frac{Xe \sigma^{Xe}_a}{\Sigma_a}
$$

and the control rods being Sigmoidal:

$$
\rho_{cr} = 4000/(1+exp(-x10))-2000
$$

 with the steady state at the control rod position $x = 0$


The Xenon concentration $Xe$ is calculated using the flux

$$
\Phi = \frac{P}{\Sigma_f Q_{fission}}
$$

and the usual:

$$
\frac{d I(t)}{dt} = \gamma_I \Sigma_f \Phi -\lambda_i I(t)
$$

$$
\frac{d Xe(t)}{dt} = \gamma_X \Sigma_f \Phi +\lambda_i I(t) - \lambda_x Xe(t)- \sigma_a^X \Phi Xe(t)
$$

## Temperature balance 

The coolant is assumed to be incompressible. In a more advanced model, it should be adjusted.

The entire power of the reactor is assumed to be dumped into the coolant gas. At some point, it might be interesting to test radiative cooling.

In this model, the change in coolant temperature:

$$
\Delta T = \frac{P(t)}{C_R \dot m_R}
$$

with the specific heat capacity of the reactor coolant of $c_R$ and the reactor coolant mass flow rate $\dot m_R$.

So that the temperature at point $R_{out}$
$$
T_{Rout} = T_{Rin} +\Delta T = T_{Rin} + \frac{P(t)}{C_R \dot m_R}
$$


# Battery

The battery is a thermal battery. The exact composition is unknown, but scenarios like loss of battery fluid could be simulated in the future.

The change in temperature is given by

$$
C_B \frac{d T_B}{dt} = \dot Q_{in}-\dot Q_{out}
$$

With the heat capacity of the battery $C_b$

The heat flow into the battery is the difference of the in and outlet temperature (incopressible). So that 

$$
\dot Q_{in} = \dot m_R c_R (T_{Rout}-T_{Rin})
$$

And

$$
\dot Q_{out} = \dot m_T c_T (T_{Bout}-T_{Bin})
$$

with the specific heat capacity of the turbine working fluid of $c_T$ and the reactor coolant mass flow rate $\dot m_T$. 


With the NTU method, we estimate the heat transfer from the gas to the battery. 

$$
T_{Rin} = T_B +(T_{Rout}-T_B) \exp\left[-\frac{A_R}{\dot m_R c_R}\right]
$$

In the same way, we calculate the outlet temperature on the Turbine side.

$$
T_{Bout} = T_B +(T_{Bin}-T_B) \exp\left[-\frac{A_T}{\dot m_T c_T}\right]
$$

# Turbine

We assume that the efficiency is 50% of the Carnot efficiency, with the final heat sink temperature of $293.15$ K.

$$
\eta = 0.5 (1-\frac{T_c}{T_h})
$$

In this case $T_c = T_{B}$ and $T_h = T_{Bout}$.
