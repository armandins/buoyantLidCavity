# Lid Driven Cavity with Buoyancy Effects Solution via Stream-function Vorticity Approach
Lid driven cavity problem solved using stream function vorticity approach. Buoyant forces are considered and energy equation is also solved on top of transport and continuity equations. The final form of the nondimensionalized equations which are solved for here are:

$$
\frac{\partial U}{\partial X} + \frac{\partial V}{\partial Y} = 0
$$

$$
\frac{\partial U}{\partial \tau} + U\frac{\partial U}{\partial X} + V\frac{\partial U}{\partial Y} = -\frac{\partial P^\star}{\partial X} + \frac{1}{Re}(\frac{\partial^2U}{\partial X^2} + \frac{\partial^2 U}{\partial Y^2})
$$

$$
\frac{\partial V}{\partial \tau} + U\frac{\partial V}{\partial X} + V\frac{\partial V}{\partial Y} = -\frac{\partial P^\star}{\partial Y} + \frac{1}{Re}(\frac{\partial^2V}{\partial X^2} + \frac{\partial^2 V}{\partial Y^2}) + Ri\Theta
$$

$$
\frac{\partial \Theta}{\partial \tau} + U\frac{\partial \Theta}{\partial X} + V\frac{\partial \Theta}{\partial Y} = \frac{1}{Pe}(\frac{\partial^2\Theta}{\partial X^2} + \frac{\partial^2\Theta}{\partial Y^2})
$$  

Where $\Theta = \frac{T - T_c}{T_h - T_c}$ , $\tau = \frac{t U_0}{l}$ , $X = \frac{x}{l}$ , $Y = \frac{y}{l}$ , $U = \frac{u}{u_0}$ , $V= \frac{v}{u_0}$ and $P^{\star} = \frac{p}{\rho (u_0)^2}$. 
## Boundary Conditions
Boundary conditions are $U = 1$ for top wall and $U = 0$ for the rest. $\Theta_l = 1$ or $T = T_h$ for left wall, $\Theta_r = 0$ or $T = T_c$ for right wall and bottom and top walls have zero gradient conditions. 

## Results

$Ri = \frac{Gr}{Re^2}$ and $Pe = RePr$. Gr stands for Grashof number which is defined by: 

$$
Gr = \frac{g\beta L^3 (T_h - T_c)}{\nu ^2}
$$

Also for Reynolds number and Prandtl we have the definitions stated below: 

$$
Re = \frac{U_0 L}{\nu}
$$ 
$$
Pr = \frac{\nu}{\alpha}
$$ 
 
The results are written to a .VTK format and can be post-processed via Paraview. To run the code simply use: 
```
g++ -ggdb3 -O1 -std=c++23 -Wall -Wextra -pedantic -o sfvEnergyCavity.out sfvEnergyCavity.cc
```

<div align="center">
    <img src="images/thetaa.png" alt="Alt text" width="300" />
    <p style="font-size: 14px; color: #555;">Temperature variations at Re = 100 and Gr = 1000 Case</p>
</div>


<div align="center">
    <img src="images/sf.png" alt="Alt text" width="300" />
    <p style="font-size: 14px; color: #555;">Stream functions at Re = 100 and Gr = 1000 Case</p>
</div>



