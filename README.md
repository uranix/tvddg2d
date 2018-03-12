# tvddg2d
General purpose high-order Discontinunous Galerkin 2D solver for hyperbolic conservation laws problems.

Features:
* Solves homogenous hyperbolic systems with any number of equations. User needs to define the variable vector
and the flux computation routines
* Riemmann problem based. Used needs to provide an (approximate) generalized Riemmann solver. Generalized means that not
only state variables have a discontinuity but also may have the parameters. If parameters are constant then this problem
reduces to regular Riemmann problem.
* High order DG approximation. Each cell is subdivided into several subcells (up to 6x6) to raise the order of the spatial discretization
* Strong stability preserving Runge-Kutta integrators (orders 1-3).
* The scheme is monotonized with a flux limiter technique and is total variation diminishing. It is nessesary to provide routines to
define the complete eigendecomposition for the Jacobi matrix of the system.

Here's an example of a droplet simulation inside a basin (left - low order, right - high order TVD). A version without limiting is 
not shown since it breaks when water height is negative due to spurious numerical oscillations.
![LO vs TVD](https://i.imgur.com/oOXtDBI.png)
