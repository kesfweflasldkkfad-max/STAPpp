# Q4 Patch Test Report

## Problem

For the plane stress case,

- E = 1000
- nu = 0.3
- u(x, y) = 0.01 x
- v(x, y) = -0.003 y

The exact small-strain field is

- epsilon_x = du/dx = 0.01
- epsilon_y = dv/dy = -0.003
- gamma_xy = du/dy + dv/dx = 0

For plane stress,

sigma = D epsilon,

where

D = E / (1 - nu^2) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]]

Substituting the strain field gives

- sigma_x = 10
- sigma_y = 0
- tau_xy = 0

So this is a constant-stress patch-test field with zero body force.

## What Was Verified

The current STAP++ repository only implements the Bar element. Q4 is listed in the element-type enum, but there is no Q4 element class, no Q4 material class, and no Q4 data example. To complete the assignment, a standalone script was added:

- STAPpy test.py

The script builds a distorted 2x2 mesh of four bilinear isoparametric Q4 elements. Exact displacements are prescribed on the boundary nodes, leaving the center node free. With standard 2x2 Gauss integration, the solved center displacement and the Gauss-point stresses are checked against the exact solution.

## Numerical Result

Script output:

```text
Q4 patch test on a distorted 2x2 mesh
E = 1000.0, nu = 0.3
Exact strain = [ 0.01  -0.003  0.   ]
Expected stress = [1.0000000e+01 1.5335718e-17 0.0000000e+00]
Solved center displacement = [ 0.01  -0.003]
Exact center displacement = [ 0.01  -0.003]
Max center displacement error = 1.734723e-18
Max Gauss-point stress error = 5.329071e-15
Max free-DOF residual = 8.881784e-16
Patch test passed.
```

## Conclusion

The bilinear 4-node isoparametric quadrilateral with 2x2 Gauss integration passes the patch test for the displacement field

- u(x, y) = 0.01 x
- v(x, y) = -0.003 y

and reproduces the constant plane-stress state

- sigma_x = 10
- sigma_y = 0
- tau_xy = 0

to machine precision.