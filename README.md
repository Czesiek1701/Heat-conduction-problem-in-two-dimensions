# Heat-conduction-problem-in-two-dimensions
Project realised during course: FEM technique. Scripts contain problems of generating mesh, applying boundary conditions, agregation of matrices and presentation of results.

Used libraries:
 - numpy
 - sympy
 - matplotlib

Differential equation solved with Galerking FEM method. Shape function determined by interpolation to nodes values with Vandermonde matrix.

Inputs are file containing geometry parameters and boundary conditions, as sample below:
file geometry.txt
1 -0.00 0.02 0.06 0.03
2 0.05 0.00 0.06 0.02



