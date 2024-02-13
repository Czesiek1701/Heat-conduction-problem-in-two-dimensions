# Heat-conduction-problem-in-two-dimensions
Project realised during course: FEM technique. Scripts contain problems of generating mesh, applying boundary conditions, agregation of matrices and presentation of results.

Used libraries:
 - numpy
 - sympy
 - matplotlib

Differential equation solved with Galerking FEM method. Shape function determined by interpolation to nodes values with Vandermonde matrix.

Inputs are file containing geometry parameters and boundary conditions, as sample below:

 - file *geometry.txt* - contains points of opposite corners

`
1 -0.00 0.02 0.06 0.03 # first rectangle `<br>`
2 0.05 0.00 0.06 0.02  # second rectangle
`

 - file *boundary_conditions.txt*

`
1 #cnv - convection, hf - heatflow, T - temperature`<br>`
2 -0.00 0.02 -0.00 0.03 cnv 85 20 `<br>`
3 -0.00 0.02 0.05 0.02 cnv 85 20 `<br>`
4 -0.00 0.03 0.06 0.03 cnv 85 20 `<br>`
5 0.06 0.00 0.06 0.03 T 130 `<br>`
6 0.05 0.00 0.05 0.02 cnv 85 20 `<br>`
7 0.05 0.00 0.06 0.00 cnv 85 20 
`

This input generate results results above:

![edge_size_0_001](https://github.com/Czesiek1701/Heat-conduction-problem-in-two-dimensions/assets/157902583/1ccd8295-e638-4c37-bdfc-f930a4aecae7)



