# CellListMapArticleCodes

Code blocks of the CellListMap paper explained. 

The package is available at: [http://github.com/m3g/CellListMap.jl](http://github.com/m3g/CellListMap.jl) 


Article: Martínez, L. *CellListMap.jl: Efficient and customizable cell list implementation for calculation of pairwise particle properties within a cutoff.* https://doi.org/10.48550/arXiv.2202.06427


<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock1.jl.html>Code Block 1</a><br>
Minimal working example for the computation of the sum of distances of 100k particles in a cubic periodic box of side 1 and a cutoff of 0.05.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock2.jl.html>Code Block 2</a><br>
General format of the function to be evaluated for each pair of particles closer than the cutoff distance, to be passed to the map_pairwise function.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock3.jl.html>Code Block 3</a><br>
Initialization of the system Box, with orthorhombic or triclinic periodic boundary conditions. The system’s geometry is defined by the type of unit cell matrix, and an orthorhombic cell is assumed if a vector of box sides is supplied.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock4.jl.html>Code Block 4</a><br>
Computing the cell lists from the coordinates, x, and the system box. Particles are replicated at the boundaries to avoid coordinate wrapping in the function mapping step.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock5.jl.html>Code Block 5</a><br>
Mapping the computation of a property into the pairs of particles which, according to the cell lists and system box properties, are within the desired cutoff.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock6.jl.html>Code Block 6</a><br>
Calculation of a simple Lennard-Jones potential energy of 3 million Neon atoms in 3 dimensions, in a periodic cubic box with sides of 31.034 nm and a cutoff of 1.2 nm.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock7.jl.html>Code Block 7</a><br>
Example code for the calculation of a vector of forces between particles. The function will update the `f` vector.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock8.jl.html>Code Block 8</a><br>
Example code for the calculation of a nearest-neighbor list between two independent sets of particles. A custom reduction function is required to merge lists, keeping the minimum distances.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock9.jl.html>Code Block 9</a><br>
Units, uncertainties, and automatic differentiation propagating through pairwise computations with CellListMap.jl.

<a href=http://m3g.github.io/CellListMapArticleCodes/CodeBlock10.jl.html>Code Block 10</a><br>
Computing a histogram of average pairwise velocities between galaxies, as a function of their relative distances, a typical calculation in astrophysical simulations.

