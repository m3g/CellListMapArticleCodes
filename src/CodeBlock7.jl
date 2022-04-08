### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 562abb8e-37ed-4fac-8293-c47e7d1d6203
using CellListMap

# ╔═╡ 56744cec-3cb1-45a1-b681-c99f0ae2eaae
using FastPow

# ╔═╡ c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
md"
### Code Blocks of the CellListMap.jl paper explained
"

# ╔═╡ 7dbfadf4-2cb8-4145-948f-df1e7ea06b14
md"
# Code Block 7
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ dc0e8a10-30ba-44e5-be65-61b0402ed2ed
md"
This example involves the computation of Lennard-Jones forces. Many things are similar to Code Block 6, so we recommend reading that one first. We use the `FastPow` to accelerate the computation of high powers:
"

# ╔═╡ 8d685c45-f9d7-42fd-a6e5-28361d4a7de3
md"
The Lennard-Jones **force** between a pair of particles has the form

$\vec{f}_{ij} = -12\varepsilon\left(
\frac{\sigma^{12}}{d^{13}} - \frac{\sigma^6}{d^7}
\right)\frac{\vec{r}}{d}$

where `\vec{r}` is the vector connecting the two particles, and $d = ||\vec{r}||$ is the distance between them. Considering the fact that the interface in `CellListMap` provides the squared distance, we can rewrite this equation as

$\vec{f}_{ij} = -12\varepsilon\left(
\frac{\sigma^{12}}{(d^2)^{7}} - \frac{\sigma^6}{(d^2)^4}
\right)\vec{r}$

This computes the force vector $\vec{f}_{ij}$ for a pair of particles, which, for each particle, must be summed to the contributions of the other interactions, to obtain the resulting force. 

The following function will, for one pair of particles, compute the above force contribution and add it to the array of forces in the positions associated to each particle:
"

# ╔═╡ f2178013-1834-4f7d-ab11-5409799b5847
function flj(x,y,i,j,d2,ε,σ,f)
	r = y - x
	@fastpow dudr = -12*ε*(σ^12/d2^7 - σ^6/d2^4)*r
	f[i] = f[i] + dudr
	f[j] = f[j] - dudr
	return f
end

# ╔═╡ 64b75ebb-5c70-4464-b6d1-c4ba15df7792
md"
As mentioned above, `r = y - x` is the vector connecting the two particles. One important detail is involved here, which is that we are writting the functions with the coordinates, and the forces, are given by arrays of 3D vectors. Typically, in Julia, these arrays will be vectors of static-vectors, or more specifically a `Vector{SVector{3,Float64}}` container. 

Thus, in the code above `dudr` is a 3D vector, because `r` is a 3D vector. And the lines that update the forces associated to each particle are updating the vector of forces in a one-line operation (`f[i] = f[i] + dudr`, for example).

The computation of the forces require the relative position of the particles, thus here we explicitly use `x` and `y` from the standard interface, and we need the indexes of the particles to update the force vector accordingly. 

We also need to provide the Lennard-Jones parameters and the force vector to be updated. The output variable here is the vectors of forces, thus it will be mutated.
"

# ╔═╡ 58d0db1e-e861-495a-95b8-5d4bd2b590cd
md"
The computation of the forces, including the construction of the cell lists, can be done with the following function:
" 

# ╔═╡ 483fcbc9-a493-4e2e-a6cb-be22888fb361
function computef(x,box)
    cl = CellList(x,box)
    f = zero(x) # vector similar to x but with zeros
    ε, σ = 0.0442, 3.28 # Neon
    map_pairwise!((x,y,i,j,d2,f) -> flj(x,y,i,j,d2,ε,σ,f), f, box, cl)
    return f
end

# ╔═╡ 7053c861-a2eb-4a24-b863-46eb2abd99ab
md"
That is, given the coordinates and the box, we map the `flj` function through the pairs of particles that are within the cutoff, updating the `f` array of forces. 

Note that, since `f` is mutale, we opt to use the `map_pairwise!` function, with the `!`, following the convention that, in Julia, functions ending with the `!` mutate one or more of its arguments. There is no need to reassign `f` from the output of `map_pairwise!`, since `f` is passed by sharing.

To generate a system with atomic dimensions, we use a convenient test-constructor available within `CellListMap`:
"

# ╔═╡ df9d4826-1b69-4671-8ad2-3de4ff0cca38
x, box = CellListMap.xatomic(10^6)

# ╔═╡ 904832a7-bf6f-4295-a6f2-e9eed0735762
md"
The above command returns random positions for `10^6` particles, and a `Box` with dimensions such that the atomic density of the system will be that of water. 

The coordinates are returned as a vector of 3D static vectors, thus, for example, the first element of `x` is the vector of coordinates of the first atom:
"

# ╔═╡ 4985bf71-739e-45e6-83bd-df7f3bb39d6e
x[1]

# ╔═╡ f07051c5-55c8-4f46-9b8d-0c179af55b1a
md"
Now, we can call `computef` with this input, to obtain the vector of forces:
"

# ╔═╡ adf92db5-a120-4e1a-88af-d4978f82cfdc
f = computef(x,box)

# ╔═╡ 87875721-cb86-4d16-aaa1-4957fa16cf73
md"
Note, again, that the output is a vector of 3D static vectors. Thus, each element of `f` is the vector of forces acting on each atom. For example, the forces acting on the first atom are:
" 

# ╔═╡ 23e4b8d6-fd6d-4f02-a46a-46f675ec23f5
f[1]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CellListMap = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
FastPow = "c0e83750-1142-43a8-81cf-6c956b72b4d1"

[compat]
CellListMap = "~0.7.13"
FastPow = "~0.1.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CellListMap]]
deps = ["DocStringExtensions", "LinearAlgebra", "Parameters", "ProgressMeter", "Random", "Setfield", "StaticArrays"]
git-tree-sha1 = "d87330f6e028ea02206bcd7491001caae23822ae"
uuid = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
version = "0.7.13"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.FastPow]]
git-tree-sha1 = "7d961335144dad74de0e1b3a9b60e4d114a78dc2"
uuid = "c0e83750-1142-43a8-81cf-6c956b72b4d1"
version = "0.1.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╟─c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
# ╟─7dbfadf4-2cb8-4145-948f-df1e7ea06b14
# ╟─b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
# ╠═562abb8e-37ed-4fac-8293-c47e7d1d6203
# ╟─dc0e8a10-30ba-44e5-be65-61b0402ed2ed
# ╠═56744cec-3cb1-45a1-b681-c99f0ae2eaae
# ╟─8d685c45-f9d7-42fd-a6e5-28361d4a7de3
# ╠═f2178013-1834-4f7d-ab11-5409799b5847
# ╟─64b75ebb-5c70-4464-b6d1-c4ba15df7792
# ╟─58d0db1e-e861-495a-95b8-5d4bd2b590cd
# ╠═483fcbc9-a493-4e2e-a6cb-be22888fb361
# ╟─7053c861-a2eb-4a24-b863-46eb2abd99ab
# ╠═df9d4826-1b69-4671-8ad2-3de4ff0cca38
# ╟─904832a7-bf6f-4295-a6f2-e9eed0735762
# ╠═4985bf71-739e-45e6-83bd-df7f3bb39d6e
# ╟─f07051c5-55c8-4f46-9b8d-0c179af55b1a
# ╠═adf92db5-a120-4e1a-88af-d4978f82cfdc
# ╟─87875721-cb86-4d16-aaa1-4957fa16cf73
# ╠═23e4b8d6-fd6d-4f02-a46a-46f675ec23f5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
