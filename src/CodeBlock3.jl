### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 562abb8e-37ed-4fac-8293-c47e7d1d6203
using CellListMap

# ╔═╡ a5c65b9c-4cc3-4c68-be5f-e2ab89b18581
using Unitful

# ╔═╡ c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
md"
### Code Blocks of the CellListMap.jl paper explained
"

# ╔═╡ 7dbfadf4-2cb8-4145-948f-df1e7ea06b14
md"
# Code Block 3
**Initialization of the system Box, with orthorhombic or triclinic periodic boundary conditions. The system’s geometry is defined by the type of unit cell matrix, and an orthorhombic cell is assumed if a vector of box sides is supplied.**
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ c0027b3b-d932-410d-b64b-ea63689be612
md"
Here we show the types of periodic boundary conditions that `CellListMap` supports, and how the system box can be initialized.

The relevant function here is the `Box` constructor, which will create an data structure of type `CellListMap.Box`, specific for unit cell type, type of data, and dimension of the system. 

The first example is:
"

# ╔═╡ a1171f7c-f044-45da-8f35-dc576d9ed901
box1 = Box([10,20,15],1.2)

# ╔═╡ 1b09e4cf-a776-46a0-bad5-14fa5f21bf14
md"
The input of the `Box` constructor contained two parameters: 

- The box sides, given here by the three-dimensional vector `[10,20,15]`
- The cutoff, `1.2`

Thus, it is implicit here that the box is orthorhombic (because only sides are given), and that the system is three-dimensional (otherwise a 2D vector of sides would be given). It is implicit given that `Float64` numbers are used, from the type of the cutoff variable, `typeof(1.2) == Float64`.

Other properties of the `Box` are provided as output: the number of cells of the grid in each dimension, the computing cell sizes, the `lcell` parameter ([which can be tuned for efficiency](https://m3g.github.io/CellListMap.jl/stable/performance/#Optimizing-the-cell-grid)), and the final total number of cells.

A general triclinic unit cell can be used by providing instead of a vector of box sides, a complete unit cell matrix, where the columns of the matrix are the unit cell vectors. This is illustrated in the second example of the code block:
"

# ╔═╡ b1643bac-3de5-4afb-a370-990d0ec21c77
box2 = Box([ 10  0  0
              0 10  0
              0 20 10 ], 1.2)

# ╔═╡ 981183bb-abfc-45ef-8831-e54dec0d078c
md"
Note that the type of `Box` is now `Box{TriclinicCell}`. This is currently important because some specializations for `OrthorhombicCells` can be done to improve performance. 

To provide a constrasting example, here we construct a box with different dimensions and types of variables:
" 

# ╔═╡ 90070e77-4e4f-4db1-b7a1-af7fcf233774
box3 = Box( [10.f0, 20.f0 ]u"nm", 1.2f0u"nm")

# ╔═╡ 4ae9154e-c461-4e55-b09a-41b55160a963
md"
In this example we have changed the type of input variable, which are now `Unitful` quantities represented by `Float32` reals with `nm` units. Additionaly, the example is two-dimensional. And this is reflected in the (now somewhat extensive) `Box` type signature. 
" 

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CellListMap = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
CellListMap = "~0.7.13"
Unitful = "~1.11.0"
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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

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

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╟─c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
# ╟─7dbfadf4-2cb8-4145-948f-df1e7ea06b14
# ╟─b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
# ╠═562abb8e-37ed-4fac-8293-c47e7d1d6203
# ╟─c0027b3b-d932-410d-b64b-ea63689be612
# ╠═a1171f7c-f044-45da-8f35-dc576d9ed901
# ╟─1b09e4cf-a776-46a0-bad5-14fa5f21bf14
# ╠═b1643bac-3de5-4afb-a370-990d0ec21c77
# ╟─981183bb-abfc-45ef-8831-e54dec0d078c
# ╠═a5c65b9c-4cc3-4c68-be5f-e2ab89b18581
# ╠═90070e77-4e4f-4db1-b7a1-af7fcf233774
# ╟─4ae9154e-c461-4e55-b09a-41b55160a963
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
