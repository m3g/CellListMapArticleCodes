### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 562abb8e-37ed-4fac-8293-c47e7d1d6203
using CellListMap

# ╔═╡ c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
md"
### Code Blocks of the CellListMap.jl paper explained
"

# ╔═╡ 7dbfadf4-2cb8-4145-948f-df1e7ea06b14
md"
# Code Block 5
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ 6ec7d8b3-d8d2-4b33-af9e-d9e01c58988a
md"
This code block shows the syntax necessary for calling the `map_pairwise` function, in minimal examples. Wee need the box and the cell list first, and we generate them here similarly to what was done in Code Block 4. First, we generate the coordinates:
"

# ╔═╡ 2744f768-1ecb-4b86-85c5-c2a54c573abf
x = rand(3,10^5)

# ╔═╡ b7089186-a3ca-4339-a837-5cc7897e7855
md"
And a `box`, orthorhombic in this case for simplicity, associated with a cutoff of `0.05`:
"

# ╔═╡ 778f8336-e6cc-4c4c-beab-de2f29c52e89
box = Box([1,1,1],0.05)

# ╔═╡ aa0a7c25-8109-464e-943d-4131824edf5f
md"
And finally the cell lists:
"

# ╔═╡ a1171f7c-f044-45da-8f35-dc576d9ed901
cl = CellList(x,box)

# ╔═╡ 69250c96-c9b8-4076-9d77-f720337d1538
md"
In the first example of the block we compute the sum of the inverse of the distance of all pairs of particles within the cutoff.

The result is an immutable value, to which the label `u` will be assigned. The initial value of `u` is `0.`, and the function to be added is provided as an anonymous function that adapts the function call to the interface required:
"

# ╔═╡ a0495ede-8aaa-4f2e-9e15-f133142bbe01
u = map_pairwise(
	(x,y,i,j,d2,u) -> u += 1 / sqrt(d2),
	0., box, cl
)

# ╔═╡ 87bafd9b-6b97-448a-b669-407e1b8aad83
md"
The inner interface of the function to be mapped presuposes a call with six arguments, `(x,y,i,j,d2,u)`, but here we only need the squared distance `d2` and the output to be updated, `u`. The function to be mapped, therefore, just ignores the additional parameters and sums to `u` the inverse of the distance between the particles.

Note that the `cutoff` is contained in the `box` structure, and that function will only be evaluated for the pairs of particles with distances smaller than the cutoff.  
" 

# ╔═╡ 5d25e653-5d5c-444b-a04a-045d3e7736ca
md"
In the second example, we provide an additional parameter, the \"masses\" of the particles, which are not part of the inner interface of the mapped function:
"

# ╔═╡ b17a2f6d-0f17-435a-98af-241dba0cad05
const mass = rand(10^5)

# ╔═╡ 5bcf60eb-37e4-4437-b5da-a701d400c5be
u2 = map_pairwise(
	(x,y,i,j,d2,u) -> u += mass[i]*mass[j] / sqrt(d2),
	0., box, cl
)

# ╔═╡ a1ba128f-a2eb-4070-b508-8449ae8bb61a
md"
The important difference relative to the previous example is that now the annonymous function \"closes-over\" the masses. That is, it uses the masses from the outer scope of the `map_pairwise` function, to compute the mapped function. This allows for general property implementations. 

Some care has to be taken to guarantee that the closed-over variables are type-stable, by keeping everying in a function, for example, or declaring the parameters constant, as we did here with the masses. This is important performance-wise, in Julia.
" 

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CellListMap = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"

[compat]
CellListMap = "~0.7.13"
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
# ╟─6ec7d8b3-d8d2-4b33-af9e-d9e01c58988a
# ╠═2744f768-1ecb-4b86-85c5-c2a54c573abf
# ╟─b7089186-a3ca-4339-a837-5cc7897e7855
# ╠═778f8336-e6cc-4c4c-beab-de2f29c52e89
# ╟─aa0a7c25-8109-464e-943d-4131824edf5f
# ╠═a1171f7c-f044-45da-8f35-dc576d9ed901
# ╟─69250c96-c9b8-4076-9d77-f720337d1538
# ╠═a0495ede-8aaa-4f2e-9e15-f133142bbe01
# ╟─87bafd9b-6b97-448a-b669-407e1b8aad83
# ╟─5d25e653-5d5c-444b-a04a-045d3e7736ca
# ╠═b17a2f6d-0f17-435a-98af-241dba0cad05
# ╠═5bcf60eb-37e4-4437-b5da-a701d400c5be
# ╟─a1ba128f-a2eb-4070-b508-8449ae8bb61a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
