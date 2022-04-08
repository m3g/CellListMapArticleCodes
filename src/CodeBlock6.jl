### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 562abb8e-37ed-4fac-8293-c47e7d1d6203
using CellListMap

# ╔═╡ 56744cec-3cb1-45a1-b681-c99f0ae2eaae
using FastPow

# ╔═╡ f855ad3c-c009-4077-a3bb-a7c18f790eab
using BenchmarkTools

# ╔═╡ c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
md"
### Code Blocks of the CellListMap.jl paper explained
"

# ╔═╡ 7dbfadf4-2cb8-4145-948f-df1e7ea06b14
md"
# Code Block 6
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ dc0e8a10-30ba-44e5-be65-61b0402ed2ed
md"
This example involves the computation of Lennard-Jones potential, in which exponential terms with high powers have to be computed. We use the `FastPow` package to automatically decompose these into smaller powers, accelerating the computation:
"

# ╔═╡ 8d685c45-f9d7-42fd-a6e5-28361d4a7de3
md"
The Lennard-Jones potential for a pair of particles has the form

$u = 4\varepsilon \left[\left(\frac{\sigma}{d}\right)^{12} 
                           - \left(\frac{\sigma}{d}\right)^6\right]$

where $d$ is the distance between the particles, $\epsilon$ is the depth of the potential well, and $\sigma$ is the combined radii of the particles. 

Since the interface provided by `CellListMap` gives us the squared distance, we rewrite this expression as:

$u = 4\varepsilon \left[\left(\frac{\sigma}{d^2}\right)^{6} 
                           - \left(\frac{\sigma}{d^2}\right)^3\right]$

which translates almost identically to the corresponding function that will be mapped:
"

# ╔═╡ f2178013-1834-4f7d-ab11-5409799b5847
ulj(d2,ε,σ,u) = @fastpow u += 4ε*((σ^2/d2)^6 - (σ^2/d2)^3)

# ╔═╡ 4cfdc28c-66c8-4924-9069-b955ccc22c98
md"
In the example of this block, we construct a system which has a physical meaning: it is a fluid of Neon atoms, with an atomic density similar to that of water. With that, the computational cost of the computation is of the order of what will be observed for most condensed-phase systems. 

We choose to evaluate the potential for 3 million atoms, implying then a cubic box of `31.034` nm. The cutoff chosen for the computation here is `1.2`nm:
"

# ╔═╡ 2ee76786-d6a2-404e-86aa-233c9dc2e73e
side = 31.034

# ╔═╡ fdeae80b-7b8c-4172-9013-b629d3fd71c6
cutoff = 1.2

# ╔═╡ 366f7f16-e79e-450c-a8d0-2f96325db22e
md"
With those dimensions in hand, we generate 3 million 3D coordinates distributed randomly in the cube of side `side`:
" 

# ╔═╡ 297ffbe6-ad64-4cd8-b15a-907534934300
x = side*rand(3,3_000_000)

# ╔═╡ d715aeeb-ea5a-4ebd-b375-dc2ec365fd2e
md"
Next, we set the `Box` and compute the cell lists:
" 

# ╔═╡ 34d7a971-277d-4a1e-8f00-cdd273258d56
box = Box([side,side,side],cutoff)

# ╔═╡ 31ea8abd-9ad1-45a4-8ff5-7446dcaab9eb
cl = CellList(x,box)

# ╔═╡ 06c63e49-803b-40f5-a093-2b5630540b8f
md"
We define the parameters of the Lennard-Jones potential for Neon:
"

# ╔═╡ f3b9e62d-5c60-4964-b0aa-25ebe3e27c96
const ε, σ = 0.0442, 3.28

# ╔═╡ 03b9b4de-7191-4759-bacf-8c69022b6d8d
md"
With the system set up, and the mapped function properly defined, we can now compute the Lennard-Jones potential for the complete set of 3 million particles, with periodic boundary conditions, using:
"

# ╔═╡ 17029fc8-3219-4614-bf90-2967a3101a53
u = map_pairwise((x,y,i,j,d2,u) -> ulj(d2,ε,σ,u), 0., box, cl)

# ╔═╡ 14a10777-fed7-41ef-b109-c104e9a98235
md"
We can benchmark properly this computation with the `BenchmarkTools` package. Here we have executed the notebook in a laptop with 6 threads:
" 

# ╔═╡ 40248541-2779-4024-9952-3e530644703c
Threads.nthreads()

# ╔═╡ c3149400-bd1e-4501-be02-ef6af6324235
@benchmark map_pairwise((x,y,i,j,d2,u) -> ulj(d2,$ε,$σ,u), 0., $box, $cl)

# ╔═╡ 1676a0c0-db3d-444b-8205-aa73081d9379
md"
Thus, the computation required something of the order of 2 seconds to finish. 
" 

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CellListMap = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
FastPow = "c0e83750-1142-43a8-81cf-6c956b72b4d1"

[compat]
BenchmarkTools = "~1.3.1"
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

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

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

[[deps.FastPow]]
git-tree-sha1 = "7d961335144dad74de0e1b3a9b60e4d114a78dc2"
uuid = "c0e83750-1142-43a8-81cf-6c956b72b4d1"
version = "0.1.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

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

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

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
# ╟─4cfdc28c-66c8-4924-9069-b955ccc22c98
# ╠═2ee76786-d6a2-404e-86aa-233c9dc2e73e
# ╠═fdeae80b-7b8c-4172-9013-b629d3fd71c6
# ╟─366f7f16-e79e-450c-a8d0-2f96325db22e
# ╠═297ffbe6-ad64-4cd8-b15a-907534934300
# ╟─d715aeeb-ea5a-4ebd-b375-dc2ec365fd2e
# ╠═34d7a971-277d-4a1e-8f00-cdd273258d56
# ╠═31ea8abd-9ad1-45a4-8ff5-7446dcaab9eb
# ╟─06c63e49-803b-40f5-a093-2b5630540b8f
# ╠═f3b9e62d-5c60-4964-b0aa-25ebe3e27c96
# ╟─03b9b4de-7191-4759-bacf-8c69022b6d8d
# ╠═17029fc8-3219-4614-bf90-2967a3101a53
# ╟─14a10777-fed7-41ef-b109-c104e9a98235
# ╠═f855ad3c-c009-4077-a3bb-a7c18f790eab
# ╠═40248541-2779-4024-9952-3e530644703c
# ╠═c3149400-bd1e-4501-be02-ef6af6324235
# ╟─1676a0c0-db3d-444b-8205-aa73081d9379
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
