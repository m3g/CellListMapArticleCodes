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
# Code Block 2
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ c0027b3b-d932-410d-b64b-ea63689be612
md"
In this block we just exemplify the general format that the function to be computed for the pairs of particles within the cutoff must have to conform with the interface of the package:
"

# ╔═╡ a1171f7c-f044-45da-8f35-dc576d9ed901
function f(x,y,i,j,d2,output,args...)
    # evaluate property for pair i,j and update output variable
    return output
end

# ╔═╡ 247a33ae-2796-4920-859f-869560164632
md"
The function interface provides six arguments `(x,y,i,j,d2,output)`, which are:

- `x`: The coordinates of the first particle.
- `y`: The coordinates of the second particle (minimum image relative to `x`).
- `i`: The index in the array of coordinates of the first particle.
- `j`: The index in the array of coordinates of the second particle.
- `d2`: The squared distance between the particles (i. e. $||y-x||^2$)
- `output`: the value of the output variable.

Additional arguments can be provided if the function requires them, for examples particle masses, etc. 

Let us generate some coordinates, a system box, and the cell lists, to be able to provide examples:
"

# ╔═╡ d10a0c9a-b669-46c5-af29-917984b65b31
x = rand(3,10^5)

# ╔═╡ bd95214b-a0c8-4b11-baf2-6fe59a54a5c4
box = Box([1,1,1],0.05)

# ╔═╡ 1dd597ee-3e67-440e-9b1e-8f4d66bd9698
cl = CellList(x,box)

# ╔═╡ 9bed3fb2-6aa2-4980-8b08-c1de9ea42c20
md"
One of the simplest computations possible is the sum of the distances. The function could be defined like this:
"

# ╔═╡ 7037d17d-7e7e-4119-bb5c-e74898dc1181
function sum_distances(d2,output)
    output += sqrt(d2)
	return output
end

# ╔═╡ f96176a5-d9d2-44be-8a21-2c4c1f045466
md"
But this function does not conform with the interface. Thus, when calling `map_pairwise` function, we define an anonymous function that adjusts the function call:
"

# ╔═╡ 79d91bc5-d862-4c8a-a691-d6c89ee8481e
map_pairwise( (x,y,i,j,d2,output) -> sum_distances(d2,output), 0., box, cl )

# ╔═╡ bbd85682-65ee-4d0c-b12c-36b1d9fc1fda
md"
Now, let us suppose that each particle has a mass, and we want to compute the product of the masses divided by the distance. Here are the mass arrays:
"

# ╔═╡ 1ff7a23d-560e-47d8-9f64-c8173681d167
mass = rand(10^5)

# ╔═╡ 52d9b1df-daef-4aef-8125-ff8cd54d253c
md"
The function to be computed is, now:
"

# ╔═╡ 59a19df0-607e-4c1f-b9ba-590f494196c5
function u(i,j,d2,output,mass)
	output += mass[i]*mass[j]/sqrt(d2)
	return output
end

# ╔═╡ 9c0a0719-1d0d-4d1a-9c67-6458dce71105
md"
And finaly, to conform to the interface of `map_pairwise`, an anonymous function is used in the call, which closes over the masses:
"

# ╔═╡ e390cdb4-94e2-4435-8014-44973d789648
map_pairwise( (x,y,i,j,d2,output) -> u(i,j,d2,output,mass), 0., box, cl )

# ╔═╡ 155dc0e8-afa0-4b57-b5bb-bd5d6aa116ce
md"
Thus, while the input interface of the function is strict, any combination of the parameters required and additional parameters can be used.
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
# ╟─c0027b3b-d932-410d-b64b-ea63689be612
# ╠═a1171f7c-f044-45da-8f35-dc576d9ed901
# ╟─247a33ae-2796-4920-859f-869560164632
# ╠═d10a0c9a-b669-46c5-af29-917984b65b31
# ╠═bd95214b-a0c8-4b11-baf2-6fe59a54a5c4
# ╠═1dd597ee-3e67-440e-9b1e-8f4d66bd9698
# ╟─9bed3fb2-6aa2-4980-8b08-c1de9ea42c20
# ╠═7037d17d-7e7e-4119-bb5c-e74898dc1181
# ╟─f96176a5-d9d2-44be-8a21-2c4c1f045466
# ╠═79d91bc5-d862-4c8a-a691-d6c89ee8481e
# ╟─bbd85682-65ee-4d0c-b12c-36b1d9fc1fda
# ╠═1ff7a23d-560e-47d8-9f64-c8173681d167
# ╟─52d9b1df-daef-4aef-8125-ff8cd54d253c
# ╠═59a19df0-607e-4c1f-b9ba-590f494196c5
# ╟─9c0a0719-1d0d-4d1a-9c67-6458dce71105
# ╠═e390cdb4-94e2-4435-8014-44973d789648
# ╟─155dc0e8-afa0-4b57-b5bb-bd5d6aa116ce
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
