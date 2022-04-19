### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 562abb8e-37ed-4fac-8293-c47e7d1d6203
using CellListMap

# ╔═╡ 23553122-2a28-452b-adfc-2fda2fb1a480
using Unitful

# ╔═╡ 0fac3c08-a0a5-42f9-846e-ee6e3fef677a
using Measurements

# ╔═╡ aa105270-a04a-4952-90fd-e94777d9ec41
using ForwardDiff

# ╔═╡ a509b035-25fb-4198-b314-03c21489a7d0
using BenchmarkTools

# ╔═╡ c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
md"
### Code Blocks of the CellListMap.jl paper explained
"

# ╔═╡ 7dbfadf4-2cb8-4145-948f-df1e7ea06b14
md"
# Code Block 9
**Units, uncertainties, and automatic differentiation propagating through pairwise computations with CellListMap.jl.**
"

# ╔═╡ b68cfd92-d1b3-4f33-b9bf-8f72b101a5b9
md"
This code block illustrates the use of units, measurements, and automatic differentation, in combination with `CellListMap`.
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ 1de0a778-1222-4791-8012-7cd51b3372cd
md"
### Units

A popular package to deal with units in Julia is `Unitful`. Here, we reproduce the simplest complete running code of Code Block 2, but using units for the variables.

First, we load the `Unitful` package:
 "

# ╔═╡ b17481f8-b5d0-41fb-8223-9c8704ac94b0
md"
and, now, we reproduce the computation of the sum of the distances of 1000 particles, but where the units of the coordinates are nanometers. The particle coordinates are represented as a `3 x 1000` matrix, where the columns are associated to each particle. The complete code is, in this notebook, wrapped in a `let` block because `Pluto` (this notebook) requires on single instruction in each block:
"

# ╔═╡ 21b40e29-f829-4675-8bba-9ff00e664727
let
	x = rand(3,1000)u"nm";
	box = Box([1.,1.,1.]u"nm",0.05u"nm")
	cl = CellList(x,box)
	map_pairwise((x,y,i,j,d2,out) -> out += sqrt(d2), 0.0u"nm", box, cl)
end

# ╔═╡ 4d6f9c0f-84d9-41a3-9c63-de98d67899bd
md"
It is important to be careful in setting correctly the units for all quantities involved, meaning the coordinates, the box sides and cutoff, and the initial value of the output variable `0.0u\"nm\"` in this case.
"

# ╔═╡ 398359f8-dde3-47df-8716-d722063f12c3
md"
### Propagation of Uncertainties

The `Measurements` package provides the support for propagating uncertainties:

"

# ╔═╡ 66e8272c-d814-4f8d-999b-e348450ef87d
md"
With this package, one can propagate measurements in operations with quantities with uncertainties, for example:
"

# ╔═╡ 3a984a74-61a9-49d1-bb19-c935c2de6a49
(1.0 ± 0.1) * (1.0 ± 0.2)

# ╔═╡ 296616d8-061b-4ad1-9333-be34c25fe8b3
md"
We can do the same with the complete calculation mapped through `CellListMap`. Here, we illustrate the calculation by providing coordinates as a vector of 3D vectors, with random coordinates containing uncertainties. Other quantities have to be defined with uncertainties as well, and the result will have the uncertainties propagated to the final result:
"

# ╔═╡ e586aedc-06b1-492c-8aa3-ca0c83e9b1b3
let
	x = [ [rand() ± 0.1 for _ in 1:3 ] for _ in 1:1000 ]
	box = Box([1 ± 0; 1 ± 0; 1 ± 0], 0.05 ± 0.0);
	cl = CellList(x,box);
	map_pairwise((x,y,i,j,d2,out) -> out += sqrt(d2), 0. ± 0., box, cl)
end

# ╔═╡ ed4eed22-0d7e-41e4-aa19-5862f9b6214c
md"
### Automatic differentiation

Using automatic differentiation can be useful if the computation will be used in the context of an optimization problem.

To use automatic differentiation, we write the function to be computed with general types. The type `T` of the function below is the type of number representing the coordinates. Normally these would be `Float32` of `Float64` numbers, but by writting the function with generic types one allows the propagation of the dual numbers required for the computation of automatic derivatives:
"

# ╔═╡ 15b78386-01a7-412d-9080-f67d8a30d666
function sum_d(x::Matrix{T},sides,cutoff) where T
	box = Box(T.(sides),T(cutoff))
	cl = CellList(x,box)
	return map_pairwise(
		(x,y,i,j,d2,out) -> out += sqrt(d2),
		zero(T), box, cl
	)
end

# ╔═╡ 668e2a8b-2dba-4ec0-ac1a-c5131b80258f
md"
Let us define the side and the cutoff, as usual, noting that they will be converted to the appropriate type `T` inside the function above:
"

# ╔═╡ 5ef17bd3-372c-49b7-83e4-6ab4247526bb
sides = [1,1,1]

# ╔═╡ f2a0cd43-e76d-4ed4-bcae-5911a5566b7c
cutoff = 0.05

# ╔═╡ 7ef7c828-8015-4f64-87e0-1b5b6d1cb3e3
md"
And the coordinates will be represented in matricial form:
"

# ╔═╡ 29e75f6a-f1f7-46b0-8244-82a8cc5c7bef
coordinates = rand(3,1000)

# ╔═╡ e9387851-f568-4490-b60a-725aa4dcdddf
md"
Calling the `sum_d` function gives us the result:
"

# ╔═╡ 769d963e-2bd0-4f02-ac48-7bea0f1a9337
sum_d(coordinates,sides,cutoff)

# ╔═╡ e96f0d37-782f-4a91-ac04-f9b8d5272433
md"
and, finally, we can automatically differentiate this function, here, for example, with the `ForwardDiff` package:
"

# ╔═╡ fe0e4fd7-6f9c-4457-8dc2-1ede4b7f3835
ForwardDiff.gradient(x -> sum_d(x,sides,cutoff), coordinates)

# ╔═╡ 5327ace0-1b45-49a0-b056-bb6d9a589fdf
md"
In this simple example, it is easy to confirm that this result is correct, by computing the gradient by hand, also using `CellListMap`:
"

# ╔═╡ 61f27118-7aff-46ea-9604-7596915b95fd
function update_∇f!(x,y,i,j,d2,∇f)
	r = x - y
	∇f[:,i] .= r / sqrt(d2)
	∇f[:,j] .= -r / sqrt(d2)
	return ∇f
end

# ╔═╡ b2e15594-43f0-4998-9f2b-93191122a858
function ∇sum_d(x::Matrix{T},sides,cutoff) where T
	box = Box(T.(sides),T(cutoff))
	cl = CellList(x,box)
	∇f = zeros(T,size(x))
	map_pairwise!(update_∇f!, ∇f, box, cl)
	return ∇f
end

# ╔═╡ 5e06aabe-a051-4635-93c2-0f0cb7b7775c
∇sum_d(coordinates, sides, cutoff)

# ╔═╡ 64de2303-3fca-4b84-b235-7e9cb4a31da0
md"
It is expected, though, that the manual implementation is faster, in general, and this is the case here:
"

# ╔═╡ 9e53b2aa-88f5-4634-8d25-c2ad1b45fd05
@benchmark ForwardDiff.gradient(x -> sum_d(x,$sides,$cutoff), $coordinates)

# ╔═╡ 98418030-ece5-40b4-81c1-e967e79170c8
@benchmark ∇sum_d($coordinates, $sides, $cutoff)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CellListMap = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
BenchmarkTools = "~1.3.1"
CellListMap = "~0.7.13"
ForwardDiff = "~0.10.25"
Measurements = "~2.7.1"
Unitful = "~1.11.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CellListMap]]
deps = ["DocStringExtensions", "LinearAlgebra", "Parameters", "ProgressMeter", "Random", "Setfield", "StaticArrays"]
git-tree-sha1 = "d87330f6e028ea02206bcd7491001caae23822ae"
uuid = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
version = "0.7.13"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "88cd033eb781c698e75ae0b680e5cef1553f0856"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.7.1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

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

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
# ╟─7dbfadf4-2cb8-4145-948f-df1e7ea06b14
# ╟─b68cfd92-d1b3-4f33-b9bf-8f72b101a5b9
# ╟─b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
# ╠═562abb8e-37ed-4fac-8293-c47e7d1d6203
# ╟─1de0a778-1222-4791-8012-7cd51b3372cd
# ╠═23553122-2a28-452b-adfc-2fda2fb1a480
# ╟─b17481f8-b5d0-41fb-8223-9c8704ac94b0
# ╠═21b40e29-f829-4675-8bba-9ff00e664727
# ╟─4d6f9c0f-84d9-41a3-9c63-de98d67899bd
# ╟─398359f8-dde3-47df-8716-d722063f12c3
# ╠═0fac3c08-a0a5-42f9-846e-ee6e3fef677a
# ╟─66e8272c-d814-4f8d-999b-e348450ef87d
# ╠═3a984a74-61a9-49d1-bb19-c935c2de6a49
# ╟─296616d8-061b-4ad1-9333-be34c25fe8b3
# ╠═e586aedc-06b1-492c-8aa3-ca0c83e9b1b3
# ╟─ed4eed22-0d7e-41e4-aa19-5862f9b6214c
# ╠═15b78386-01a7-412d-9080-f67d8a30d666
# ╟─668e2a8b-2dba-4ec0-ac1a-c5131b80258f
# ╠═5ef17bd3-372c-49b7-83e4-6ab4247526bb
# ╠═f2a0cd43-e76d-4ed4-bcae-5911a5566b7c
# ╟─7ef7c828-8015-4f64-87e0-1b5b6d1cb3e3
# ╠═29e75f6a-f1f7-46b0-8244-82a8cc5c7bef
# ╟─e9387851-f568-4490-b60a-725aa4dcdddf
# ╠═769d963e-2bd0-4f02-ac48-7bea0f1a9337
# ╟─e96f0d37-782f-4a91-ac04-f9b8d5272433
# ╠═aa105270-a04a-4952-90fd-e94777d9ec41
# ╠═fe0e4fd7-6f9c-4457-8dc2-1ede4b7f3835
# ╟─5327ace0-1b45-49a0-b056-bb6d9a589fdf
# ╠═61f27118-7aff-46ea-9604-7596915b95fd
# ╠═b2e15594-43f0-4998-9f2b-93191122a858
# ╠═5e06aabe-a051-4635-93c2-0f0cb7b7775c
# ╟─64de2303-3fca-4b84-b235-7e9cb4a31da0
# ╠═a509b035-25fb-4198-b314-03c21489a7d0
# ╠═9e53b2aa-88f5-4634-8d25-c2ad1b45fd05
# ╠═98418030-ece5-40b4-81c1-e967e79170c8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
