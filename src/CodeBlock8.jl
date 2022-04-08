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
# Code Block 8
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ 31fa7f72-0111-413c-97c6-dd3d1926527f
md"
In this block we implement the computation of the computation of a list of k nearest neighbors between two sets of particles. The purpose of example is to describe how to implement custom functions which require the implementation of custom parallel reduction functions.

The list of neighbors will be a list of named tuples, for the form `(i=...,j=...,dsq=...)` where `i` and `j` are the indexes of the two atoms involved in the contact, and `dsq` is the squared distance between them.
"

# ╔═╡ 9b4c50bf-f8b8-4d07-bc4e-792072a74891
md"
### Updating a list of neighbors

To start, we need to define a function that, given a current `list` of neighbors, updates the list of the distance between two particles is found to be shorter for one pair. 

For example, let us suppose that we have the following list:
"

# ╔═╡ f313f454-9352-4dbf-89c8-3d588cae332d
list0 = [ (i=1,j=2,dsq=5.), (i=1,j=3,dsq=7.) ]

# ╔═╡ 2ca6d0b3-2348-4bfb-84d2-30f738f86e98
md"
If we find a new pair for which the distance is, for example, greater than `5.` but smaller than `7.`, the second pair of the list must be replaced. This new pair will be the pair `(i=1,j=4,dsq=6.)`.

The function below implements this updating of the list of pairs. It receives the `list`, the indexes of the particles `i`, and `j`, and the squared distance between them, `d2`:
"

# ╔═╡ 5fff0a3c-248e-41bd-9709-d359d0ddf903
function replace_pair!(list,i,j,d2)
    pair = (i=i,j=j,dsq=d2)
    ipos = searchsortedfirst(list, pair, by = p -> p.dsq)
    if ipos <= length(list)
        list[ipos+1:end] = list[ipos:end-1]
        list[ipos] = pair
    end
    return list
end

# ╔═╡ 779c12fe-8b73-4c29-ad2e-d1775d0d0ae2
md"
The function first creates the named named tuple from the input indexes and distances. Then it finds which is the index of the first pair which has a `dsq` field greater than the distance `d2` of the input. This is the output of the `searchsortedfirst` function call. 

Next, if the position is within the list (`ipos <= length(list)`), the list is modified by the insertion of the new pair at the specified position, and shifting all other pairs (with distances greater than the current one) one position further in the list, eliminating the last pair (the order of these operations is the opposite).

Finally, the list is returned, updated according to the new pair found.

Thus, given `list0` above and the new pair, we do:
"

# ╔═╡ 78efc6d8-fd52-44cb-8976-797cf3ed87f8
replace_pair!(list0, 1, 4, 6.)

# ╔═╡ 82a04ec5-9ed2-4847-8a24-6055cdf880a2
md"
obtaining a new list, in which the last pair was replaced by the new one. 

The function `replace_pair!` will be called for every pair of atoms found to be within the desired cutoff, promoting the list update if the distance between the particles is shorter than the greatest distance already available in the list.
"

# ╔═╡ 7c055c92-0163-41a5-80d4-525b76fa1b53
md"
### Reducing lists computed in parallel

To compute the list of neighbors in parallel we need to implement a custom reduction function. Thus, let us supose that we have computed the lists in parallel for independent subsets of the particles, in two threads. Two independent lists would have been obtained, for example:
"

# ╔═╡ 739180a6-b603-4968-9464-504aafb7534c
list_threaded = [ 
	[ (i=1,j=2,dsq=1.), (i=2,j=3,dsq=10.)], # list obtained in the first thread
	[ (i=1,j=3,dsq=5.), (i=2,j=4,dsq=11.)]  # list obtained in the second thread
]

# ╔═╡ 1ffb4217-a9ac-4385-b237-c2122681e453
md"
The combination of the two lists above should preserve the first pair of each list, because these are the two shorter distances found. 

The following function will produce the correct result, by updating a third input list with the data of the two lists above. 
"

# ╔═╡ 4faa0833-f61f-4c26-b0a2-0e630fabdbc5
function reduce_list(list,list_threaded)
    for lst in list_threaded, pair in lst
        replace_pair!(list,pair...)
    end
    return list
end

# ╔═╡ 4b1d2cde-58a6-40fb-bf99-0b1e9cf37003
md"
The function will update the input list_reduced, which is initialized with +Inf values for the distances:
"

# ╔═╡ b6585abe-e416-4ca5-911c-8abe4c1475f7
list_reduced = [(i=0,j=0,dsq=+Inf),(i=0,j=0,dsq=+Inf)]

# ╔═╡ 9bf00583-3801-4181-ab74-efcd91a927ae
md"
Applying the function above to the threaded list of lists, we get:
"

# ╔═╡ 3d648dc9-78ac-4504-9de1-08c331282892
reduce_list(list_reduced, list_threaded)

# ╔═╡ fa07cbd1-7e05-48a3-93b3-3ba7b0863d13
md"
which, as expected, retains the two pairs with the smaller distances. 
"

# ╔═╡ d957bc3e-ef5b-4c9b-9538-af2c1ce5a961
md"
### Application

Now, let us apply this complete solution to two sets of particles. We will create two sets of particles with 100 and 1000 points, as `3 x N` matrices, where the columns represent the coordinates of each particle:

"

# ╔═╡ c199024a-e6c5-48a1-bf64-512985678ecc
x = rand(3,100)

# ╔═╡ d8bd5a83-17fe-469a-aa6f-12b3a997386e
y = rand(3,1000)

# ╔═╡ 2a46bc74-1163-44a7-b0fc-69ad4b417b0c
md"
Our goal will be to obtain the 5 pairs of particles with shorter distances. Thus, we define an initial list of pairs with 5 elements, and assign `dsq=+Inf` to each:
"

# ╔═╡ 69ae6c59-57b5-426d-83db-16a0ae102cec
list = [ (i=0,j=0,dsq=+Inf) for _ in 1:5 ]

# ╔═╡ 17bf48e9-1cdb-4e16-8148-9c67ac13a881
md"
In this toy example, we will use a simple cubic box of side `1.0`, and a cutoff of `0.1`:
"

# ╔═╡ 1c5b099c-c61d-4fe1-a11e-99a1bf534f48
box = Box([1,1,1],0.1)

# ╔═╡ e4629f37-2863-4fc4-bdad-35e4b43e245a
md"
We construct the cell lists, and use this example to show how to call the constructor when the calculation involves two sets of particles:
"

# ╔═╡ 32e2d69f-df59-4548-8534-244e87479b06
cl = CellList(x,y,box)

# ╔═╡ bd2ec98e-fdb8-4d91-84d6-6cfb1cb66908
md"
Finally, we map the `replace_pair!` function to the pairs of particles within the cutoff. We need to inform that a custom reduction function was defined. The result will contain the 5 pairs of particles with the smaller distances between the two sets:
"

# ╔═╡ 19cc9c02-4b7c-4970-b82e-544bc1840eae
map_pairwise!(
    (x,y,i,j,d2,list) -> replace_pair!(list,i,j,d2), 
    list, box, cl,
    reduce=reduce_list
)

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
# ╟─31fa7f72-0111-413c-97c6-dd3d1926527f
# ╟─9b4c50bf-f8b8-4d07-bc4e-792072a74891
# ╠═f313f454-9352-4dbf-89c8-3d588cae332d
# ╟─2ca6d0b3-2348-4bfb-84d2-30f738f86e98
# ╠═5fff0a3c-248e-41bd-9709-d359d0ddf903
# ╟─779c12fe-8b73-4c29-ad2e-d1775d0d0ae2
# ╠═78efc6d8-fd52-44cb-8976-797cf3ed87f8
# ╟─82a04ec5-9ed2-4847-8a24-6055cdf880a2
# ╟─7c055c92-0163-41a5-80d4-525b76fa1b53
# ╠═739180a6-b603-4968-9464-504aafb7534c
# ╟─1ffb4217-a9ac-4385-b237-c2122681e453
# ╠═4faa0833-f61f-4c26-b0a2-0e630fabdbc5
# ╟─4b1d2cde-58a6-40fb-bf99-0b1e9cf37003
# ╠═b6585abe-e416-4ca5-911c-8abe4c1475f7
# ╟─9bf00583-3801-4181-ab74-efcd91a927ae
# ╠═3d648dc9-78ac-4504-9de1-08c331282892
# ╟─fa07cbd1-7e05-48a3-93b3-3ba7b0863d13
# ╟─d957bc3e-ef5b-4c9b-9538-af2c1ce5a961
# ╠═c199024a-e6c5-48a1-bf64-512985678ecc
# ╠═d8bd5a83-17fe-469a-aa6f-12b3a997386e
# ╟─2a46bc74-1163-44a7-b0fc-69ad4b417b0c
# ╠═69ae6c59-57b5-426d-83db-16a0ae102cec
# ╟─17bf48e9-1cdb-4e16-8148-9c67ac13a881
# ╠═1c5b099c-c61d-4fe1-a11e-99a1bf534f48
# ╟─e4629f37-2863-4fc4-bdad-35e4b43e245a
# ╠═32e2d69f-df59-4548-8534-244e87479b06
# ╟─bd2ec98e-fdb8-4d91-84d6-6cfb1cb66908
# ╠═19cc9c02-4b7c-4970-b82e-544bc1840eae
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
