### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 562abb8e-37ed-4fac-8293-c47e7d1d6203
using CellListMap

# ╔═╡ 3e39c829-f965-46a0-86e5-d08e4a101fdb
using LinearAlgebra: norm

# ╔═╡ 1ed49e14-1cca-4aae-b542-842a90e05ab9
using StaticArrays

# ╔═╡ c319e34e-b6a6-11ec-0b3e-c760ba59e0e0
md"
### Code Blocks of the CellListMap.jl paper explained
"

# ╔═╡ 7dbfadf4-2cb8-4145-948f-df1e7ea06b14
md"
# Code Block 10
"

# ╔═╡ b68cfd92-d1b3-4f33-b9bf-8f72b101a5b9
md"
In this block we implement the computation of a histogram of \"galactic\" pair velocities, which is a common computation in astrophysical simulation analsys. The computation consists in obtainin an histogram of the relative velocities of the galaxies, as a function of their distances. 

Thus, for each distance between the galaxies, we want to obtain the average relative speed between the galaxies. This involves counting the number of pairs of galaxies at each distance, summing up their relative velocities, and averaging at the end.

Besides being an example which implements a practical application of the package outside the field of molecular simulations, this codes illustrates how to use input properties that are not needed for the construction of the cell lists (the velocities), and updating mutable outputs.
"

# ╔═╡ b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
md"
Load the CellListMap package
"

# ╔═╡ 34f87c93-ef98-41be-97f1-7fda2bfb5126
md"
And here we will use the `norm` function from `LinearAlgebra`, for simplicity:
"

# ╔═╡ 26b00d28-9884-4e0d-b4b3-be4e06edb6e4
md"
The purpose of the code is to build an histogram of relative velocities, as a function of the distances between the galaxies. Thus, our output will be a matrix with 2 columns, the first will contain the number of pairs of galaxies found at each distance. The second column will, initially, sum the relative velocities of the galaxies, and finally return the average relative velocity (by dividing the second column by the first column).
"

# ╔═╡ fd5d7c8b-653c-4a49-b62e-4b8610cce854
md"
### Updating the histogram

First, we need to define the function that will be mapped for each pair of galaxies found to be within the cutoff. This function has to update the histogram, by adding the to the first column of the histogram the pair count, and the relative velocity to the corresponding element of the second column. 

This is the function:
"

# ╔═╡ 3b4b4e0c-7616-457f-b538-f7f739588988
function up_histogram!(i,j,d2,vel,binstep,hist)
    bin = Int(div(sqrt(d2),binstep,RoundUp))
    if bin <= size(hist,1)
        hist[bin,1] += 1
        hist[bin,2] += norm(vel[i] - vel[j])
    end
    return hist
end

# ╔═╡ dd89f42b-b3b2-4989-9095-be09269da22b
md"
The function receives the indexes of the particles, which are important because they will be used to index the velocity array. 

It receives the squared distance, `d2`, which comes from the `CellListMap` interface. Additionally, two parameters will be closed over by the function: the distance step of the histogram (`binstep`) and the array of velocities `vel`, containing the velocity of each galaxy. 

Finally, the `hist` variable contains the `2 x N` histogram, which is the output to be updated, as described above.

Let us experiment with this function. First, we define an empty histogram, similar to the one we will use in the actual calculation:

"

# ╔═╡ 725a1663-31e5-47a4-83eb-fd41fc69fb25
hist_test = zeros(10,2)

# ╔═╡ b626a93c-6ced-42c9-9a5f-f4c1645c7d6e
md"
Let us imagine that we have 5 galaxies, such that the velocity array has 5 components, each of which is a vector of 3 dimensions:
"
 

# ╔═╡ 0b3a5604-079a-4ec2-abee-6acfb7d6511c
vel_test = rand(SVector{3,Float64},5)

# ╔═╡ b5572781-0ac0-4c3f-875f-20aa57c010ef
md"
Using a `binstep` of `0.5`, and assuming that the pair of galaxies with `i=2` and `j=5` have a distance of `3.0`, we have:
"

# ╔═╡ acdb8062-5edb-477a-8baf-b07345f630b9
up_histogram!(2,5,3.0^2,vel_test,0.5,hist_test)

# ╔═╡ f603c345-c0f3-4475-b6f2-8c0971a8e503
md"
Which gives us the updated `hist_test` histogram, with the count of the number of pairs of galaxies on each bin, and the sum of the corresponding relative velocities.
"

# ╔═╡ 0ecf56d5-71d6-41da-bdfa-93c377b392ee
md"
### Generating a galactic system

The internal `CellListMap.xgalactic` function generates positions and a box with the typical density and cutoff used in this type of calculations. Here, we will generate a system with 1 million \"galaxies\" (with random positions and velocities):
"

# ╔═╡ 31f3e90a-ac38-4c9d-885d-3c0a84d5307d
pos, box = CellListMap.xgalactic(10^6);

# ╔═╡ c8de9825-6a2d-44c2-ae22-749a544da7a4
md"
The box has 43.7 mega-parsecs of side, and a typical cutoff of 5 mega-parsecs is used:
"

# ╔═╡ 6694c59d-d165-4f1b-b07d-a7fc55272c8f
box

# ╔═╡ 4407e6e9-3ffa-42cc-aef0-2ff0b2a923e0
md"
The positions are just randomly generated in such a box:
"

# ╔═╡ 0cc6814b-19cd-40c4-8d03-fb77b5ae794c
pos

# ╔═╡ e370be09-3d13-4a4b-b51d-c7ebdff8fad8
md"
Let us generate random velocities:
" 

# ╔═╡ eac11d86-54f1-4ae0-b451-e3de9e9648c7
vel = [ rand(SVector{3,Float64}) for _ in pos ]

# ╔═╡ 60609fdf-e96c-4b46-ae52-ca9014ee2622
md"
### Computing the histogram

With thge above that, we only need to initialize the histogram, compute the cell lists and map the `up_histogram!` function:

"

# ╔═╡ d4b5f4df-3200-46da-bcf9-a5adf6ece67f
hist = zeros(10,2)

# ╔═╡ 75af8321-b8df-46e2-9ae5-a651a790e4a8
md"
Since our cutoff is 5.0, and our histogram has 10 bins, we use:
"

# ╔═╡ be2c146b-53ac-4b7a-96fb-d752c5c57c42
binstep = 0.5 # cutoff is 5 (thus the histogram has 10 positions)

# ╔═╡ 804823db-8c12-4cd6-b2ed-928b96de61f1
md"
Let us construct the cell lists:
"

# ╔═╡ bc9d6fa4-2617-4161-915d-ca869a57120f
cl = CellList(pos,box)

# ╔═╡ 34a0edbf-00a7-4b73-aeb2-c56bee01e4ec
md"
And, finally, map the `up_histogram!` through the cell lists:
"

# ╔═╡ 77363ec5-df4e-4c7c-8406-e07587948f97
    map_pairwise!(
        (x,y,i,j,d2,out) -> up_histogram!(i,j,d2,vel,binstep,hist),
        hist, box, cl,
    )

# ╔═╡ 3667b3e4-51d5-4ac3-988d-94a69a5b5b15
md"
Since what we want is the *average* relative velocity in each bin, we normalize the second column of the histogram by the first column:
"

# ╔═╡ a6ea164f-4026-4d19-adc6-a9d730218898
hist[:,2] ./ hist[:,1]

# ╔═╡ 41d778c5-f319-4152-8299-bcff8769b7ea
md"
Since the velocities were randomly generated, the histogram is simply flat.
"

# ╔═╡ 3da1c372-bc51-4970-9424-83302ec0909d
md"
In the code block, we have put all these steps inside a single function, which is actually recommended to avoid type instabilities. Thus, the function of the code block is:
"

# ╔═╡ c5346c37-c977-4355-8415-a00786b2c3cd
function pairwise_velocities(N)
    hist = zeros(10,2)
    pos, box = CellListMap.xgalactic(N)
    vel = [ rand(SVector{3,Float64}) for _ in pos ]
    binstep = 0.5 # cutoff is 5.0
    cl = CellList(pos,box)
    map_pairwise!(
        (x,y,i,j,d2,out) -> up_histogram!(i,j,d2,vel,binstep,hist),
        hist, box, cl,
    )
    return hist[:,2] ./ hist[:,1]
end

# ╔═╡ 94ad5c00-f119-4429-8f0f-2c8c1dfd6b3d
md"
Which is simply run with:
"

# ╔═╡ eb59bf39-fa68-40a5-911a-7664057ce0fe
pairwise_velocities(10^6)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CellListMap = "69e1c6dd-3888-40e6-b3c8-31ac5f578864"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
CellListMap = "~0.7.13"
StaticArrays = "~1.4.3"
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
# ╟─b68cfd92-d1b3-4f33-b9bf-8f72b101a5b9
# ╟─b8f7a5e0-2ed6-4df5-8c76-aeff9c113ab6
# ╠═562abb8e-37ed-4fac-8293-c47e7d1d6203
# ╟─34f87c93-ef98-41be-97f1-7fda2bfb5126
# ╠═3e39c829-f965-46a0-86e5-d08e4a101fdb
# ╟─26b00d28-9884-4e0d-b4b3-be4e06edb6e4
# ╟─fd5d7c8b-653c-4a49-b62e-4b8610cce854
# ╠═3b4b4e0c-7616-457f-b538-f7f739588988
# ╟─dd89f42b-b3b2-4989-9095-be09269da22b
# ╠═725a1663-31e5-47a4-83eb-fd41fc69fb25
# ╟─b626a93c-6ced-42c9-9a5f-f4c1645c7d6e
# ╠═1ed49e14-1cca-4aae-b542-842a90e05ab9
# ╠═0b3a5604-079a-4ec2-abee-6acfb7d6511c
# ╟─b5572781-0ac0-4c3f-875f-20aa57c010ef
# ╠═acdb8062-5edb-477a-8baf-b07345f630b9
# ╟─f603c345-c0f3-4475-b6f2-8c0971a8e503
# ╟─0ecf56d5-71d6-41da-bdfa-93c377b392ee
# ╠═31f3e90a-ac38-4c9d-885d-3c0a84d5307d
# ╟─c8de9825-6a2d-44c2-ae22-749a544da7a4
# ╠═6694c59d-d165-4f1b-b07d-a7fc55272c8f
# ╟─4407e6e9-3ffa-42cc-aef0-2ff0b2a923e0
# ╠═0cc6814b-19cd-40c4-8d03-fb77b5ae794c
# ╟─e370be09-3d13-4a4b-b51d-c7ebdff8fad8
# ╠═eac11d86-54f1-4ae0-b451-e3de9e9648c7
# ╟─60609fdf-e96c-4b46-ae52-ca9014ee2622
# ╠═d4b5f4df-3200-46da-bcf9-a5adf6ece67f
# ╟─75af8321-b8df-46e2-9ae5-a651a790e4a8
# ╠═be2c146b-53ac-4b7a-96fb-d752c5c57c42
# ╟─804823db-8c12-4cd6-b2ed-928b96de61f1
# ╠═bc9d6fa4-2617-4161-915d-ca869a57120f
# ╟─34a0edbf-00a7-4b73-aeb2-c56bee01e4ec
# ╠═77363ec5-df4e-4c7c-8406-e07587948f97
# ╟─3667b3e4-51d5-4ac3-988d-94a69a5b5b15
# ╠═a6ea164f-4026-4d19-adc6-a9d730218898
# ╟─41d778c5-f319-4152-8299-bcff8769b7ea
# ╟─3da1c372-bc51-4970-9424-83302ec0909d
# ╠═c5346c37-c977-4355-8415-a00786b2c3cd
# ╟─94ad5c00-f119-4429-8f0f-2c8c1dfd6b3d
# ╠═eb59bf39-fa68-40a5-911a-7664057ce0fe
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
