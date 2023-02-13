
function gen_data()
    # We generate some test data similar to the input files 
    
    # The first one will have the paths for each SeqX: nodes and a +/- sign indicating their orientation
    open("paths.txt", "w+") do h1
        paths = "Seq1\t100+ 30+ 5+ 1- 5+\nSeq2\t200- 400- 100+ 30+ 18+ 9-\nSeq3\t100+ 30+ 5+ 8-\nSeq4\t100+ 13 14 5+ 1- 5+\n"
        write(h1, paths)
    end

    # The other file has the lengths of each of the nodes
    open("lens.txt", "w+") do h2
        lens = "100 10\n30  100\n5   1\n1   4\n200  80\n400 60\n18  2\n9    3\n8    6\n13   11\n14  9\n"
        write(h2, lens)
    end
end

function index_lengths(len_file::String)
    # This will create a mapping of the node_id > node size
    len_index = Dict{Int64, Int64}()
    for line in eachline(len_file)
        node, len = split(line)
        len_index[parse(Int64, node)] = parse(Int64, len)
    end 
    return len_index
end

function parse_path(path::SubString{String})
    # To convert the path to a numerical vector 
    numbers = Vector{Int64}()
    for node in split(path, " ")
        #println("N:", node)
        multiply = node[end] == '+' ? 1 : -1 
        # We store the sign to keep track of the traversal direction 
        # see "align"
        node = parse(Int64, view(node, 1:length(node)-1)) * multiply
        push!(numbers, node)
    end
    return numbers
end

function get_sizes(nodes::Vector{Int64}, len_index::Dict{Int64,Int64})
    # This will return an array the size of nodes with the sizes for each node
    sizes = Vector{Int64}(undef, length(nodes))
    for (i, node) in enumerate(nodes)
        # As we stored the traversal direction as the integer sign 
        # we have to take the "abs" to query the length in the index
        sizes[i] = len_index[abs(node)]
    end 
    return sizes
end

function extract_queries(f::String, queries::Vector{String}, len_index::Dict{Int64, Int64})
    # This will just iterate through the input file to find the queries. Then return 
    # the query ids, query paths and query path sizes (i.e. sizes for each node)
    query_ids = Vector{String}()
    query_paths = Vector{Vector{Int64}}()
    query_node_lens = Vector{Vector{Int64}}()

    h = open(f, "r")
    for line in eachline(h)
        # Split the line
        identifier, path = split(line, "\t")

        # Check if this is a query
        if identifier in queries

            # Extract the numbers (node ids) from the line
            node_path = parse_path(path)
            # Also get the sizes for each of the nodes
            node_path_sizes = get_sizes(node_path, len_index)
            
            # Store them 
            push!(query_node_lens, node_path_sizes)
            push!(query_ids, identifier)
            push!(query_paths, node_path)
        end
    end
    return query_ids, query_paths, query_node_lens
end

function index_ref_path(number_path::Vector{Int64})
    # This will create an hashmap storing the location of each node in the path 
    # We want to map the absolute node ids to the postion
    # while preserving its traversal direction (i.e. sign of node)
    # therefore we move the sign from the node id to the position (i * sign(node))
    d = Dict{Int64, Vector{Int64}}()
    for (i, node) in enumerate(number_path)
        if haskey(d, abs(node))
            push!(d[abs(node)], i * sign(node))
        else 
            d[abs(node)] = [i * sign(node)]
        end
    end
    return d
end 

function extend_seed(q::Array{Int64}, r::Array{Int64}, qi::Int64, ri::Int64, init_sign::Bool)
    # This will extend a seed (see description above the function below)
    q_len = length(q)
    r_len = length(r)
    r_dir = if init_sign 1 else -1 end
    
    # Extend to the left
    lqi = qi
    lri = ri
    @inbounds while lqi > 0 && lqi <= q_len && lri > 0 && lri <= r_len
        q_node = q[lqi]
        r_node = r[lri]
        # We abort the walk when:
        # - the node between the query and ref are not the same anymore 
        # - the signs become differnet (as the nodes are traverssed differently then)
        if abs(q_node) != abs(r_node) || (q_node ⊻ r_node >= 0) != init_sign
            break
        end
         lri -= r_dir
         lqi -=1
    end
    lqi +=1
    
    # Extend to the right
    rqi = qi
    rri = ri
    @inbounds while rqi > 0 && rqi <= q_len && rri > 0 && rri <= r_len
        q_node = q[rqi]
        r_node = r[rri]
        if abs(q_node) != abs(r_node) || (q_node ⊻ r_node >= 0) != init_sign
            break
        end
         rri += r_dir
         rqi +=1
    end
    rqi -=1
    
    return lqi, rqi
end

#=
For extension, we have to define when a walk between a query q and reference r cohere Let Wq = (q1, q2, ... , qm) 
be the node walk for q,  Wr = (r1, r2, ... , rn) be the node walk for r, and S(n) = n < 0 extract the sign, 
i.e. traversal direction, of a node n.  Then we define a match of l + 1 nodes between Wq and Wr if for some i,j
    (qi, qi +1, ... qi+l) = (rj, rj+1, ... , rj+l)
and
	(s(qi), s(qi + 1), ... , s(qi + l) = (s(rj), s(rj + 1), ... , s(rj +l)
Or, when a sequence is present as reverse complement compared to the the query if
    (qi, qi +1, ... qi+l) = (rj, rj-1, ... , rj-l)
and
	(s(qi), s(qi+1), ... , s(qi +l )) = (-s(rj), -s(rj-1), ... , -s(rj-l))

Given any seed (qi=rj), we can apply the above rules to bidirectionally extend it (compare function)
Extending every seed would be equivalent to finding all Maximum Exact Matches(MEM) between q and r. 
As we are only interested in LMEMs we apply a seed selection criteria. Let M(i,i+l) be a MEM of length l+1 spanning nodes i to i+l in Wq.
 Every seed at position i <pi +l will result in a shorter or equally long MEM. The next opportunity to find a longer MEM is seeding at i + l + 1.  
 While extending a seed at i + l + 1 can never yield an alignment over the entirety of the previous MEM, it can align back up to i + 1. 
 When a new MEM overlaps with the previous MEM we have to assign the overlapping nodes to the longest of the two, the LMEM. 
 Longest refers to the covered nucleotides, not the number of nodes. Using the seed selection criteria, we move along Wq, bidirectionally extend the seeds, and resolve overlaps. 
 The result is a vector the same length of Wq recording the LMEM size for each node in the query
=#
function pairwise_aligner(q::Array{Int64}, qls::Array{Int64}, r::Array{Int64})
    ref_index = index_ref_path(r)
    
    qi = 1
    q_len = length(q)
    coverage =  zeros(Int32, q_len)
    
    while qi <= q_len
        # The current query node
        q_node = q[qi]
       
        # How far can we extend the node to the left and right
        max_left = qi
        max_right = qi
        max_extend = 0

        # Check if this query node is in our index
        if haskey(ref_index, abs(q_node))
            # Try to extend each seed 
            for ri in ref_index[abs(q_node)]
                # Are the signs the same or not 
                init_sign = q_node ⊻ ri >= 0 
                # Extend the seed to the left and right as far as possible 
                extend_left, extend_right = extend_seed(q, r, qi, abs(ri), init_sign)
                # Is this extenstion longer than the previous max
                if (extend_right - extend_left > max_extend) 
                    max_left = extend_left
                    max_right = extend_right
                end 
            end
            # Now we check how many nucloetides the segment is
            # translating the #nodes to the length in nucleotides
            nts = sum(view(qls, max_left:max_right))

            # Fil this in in the overal matches, we replace a previous 
            # match if the current one is longer regarding the nucleotides
            for jj in max_left:max_right
                if coverage[jj] < nts
                    coverage[jj] = nts
                end
            end
        end

        # Now we to the end of the maximum extension we previously did
        qi = max_right 
        qi +=1
    end 
    return coverage                
end

function align_dispatch(query_ids::Vector{String}, query_paths::Vector{Vector{Int64}}, query_node_lens::Vector{Vector{Int64}}, path_file::String)
    
    # Iterate over each query 
    for i in eachindex(query_ids)
        q_id = query_ids[i]
        q_path = query_paths[i]
        q_node_lens = query_node_lens[i]
        println("Alignments for ", q_id)
        println("Q: ", q_path)

        # To decide the best over all pairwise ones 
        overall_sizes = zeros(Int64, length(q_path))
        origin = Vector{String}(undef, length(q_path))

        # Iterate over each reference
        h = open(path_file, "r")
        for line in eachline(h)
            r_id, path = split(line, "\t")
            r_path = parse_path(path)
            if r_id != q_id # Exclude self-self
                pairwise_res = pairwise_aligner(q_path, q_node_lens, r_path)
                println("> ",r_id)
                println(pairwise_res)
                # Check if we have better matches than before, per position 
                for i in eachindex(q_path)
                    if pairwise_res[i] > overall_sizes[i]
                        overall_sizes[i] = pairwise_res[i]
                        origin[i] = r_id
                    end
                end
            end
        end

        println("Merged results")
        print("Sizes: ")
        println(overall_sizes)
        print("Origins: ")
        println(origin)

    end

end

gen_data()
length_index = index_lengths("lens.txt")
query_ids, query_paths, query_node_lens = extract_queries("paths.txt", ["Seq1"], length_index)
align_dispatch(query_ids, query_paths, query_node_lens, "paths.txt")
