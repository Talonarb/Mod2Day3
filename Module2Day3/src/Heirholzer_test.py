def pairs_to_dict(pair_list, dictionary_to_be_made): 
    for a, b in pair_list: # loops through list of pairs 
        dictionary_to_be_made.setdefault(a, []).append(b) # makes a key of the "a" position in pairs with value "b" 
    return dictionary_to_be_made # and now we have our dictionary

def printCircuit(adj, starting_node):
   global circuit_for_fncs

    # adj represents the adjacency list of
    # the directed graph
   if len(adj) == 0:
        return(1)
   
# Maintain a stack to keep vertices
# We can start from any vertex, here we start with a given value. 
   curr_path = [starting_node]
   #print("This is the function's call of starting location. It should change.")
   #print(starting_node)
# list to store final circuit
   circuit = []
   
   while curr_path:
   
        curr_v = curr_path[-1]
           
# If there's remaining edge in adjacency list of the current vertex 
        if adj[curr_v]:
  
# Find and remove the next vertex that is adjacent to the current vertex
            next_v = adj[curr_v].pop()
# Push the new vertex to the stack
            curr_path.append(next_v)
# back-track to find remaining circuit, pop function of a recursive element 
        else:
# Remove the current vertex and put it in the circuit
            circuit.append(curr_path.pop())

# we've got the circuit, now print it in reverse
   reversed_circuit = circuit[::-1]
   circuit_for_fncs = reversed_circuit
   return(circuit_for_fncs)
    ### This is the end of the function

def visualize_Eulerian_tour(edges, edge_dict, starting_node):
# Initialize a list of empty "[]", which act as a node, outgoing connections symbolized by any integers between the brackets, index of the 
# closed brackets corresponding to the node list 
    node_list = []
# As the node list generated from kmerization is not necessarily complete, fill the nodes from the edge list, which will always have
# every node visited in at least pair as either source or destination. No need for redudant nodes, but want list traversal, so using set
# followed by list converts the list to only unique values, and from there is converted back into a list. 
    for src, dest in edges:
        node_list.append(src)
        node_list.append(dest)
    # Remove duplicates vs set function, is it faster? Will time it next time. 
    node_list = list(set(node_list))

    # This sorts the node list, helpful for making sure you can always reproduce results with same starting node 
    node_list = sorted(node_list)

    circuit_for_test = [ [] for _ in node_list ]

    for i in node_list:

    # Raise an exception for when a node has no outgoing connections, such that starting node can be set to a different node
        try:
            index_to_append = int(node_list.index(i))
            for j in edge_dict[i]:
    # CNCT in this case is the connection 
                cnct = int(node_list.index(j))
                circuit_for_test[index_to_append].append(cnct)
        except KeyError: 
                print("This message occurs when there are no outgoing links from a node in a graph, the function will continue iterating for other nodes.")
    circuit_for_fncs = printCircuit(circuit_for_test, starting_node)
    start_node = node_list[(circuit_for_fncs[int(starting_node)]+1)]
    start_node = str(start_node)

    list_of_nodes_to_be_converted = []
    for i in circuit_for_fncs[1:]:
        list_of_nodes_to_be_converted.append(node_list[i])
    list_of_converted_nodes = [x[-1] for x in list_of_nodes_to_be_converted]
    list_of_converted_nodes_to_string = ''.join(map(str, list_of_converted_nodes))
    contig_from_path = start_node + list_of_converted_nodes_to_string
    print("This is the assembled contig: ", contig_from_path)
    return(contig_from_path)