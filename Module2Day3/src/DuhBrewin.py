# De Bruijn graph creation, function named with the intent to describe the phoenetics, while also being mildly humorous.
# Majority of this function by https://www.youtube.com/watch?v=f5kgmqcwb8M , with minor edits and comments. 

def DuhBrewin(txt_to_be_kmerized, k):
    edges = []
    nodes = set()
# Loop through the text to identify individual sequences
    for seq in txt_to_be_kmerized :
        for i in range(len(seq) - k + 1): # Loop through sequences to create kmers
            edges.append((seq[i : i + k - 1], seq[i + 1: i + k])) # append the list of edges using the connections between k-1mers as kmers
# Left and right k-1mers  added to edge list #
        nodes.add(seq[i:i+k-1]) # and establish the nodes of k-1mers
        nodes.add(seq[i+1:i+k]) # still k-1mer, just moved over one position 
    edges.sort
    return nodes, edges #no sense in doing this without returning the values

def graph_representation(sequences, k):
    # This function creates a graphical representation of the De Bruijn graph created by the DuhBrewin function, and outputs it for easier
    # interpretation.
    nodes, edges = DuhBrewin(sequences, k)  #Inputs for the nodes/edges in the graph to be ##
    dot_str = 'digraph "DeBruijn Graph" {\n' # Specifies directed graph, DeBruijn in this case #  
    for node in nodes:
        dot_str += ' %s [label="%s"] ;\n ' % (node, node) # Adds nodes and their respectives labels #
    for src, dst in edges:
        dot_str += ' %s -> %s ;\n ' % (src, dst) # As above, no labels necessary in this case
    return dot_str + '}\n' 
    