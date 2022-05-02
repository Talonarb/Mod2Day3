# This is a function that takes the input of an edge dictionary and creates an adjacency matrix. No longer used, but useful code. 
def dict_to_matrix(edge_dict):
    global edge_matrix
    keys=sorted(edge_dict.keys())
    size=len(keys)
    edge_matrix = [ [0]*size for i in range(size) ]
    print(edge_dict.keys())
    print(size)
    for a,b in [(keys.index(a), keys.index(b)) for a, row in edge_dict.items() for b in row]:
        print(a, b)
        edge_matrix[a][b] = 2 if (a==b) else 1
    print(edge_matrix)
    print(keys)
    for i in range(size): 
        for j in range(size): 
            print(edge_matrix[i][j], end = " ")
        print() 

