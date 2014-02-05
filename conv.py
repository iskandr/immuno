"""
Also, the basic sketch of my idea for a convolutional-neural-net inspired approach to this immunogenicity problem:

Each amino acid has some relevant chemical / physical properties which affects its binding properties, but we don't know which properties those are (a wide range all seem to work equally well). So, learn a stack of encoding dictionaries that transform each amino acid in a sequence to several numbers. 

For example (with two encoding dictionaries): 
AAGY -> [0.1 0.1 -0.3 -0.9]
              [ 1.7 1.7 1.2   1.9]

From this point, use small spatial  filters (2-4 bases wide) that can see all the encoded channels.

For example, if there's only one spatial filter of width 3, called f here, the output would be:
               [0.1 0.1 -0.3 -0.9] --> [f([[0.1, 0.1, -0.3], [1.7 1.7 1.2]], f([0.1, -0.3, -0.9], [1.7, 1.2, 1.9]]) -> [0.8, 2.7]
              [ 1.7 1.7 1.2   1.9] 

Downsample this spatially filtered result by taking the max of a small window (i.e. max of 2 neighboring values). 
"""
