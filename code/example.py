from utils import get_intcomplex
from IntComplex import IntComplex





def main():
    point = [ 
              [0,0],
              [2,0],
              [2,1],
              [0,1]
     ]
    weight = [1,1,3,2  ]
    cutoff = 4
    
    # construct intcomplex
    intc = get_intcomplex(point, weight, cutoff)
    
    # persistent homology
    dim = 3
    complex1 = IntComplex(intc,dim)
    complex1.get_persistence()
    h1 = complex1.diagram['1']
    h2 = complex1.diagram['2']
    print('persistent homology:\n',complex1.diagram)
    
    # persistent layer homology
    P = 2
    layer = complex1.get_layer_persistence(P)
    layer1 = layer[1]
    print('persistent layer homology:\n',layer1)
        
        
main()
        
