from IntComplex import IntComplex
from scipy.spatial.distance import cdist



def get_interaction_vertex(inter):
    res = []
    for i in inter:
        if i!=-1 and i!=-2:
            res.append(i)
    return res

def get_intcomplex(point,weight,cutoff,add_double=False):
    '''
    point : point coordinate
    weight : point weight
    '''
    D = cdist(point, point, metric="euclidean")
    
    C1,C2,C3,C4 = [],[],[],[]
    
    # 1-interaction
    for i in range(len(point)):
        tmp = [ [i],0 ]
        C1.append(tmp)
    
    # 2-interaction
    for i in range(len(point)):
        for j in range(len(point)):
            if i!=j:
                dis = D[i,j]
                if dis<=cutoff:
                    if weight[i]<weight[j]:
                        tmp = [ [-2,j,i,-1],dis ]
                        C2.append(tmp)
                    elif weight[i]>weight[j]:
                        tmp = [ [-2,i,j,-1],dis ]
                        C2.append(tmp)
                    else:
                        if add_double:
                            tmp1 = [ [-2,j,i,-1],dis ]
                            C2.append(tmp1)
                            tmp2 = [ [-2,i,j,-1],dis ]
                            C2.append(tmp2)
    
    # 3-interaction
    for i in range(len(C1)):
        int1,d1 = C1[i]
        p1 = int1[0]
        w1 = weight[p1]
        for j in range(len(C2)):
            int2,d2 = C2[j]
            p2,p3 = int2[1],int2[2]
            if p1==p2 or p1==p3:
                continue
            w2 = (weight[p2]+weight[p3])/2
            if D[p1,p2]<=cutoff and D[p1,p3]<=cutoff:
                d = max(D[p1,p2],D[p1,p3],d2)
                if w1<w2:
                    tmp = [ [-2]+int2+int1+[-1],d ]
                    C3.append(tmp)
                elif w1>w2:
                    tmp = [ [-2]+int1+int2+[-1],d ]
                    C3.append(tmp)
                else:
                    if add_double:
                        tmp1 = [ [-2]+int2+int1+[-1],d ]
                        C3.append(tmp1)
                        tmp2 = [ [-2]+int1+int2+[-1],d ]
                        C3.append(tmp)
    '''
    # 4-interaction
    # (1,3),(3,1)
    for i in range(len(C1)):
        int1,d1 = C1[i]
        p1 = int1[0]
        w1 = weight[p1]
        for j in range(len(C3)):
            int3,d3 = C3[j]
            p2,p3,p4 = get_interaction_vertex(int3)
            if p1==p2 or p1==p3 or p1==p4:
                continue
            w2 = (weight[p2]+weight[p3]+weight[p4])/3
            if D[p1,p2]<=cutoff and D[p1,p3]<=cutoff and D[p1,p4]<=cutoff:
                d = max(D[p1,p2],D[p1,p3],D[p1,p4],d3)
                if w1<w2:
                    tmp = [ [-2]+int3+int1+[-1],d ]
                    C4.append(tmp)
                elif w1>w2:
                    tmp = [ [-2]+int1+int3+[-1],d ]
                    C4.append(tmp)
                else:
                    tmp1 = [ [-2]+int3+int1+[-1],d ]
                    C4.append(tmp1)
                    tmp2 = [ [-2]+int1+int3+[-1],d ]
                    C4.append(tmp2)
    
    # (2,2)
    for i in range(len(C2)):
        int1,d1 = C2[i]
        p1,p2 = int1[1],int1[2]
        w1 = (weight[p1]+weight[p2])/2
        for j in range(len(C2)):
            if i!=j:
                int2,d2 = C2[j]
                p3,p4 = int2[1],int2[2]
                w2 = (weight[p3]+weight[p4])/2
                if p1==p3 or p1==p4 or p2==p3 or p2==p4:
                    continue
                if D[p1,p3]<=cutoff and D[p1,p4]<=cutoff and D[p2,p3]<=cutoff and D[p2,p4]<=cutoff:
                    d = max(D[p1,p3],D[p1,p4],D[p2,p3],D[p2,p4],d1,d2)
                    if w1>w2:
                        tmp = [ [-2]+int2+int1+[-1],d ]
                        C4.append(tmp)
                    elif w1<w2:
                        tmp = [ [-2]+int1+int2+[-1],d ]
                        C4.append(tmp)
                    else:
                        tmp1 = [ [-2]+int2+int1+[-1],d ]
                        C4.append(tmp1)
                        tmp2 = [ [-2]+int1+int2+[-1],d ]
                        C4.append(tmp2)
    '''
    #print('interaction number:',len(C1),len(C2),len(C3),len(C4))
    C = C1+C2+C3+C4
    C = sorted(C,key=lambda x:( x[1]+len(x[0])/1000.0  ))
    return C

def get_intcomplex_feature(point,weight,cutoff,typ):
    intc = get_intcomplex(point, weight, cutoff)
    #print('IntComplex is ok',len(intc))
    if typ=='multi':
        dim = 3
        complex1 = IntComplex(intc,dim)
        complex1.get_persistence()
        h1 = complex1.diagram['1']
        h2 = complex1.diagram['2']
        
        P = 2
        layer = complex1.get_layer_persistence(P)
        layer1 = layer[1]['2']
        
        P = 3
        layer = complex1.get_layer_persistence(P)
        layer2 = layer[1]['1']
        layer3 = layer[1]['2']
        
        return h1,h2,layer1,layer2,layer3
    elif typ=='pair':
        dim = 2
        complex1 = IntComplex(intc,dim)
        complex1.get_persistence()
        h1 = complex1.diagram['1']
        return h1
        
