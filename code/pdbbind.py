import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import sys
from utils import get_intcomplex_feature



def atmtyp_to_ele( st ):
    # C,N,O,S,P,H,F,Cl,Br,I
    st = st.strip()
    if len(st) == 1:
        return st
    elif st[0] == 'H':
        return 'H'
    elif st[0] == 'N':
        return 'N'
    elif st[0] == 'O':
        return 'O'
    elif st[0] == 'S':
        return 'S'
    elif st[0] == 'P':
        return 'P'
    elif st[0] == 'C' and st[0:2] not in ['Cl','CL']:
        return 'C'
    elif st[0] == 'F':
        return 'F'
    elif st[0] == 'I':
        return 'I'
    elif st[0:2] in ['Cl', 'CL','cl']:
        return 'CL'
    elif st[0:2] in ['BR','Br','br']:
        return 'BR'
    elif st[1] in ['H']:
        return 'H'
    else:
        print(st, 'Not in dictionary')
        return

residue_to_one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

class Atom:
    def __init__(self,atype,resname,chain,resid,coord,charge=None):
        self.AType = atype
        self.Coord = coord
        self.ResName = resname
        self.ResId = resid
        self.Chain = chain
        self.charge = charge

class ProteinLigand:
    def __init__(self,pdb,cutoff,pdb_folder=None,pdb_feature_folder=None):
        self.pdb= pdb
        self.pdb_folder = pdb_folder
        self.pdb_feature_folder = pdb_feature_folder
        
        self.Protein_Atoms = []
        self.Protein_AtomCoord = []
        self.Protein_Weight = []
        self.Ligand_Atoms = []
        self.Ligand_AtomCoord = []
        self.Ligand_Weight = []
        self.P_Index_List = []
        self.L_Index_List = []
        
        
        self.read_atom_from_pdb()
        self.set_euclidean_index_list()
        self.get_int_feature(cutoff)
        s
        
    
    def read_atom_from_pdb(self):
        electronegativity = {
            "H": 2.20,
            "C": 2.55,
            "N": 3.04,
            "O": 3.44,
            "F": 3.98,
            "P": 2.19,
            "S": 2.58,
            "CL": 3.16,
            "BR": 2.96,
            "I": 2.66,
            "Se": 2.55,
            'B': 0.0,
        }
    
        # protein
        filename1 = self.pdb_folder + '/'+self.pdb+'/'+self.pdb+'_pocket.pdb'
        filename2 = self.pdb_folder + '/'+self.pdb+'/'+self.pdb+'_ligand.mol2'
        
        f = open(filename1)
        contents = f.readlines()
        f.close()
        for line in contents:
            if line[0:4]=='ATOM':
                atom = Atom(atype=atmtyp_to_ele(line[12:16]), resname=residue_to_one_letter[line[17:20]], chain=line[21], resid=int(line[22:26]),
                            coord=[float(line[30:38]),float(line[38:46]),float(line[46:54])],
                             )
                self.Protein_Atoms.append(atom)
                self.Protein_Weight.append(electronegativity[atmtyp_to_ele(line[12:16])]+4)
                self.Protein_AtomCoord.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        self.Protein_AtomCoord = np.array(self.Protein_AtomCoord)
        self.Protein_Weight = np.array(self.Protein_Weight)
        
        
        # ligand
        f = open(filename2)
        contents = f.readlines()
        f.close()
        
        start = 0
        end = 0
        for jj in range(len(contents)):
            if contents[jj][0:13]=='@<TRIPOS>ATOM':
                start = jj + 1
                continue
            if contents[jj][0:13]=='@<TRIPOS>BOND':
                end = jj - 1
                break
        for kk in range(start,end+1):
            if contents[kk][8:17]=='thiophene':
                print('thiophene',kk)
            line = contents[kk]
            atom = Atom(atype=atmtyp_to_ele(line[47]), resname='ligand', chain='ligand', resid='ligand',
                            coord=[float(line[16:26]),float(line[26:36]),float(line[36:46])],
                             )
            self.Ligand_Atoms.append(atom)
            self.Ligand_Weight.append(electronegativity[atmtyp_to_ele(line[47])])
            self.Ligand_AtomCoord.append([float(line[16:26]),float(line[26:36]),float(line[36:46])])
        self.Ligand_AtomCoord = np.array(self.Ligand_AtomCoord)
        self.Ligand_Weight = np.array(self.Ligand_Weight)
    
    def set_euclidean_index_list(self):
        ele2index = {'C':0, 'N':1, 'O':2, 'S':3, 'P':4,'H':5}
        
        # protein C, N, O, S
        self.P_Index_List = [ [] for _ in range(4) ]
        for idx in range(len(self.Protein_Atoms)):
            atom = self.Protein_Atoms[idx]
            if atom.AType in ['C','N','O','S']:
                self.P_Index_List[ele2index[atom.AType]].append(idx)
            
        
        # ligand C, N, O, S, P,H
        self.L_Index_List = [ [] for _ in range(6) ]
        for idx in range(len(self.Ligand_Atoms)):
            atom = self.Ligand_Atoms[idx]
            if atom.AType in ele2index:  
                self.L_Index_List[ele2index[atom.AType]].append(idx)
    
    def get_int_feature(self,cutoff):
        step = 0.5
        N = int(cutoff/step)
        
        
        
        # pair
        pair = np.zeros((24,N))
        
        for i in range(4):
            index1 = self.P_Index_List[i]
            atom1 = self.Protein_AtomCoord[index1]
            w1 = self.Protein_Weight[index1]
            for j in range(6):
                index2 = self.L_Index_List[j]
                atom2 = self.Ligand_AtomCoord[index2]
                w2 = self.Ligand_Weight[index2]
                
                atom = np.vstack([atom1,atom2])
                w = w1.tolist()+w2.tolist()
                
                
                pair_h1 = get_intcomplex_feature(atom,w,cutoff,'pair')
                now = i*6+j
                for bar in pair_h1:
                    bir,dea = bar
                    if dea == -1:
                        dea = cutoff
                    n = int(dea/step)
                    n = min(N,n)
                    pair[now,0:n] = pair[now,0:n]+1
        
                    
                
        
        # multi
        multi = np.zeros((24,N,5))
        lig = [ [0,1],[0,2],[0,5],[1,2],[1,5],[2,5] ]
        for i in range(4):
            index1 = self.P_Index_List[i]
            atom1 = self.Protein_AtomCoord[index1]
            w1 = self.Protein_Weight[index1]
            for j in range(6):
                idx1,idx2 = lig[j]
                index2 = self.L_Index_List[idx1]+self.L_Index_List[idx2]
                atom2 = self.Ligand_AtomCoord[index2]
                w2 = self.Ligand_Weight[index2]
                
                atom = np.vstack([atom1,atom2])
                w = w1.tolist()+w2.tolist()
                
                mulit_h1,multi_h2,layer1,layer2,layer3 = get_intcomplex_feature(atom,w,cutoff,'multi')
                #print(layer1['2'])
                tmp = [mulit_h1,multi_h2,layer1,layer2,layer3]
                now = i*6+j
                for k in range(5):
                    for bar in tmp[k]:
                        bir,dea = bar
                        n1 = int(bir/step)
                        if dea==-1:
                            dea = cutoff
                        n2 = int(dea/step)
                        multi[now,n1:n2,k] = multi[now,n1:n2,k]+1
                
        np.save(self.pdb_feature_folder+self.pdb+'-pair.npy',pair)
        np.save(self.pdb_feature_folder+self.pdb+'-multi.npy',multi)
        #print(pair.shape)
        #print(multi.shape)
        
    
    def get_simplicial_feature(self,cutoff):
        step = 0.5
        N = int(cutoff/step)
        
        
        
        # pair
        pair = np.zeros((24,N))
        
        for i in range(4):
            index1 = self.P_Index_List[i]
            atom1 = self.Protein_AtomCoord[index1]
            w1 = self.Protein_Weight[index1]
            for j in range(6):
                index2 = self.L_Index_List[j]
                atom2 = self.Ligand_AtomCoord[index2]
                w2 = self.Ligand_Weight[index2]
                
                atom = np.vstack([atom1,atom2])
                w = w1.tolist()+w2.tolist()
                m,n = len(atom1),len(atom2)
                
                block_dis = cdist(atom1, atom2, metric='euclidean')
                dis_m = np.full((m + n, m + n), 99.0)
                dis_m[:m,m:] = block_dis
                dis_m[m:,:m] = block_dis.T
                np.fill_diagonal(dis_m, 0.0)
                
                
                rips_complex = gudhi.RipsComplex(distance_matrix=dis_m,
                                 max_edge_length=cutoff)
                st = rips_complex.create_simplex_tree(max_dimension=1)
                diag = st.persistence()
                now = i*6+j
                for item in diag:
                    if item[0]==0:
                        bir,dea = item[1]
                        if np.isinf(dea):
                            dea = cutoff
                        n = int(dea/step)
                        n = min(N,n)
                        pair[now,0:n] = pair[now,0:n]+1
        
                    
                
        
        # multi
        multi = np.zeros((24,N))
        lig = [ [0,1],[0,2],[0,5],[1,2],[1,5],[2,5] ]
        for i in range(4):
            index1 = self.P_Index_List[i]
            atom1 = self.Protein_AtomCoord[index1]
            w1 = self.Protein_Weight[index1]
            for j in range(6):
                idx1,idx2 = lig[j]
                index2 = self.L_Index_List[idx1]+self.L_Index_List[idx2]
                atom2 = self.Ligand_AtomCoord[index2]
                w2 = self.Ligand_Weight[index2]
                
                atom = np.vstack([atom1,atom2])
                w = w1.tolist()+w2.tolist()
                
                
                rips_complex = gudhi.RipsComplex(points=atom,
                                 max_edge_length=cutoff)
                st = rips_complex.create_simplex_tree(max_dimension=2)
                diag = st.persistence()
                now = i*6+j
                for item in diag:
                    if item[0]==1:
                        bir,dea = item[1]
                        if np.isinf(dea):
                            dea = cutoff
                        n1 = int(bir/step)
                        n2 = int(dea/step)
                        multi[now,n1:n2] = multi[now,n1:n2]+1
                
                
        np.save(self.pdb_feature_folder+self.pdb+'-pair.npy',pair)
        np.save(self.pdb_feature_folder+self.pdb+'-multi.npy',multi)
        #print(pair.shape)
        #print(multi.shape)
        



def read_data(year):
    filename = './index/train_data_'+str(year)+'.txt'
    f = open(filename)
    train = eval(f.read())
    f.close()
    
    filename = './index/test_data_'+str(year)+'.txt'
    f = open(filename)
    test = eval(f.read())
    f.close()
    
    data = train+test
    print(len(data))
    
    return data




def prepare_feature(start,end):
    year = 2016
    data = read_data(year)
    
    cutoff = 8
    pdb_folder = 'data/PDBbind-'+str(year)+'/'
    pdb_feature_folder = './data/'+str(year)+'-feature/'
    
    
    for i in range(start,end):
        pdb = data[i]
        pl = ProteinLigand(pdb,cutoff,pdb_folder,pdb_feature_folder)
        
        print(i,pdb)



#start = int(sys.argv[1])
#end = int(sys.argv[2])
#prepare_feature(start,end)




def gradient_boosting(X_train,Y_train,X_test,Y_test,seed,year):
    params={'n_estimators': 4000, 'max_depth': 7, 'min_samples_split': 2,'random_state':seed,
                'learning_rate': 0.01, 'loss': 'squared_error','max_features':'sqrt','subsample':0.5}
    regr = GradientBoostingRegressor(**params)
    regr.fit(X_train,Y_train)
    pred = regr.predict(X_test)
    
    #np.save('./prediction3/' + str(year) + '-' +str(seed)+'-pred.npy',pred)
    
    pcc,_ = sp.stats.pearsonr(Y_test,pred)
    mse = mean_squared_error(Y_test, pred)
    rmse = pow(mse,0.5)
    print(f'name:{name}, seed:{seed}, PCC: {pcc:.3f}, RMSE: {rmse:.3f}')

    
    
    
def prediction(year,seed):
    train_feature,train_label,test_feature,test_label = get_feature(year)
    print(train_feature.shape,train_label.shape,test_feature.shape,test_label.shape)
    
    scaler = MinMaxScaler(feature_range=(-1, 1))
    train_feature = scaler.fit_transform(train_feature)
    test_feature = scaler.transform(test_feature)
    
    
    gradient_boosting(train_feature,train_label,test_feature,test_label,seed,year)


prediction()
