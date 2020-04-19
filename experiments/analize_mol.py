import numpy as np
import pandas as pd


def conn_type(comm):
    num_of_double = 0
    num_of_triple = 0
    for line in comm:
        if line[2] == '2':
            num_of_double += 1
        elif line[2] == '3':
            num_of_triple += 1
    if num_of_double == 0 and num_of_triple == 0:
        return 's'
    elif num_of_double == 1 and num_of_triple == 0:
        return 'd'
    elif num_of_double >= 2 and num_of_triple == 0:
        return 'w'
    elif num_of_triple > 0:
        return 't'
    else:
        return '?'

def analize_mol_file1(molekula, file):
    with open(file) as f:
        lines = f.readlines()
        num_of_atoms = int(list(lines[3].split())[0])
        num_of_edges = int(list(lines[3].split())[1])
        end = lines.index('M  END\n')
        beg = 4 + num_of_atoms
        for elem in range(1, num_of_atoms+1):
            NN = lines[3+elem].split()[3]
            num_of_neighbours = 0
            communic = [] #массив связей с конкретным элементом
            
            for i in range(end-beg):
                if int(list(lines[beg+i].split())[0]) == elem:
                    num_of_neighbours += 1
                    communic.append(list(lines[beg+i].split())[0:3])
                if int(list(lines[beg+i].split())[1]) == elem:
                    num_of_neighbours += 1
                    communic.append(list(lines[beg+i].split())[0:3])
        
            num_of_neighbours = str(num_of_neighbours)
            type_of_connection = conn_type(communic)
            name = NN + num_of_neighbours + type_of_connection
            
            if file in molekula.index:
                if name in molekula.columns:
                    if np.isnan(molekula[name][file]):
                        molekula.loc[file, name] = 1
                    else:
                        molekula.loc[file, name] += 1
                else:
                    molekula.loc[file, name] = 1
    else:
        molekula.loc[file, name] = 1


def analize_mol_file2(molekula, file):
    with open(file) as f:
        lines = f.readlines()
        num_of_atoms = int(list(lines[3].split())[0])
        num_of_edges = int(list(lines[3].split())[1])
        end = lines.index('M  END\n')
        beg = 4 + num_of_atoms
        for elem in range(1, num_of_atoms+1):
            connect = [] #массив связей с большим по номеру элементом
            
            for i in range(end-beg):
                if int(list(lines[beg+i].split())[0]) == elem and list(lines[beg+i].split())[1] > str(elem):
                    connect.append(list(lines[beg+i].split())[0:3])
                if int(list(lines[beg+i].split())[1]) == elem and list(lines[beg+i].split())[0] > str(elem):
                    connect.append(list(lines[beg+i].split())[0:3])
            for num in range(len(connect)):
                if int(connect[num][0]) == elem:
                    elem2 = int(connect[num][1])
                else:
                    elem2 = int(connect[num][0])
                NN1 = lines[3+elem].split()[3]
                NN2 = lines[3+elem2].split()[3]
                num_of_neighbours1 = 0
                num_of_neighbours2 = 0
                communic1 = []
                communic2 = []
                for i in range(end-beg):
                    if int(list(lines[beg+i].split())[0]) == elem:
                        num_of_neighbours1 += 1
                        communic1.append(list(lines[beg+i].split())[0:3])
                    if int(list(lines[beg+i].split())[1]) == elem:
                        num_of_neighbours1 += 1
                        communic1.append(list(lines[beg+i].split())[0:3])
                    if int(list(lines[beg+i].split())[0]) == elem2:
                        num_of_neighbours2 += 1
                        communic2.append(list(lines[beg+i].split())[0:3])
                    if int(list(lines[beg+i].split())[1]) == elem2:
                        num_of_neighbours2 += 1
                        communic2.append(list(lines[beg+i].split())[0:3])
                num_of_neighbours1 = str(num_of_neighbours1)
                type_of_connection1 = conn_type(communic1)
                num_of_neighbours2 = str(num_of_neighbours2)
                type_of_connection2 = conn_type(communic2)
                name1 = NN1 + num_of_neighbours1 + type_of_connection1
                name2 = NN2 + num_of_neighbours2 + type_of_connection2
                if name1 < name2:
                    name = name1 + name2
                else:
                    name = name2 + name1
                
                if file in molekula.index:
                    if name in molekula.columns:
                        if np.isnan(molekula[name][file]):
                            molekula.loc[file, name] = 1
                        else:
                            molekula.loc[file, name] += 1
                    else:
                        molekula.loc[file, name] = 1
                else:
                    molekula.loc[file, name] = 1


molekula = pd.DataFrame({})
# files = ["adrenalin.mol", "soil-1.mol"]
files = ["adrenalin.mol", "amaron.mol", "iervin.mol", 'acetona_peroxid.mol', 'aleiritic_acid.mol', "soil-1.mol"]
for file in files:
    analize_mol_file1(molekula, file)
molekula = molekula.fillna(0)
for file in files:
    analize_mol_file2(molekula, file)
molekula = molekula.fillna(0)
molekula
