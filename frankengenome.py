
with open('frankengenome_train.txt', 'r') as file:
    training_data = file.read().strip()

with open('frankengenome.fa', 'r') as file:
    genome = ''.join(file.read().splitlines()[1:])

def GetTransmission():
    total_trans = len(training_data) - 1
    trans_p = {'0': {'0': 0, '1': 0}, '1': {'0': 0, '1': 0}}
    for i in range(len(training_data) - 1):
        trans_p[training_data[i]][training_data[i+1]] += 1
        
    for row in trans_p:
        for col in trans_p[row]:
            trans_p[row][col] /= total_trans

    return trans_p

def GetEmission():
    emit_p = {'0': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
              '1': {'A': 0, 'C': 0, 'G': 0, 'T': 0}
             }
    total_e = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    for i in range(0, 50000):
        emit_p[training_data[i]][genome[i]] += 1
        total_e[genome[i]] += 1
        
    for state in emit_p:
        for e in total_e:
            emit_p[state][e] /= total_e[e]

    return emit_p

def main():
    trans_p = GetTransmission()
    emit_p = GetEmission()

main()