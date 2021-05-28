
with open('frankengenome_train.txt', 'r') as file:
    training_data = file.read().strip()

def GetTransmission():
    total_trans = len(training_data) - 1
    trans_p = {'0': {'0': 0, '1': 0}, '1': {'0': 0, '1': 0}}
    for i in range(len(training_data) - 1):
        trans_p[training_data[i]][training_data[i+1]] += 1
        
    for row in trans_p:
        for col in trans_p[row]:
            trans_p[row][col] /= total_trans

    return trans_p

def main():
    trans_p = GetTransmission()

main()