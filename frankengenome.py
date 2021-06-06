## YOUR CODE HERE
## you can borrow code from any HMM implementation, but you need to acknowledge the source
## no hmmlearn or similar ML libraries allows

## ======================== CREDITS ======================== ##
##I borrowed Ben Langmead's code for the HMM class, which is linked in the project description. The Viterbi implementation from his code was used to determine the hidden path for the genome
## because we were given the emitted string and the hidden path for the first 50k labels

import HMM

with open('frankengenome_train.txt', 'r') as file:
    training_data = file.read().strip()

with open('frankengenome_test.txt', 'r') as file:
    testing_data = file.read().strip()

with open('frankengenome.fa', 'r') as file:
    genome = ''.join(file.read().splitlines()[1:])

def CalcAccuracy(hidden_path):
    result = 0
    for i in range(len(hidden_path[1])):
        if hidden_path[1][i] == testing_data[i]:
            result += 1
    accuracy = round(result/len(testing_data), 4)
    return accuracy

def MakeMatrix2(states, pi, string, alphabet):
    trans_p = {}
    emits_p = {}
    # Make trans_p
    for state1 in states:
        # trans_p[state1] = {}
        for state2 in states:
            trans_p[state1+state2] = 0

    # Make emits_p
    for state in states:
        # emits_p[state] = {}
        for letter in alphabet:
            emits_p[state+letter] = 0
    
    return trans_p, emits_p

def GetTransmission(pi, states, itp):
    # trans_p = {'A': {'A': 0, 'B': 0, 'C': 0}, 'B': {'A': 0, 'B': 0, 'C': 0}, 'C': {'A': 0, 'B': 0, 'C': 0}}
    trans_p = itp.copy()
    total_trans = {'0': 0, '1': 0}
    for i in range(len(pi) - 1):
        trans_p[pi[i]+pi[i+1]] += 1
        total_trans[pi[i]] += 1

    for item in trans_p:
        if total_trans[item[0]] != 0:
            trans_p[item] /= total_trans[item[0]]
        else:
            trans_p[item] = 1/len(states)

    return trans_p

def GetEmission(string, pi, alphabet, iep):
    emits_p = iep.copy()
    total_e = {'0': 0, '1': 0}
    
    for i in range(len(pi)):
        emits_p[pi[i]+string[i]] += 1
        total_e[pi[i]] += 1

    for item in emits_p:
        if total_e[item[0]] != 0:
            emits_p[item] /= total_e[item[0]]
        else:
            emits_p[item] = 1/len(alphabet)

    return emits_p

def PlotData(data_size, accuracies):
    HMM.plt.axis([0, 50000, 0, 100])
    HMM.plt.plot(data_size, accuracies)

def main():
    # Preparation
    string = genome[:50000] # Training string which is the first 50k labels of frankengenome
    alphabet = ['A', 'C', 'G', 'T']
    pi = training_data # Hidden path from frankengenome_train
    states = ['0', '1']
    init_p = {'0': .5, '1': .5}
    itp, iep = MakeMatrix2(states, pi, string, alphabet)
    
    data_size = []
    accuracies = []
    
    for i in range(9):
        if i == 0:
            data_size.append(50000)
            print('Calculating Paremeters with 50,000 labels...')
            trans_p = GetTransmission(pi, states, itp)
            emits_p = GetEmission(string, pi, alphabet, iep)
        elif i == 1:
            data_size.append(25000)
            temp_pi = pi[:25000]
            temp_string = string[:25000]
            print('Calculating Parameters with 25,000 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 2:
            data_size.append(10000)
            temp_pi = pi[:10000]
            temp_string = string[:10000]
            print('Calculating Parameters with 10,000 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 3:
            data_size.append(5000)
            temp_pi = pi[:5000]
            temp_string = string[:5000]
            print('Calculating Parameters with 5,000 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 4:
            data_size.append(1000)
            temp_pi = pi[:1000]
            temp_string = string[:1000]
            print('Calculating Parameters with 1,000 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 5:
            data_size.append(500)
            temp_pi = pi[:500]
            temp_string = string[:500]
            print('Calculating Parameters with 500 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 6:
            data_size.append(100)
            temp_pi = pi[:100]
            temp_string = string[:100]
            print('Calculating Parameters with 100 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 7:
            data_size.append(10)
            temp_pi = pi[:10]
            temp_string = string[:10]
            print('Calculating Parameters with 10 labels...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
        elif i == 8:
            data_size.append(1)
            temp_pi = pi[:1]
            temp_string = string[:1]
            print('Calculating Parameters with 1 label...')
            trans_p = GetTransmission(temp_pi, states, itp)
            emits_p = GetEmission(temp_string,temp_pi, alphabet, iep)
            print('Transition probability matrix:', trans_p)
            
        # Run HMM on data
        hmm = HMM.HMM(trans_p, emits_p, init_p)
        print('Running ViterbiL on frankengenome_test...')
        hidden_path = hmm.viterbiL(genome[50000:100000])
        
        accuracy = round(CalcAccuracy(hidden_path) * 100, 4)
        accuracies.append(accuracy)
        print(f'Accuracy on frankengenome_test: {accuracy}%')
        
        cont = input('Continue?')
        
    PlotData(data_size, accuracies)

main()