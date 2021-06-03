import HMM
## YOUR CODE HERE
## you can borrow code from any HMM implementation, but you need to acknowledge the source
## no hmmlearn or similar ML libraries allows

# trans_p = {'0': {'0': 0, '1': 0}, '1': {'0': 0, '1': 0}}
# emits_p = {'0': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
#             '1': {'A': 0, 'C': 0, 'G': 0, 'T': 0}
#             }

# trans_p = {'00': 0, '01': 0, '10': 0, '11': 0}
# emits_p = {'0A': 0, '0C': 0, '0G': 0, '0T': 0, '1A': 0, '1C': 0, '1G': 0, '1T': 0}

with open('frankengenome_train.txt', 'r') as file:
    training_data = file.read().strip()

with open('frankengenome_test.txt', 'r') as file:
    testing_data = file.read().strip()

with open('frankengenome.fa', 'r') as file:
    genome = ''.join(file.read().splitlines()[1:])

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
    trans_p = itp
    total_trans = {'0': 0, '1': 0}
    for i in range(len(pi) - 1):
        trans_p[pi[i]+pi[i+1]] += 1
        total_trans[pi[i]] += 1

    # for row in trans_p:
    #     for col in trans_p[row]:
    #         if total_trans[row] != 0:
    #             trans_p[row][col] /= total_trans[row]
    #         else:
    #             trans_p[row][col] = 1/len(states)
    for item in trans_p:
        if total_trans[item[0]] != 0:
            trans_p[item] /= total_trans[item[0]]
        else:
            trans_p[item] = 1/len(states)

    return trans_p

def GetEmission(string, states, pi, alphabet, iep):
    emits_p = iep
    total_e = {'0': 0, '1': 0}
    
    for i in range(len(pi)):
        emits_p[pi[i]+string[i]] += 1
        total_e[pi[i]] += 1

    # for state in emits_p:
    #     for a in alphabet:
    #         if emits_p[state][a] != 0:
    #             emits_p[state][a] /= total_e[state]
    #         else:
    #             emits_p[state][a] = 1/len(alphabet)
    for item in emits_p:
        if total_e[item[0]] != 0:
            emits_p[item] /= total_e[item[0]]
        else:
            emits_p[item] = 1/len(alphabet)

    return emits_p


def main():
    string = genome[:50000]
    alphabet = ['A', 'C', 'G', 'T']
    pi = training_data
    states = ['0', '1']
    
    init_p = {'0': .5, '1': .5}
    itp, iep = MakeMatrix2(states, pi, string, alphabet)
    # print(itp)
    # print(iep)
    trans_p = GetTransmission(pi, states, itp)
    # print(trans_p)
    emits_p = GetEmission(string, states, pi, alphabet, iep)
    # print(emits_p)
    hmm = HMM.HMM(trans_p, emits_p, init_p)
# #     print(hmm.Alog)
    hidden_path = hmm.viterbiL(genome[50000:100000])
    
    result = 0
    # print(len(hidden_path[1]))
    for i in range(len(hidden_path[1])):
        if hidden_path[1][i] == testing_data[i]:
            result += 1
    accuracy = round(result/len(testing_data), 4)
    print('Accuracy on frankengenome_test:', accuracy * 100)

main()