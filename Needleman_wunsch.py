import numpy as np

#sequence_1 = input("Enter sequence 1:")
#sequence_2 = input("Enter sequence 2:")

sequence_2 = "CGTGGACTTGCTTACGTACCGTGCCGATGCC"
sequence_1 = "GACTGGCTAGTTCAGTCTGACTGGCTAGCTA"

main_matrix = np.zeros((len(sequence_1)+1,len(sequence_2)+1))
match_checker_matrix = np.zeros((len(sequence_1),len(sequence_2)))

# Defining the scoring matrix
match_reward = 1
mismatch_penalty = -1
gap_penalty = -2

# Initialising the match_checker_matrix which keeps tracks of all the matches
for i in range(len(sequence_1)):
    for j in range(len(sequence_2)):
        if sequence_1[i] == sequence_2[j]:
            match_checker_matrix[i][j]= match_reward
        else:
            match_checker_matrix[i][j]= mismatch_penalty

#print(match_checker_matrix)

#STEP 1 : Initialisation

# Initilising the main matrix with gap penalties
for i in range(len(sequence_1)+1):
    main_matrix[i][0] = i*gap_penalty
for j in range(len(sequence_2)+1):
    main_matrix[0][j] = j * gap_penalty

#STEP 2 : Matrix Filling

# Matrix filling by comparing all three possible cases
for i in range(1,len(sequence_1)+1):
    for j in range(1,len(sequence_2)+1):
        main_matrix[i][j] = max(main_matrix[i-1][j-1]+match_checker_matrix[i-1][j-1],
                                main_matrix[i-1][j]+gap_penalty,
                                main_matrix[i][j-1]+ gap_penalty)

print(main_matrix, '\n')

# STEP 3 : Traceback
aligned_1 = ""
aligned_2 = ""

ti = len(sequence_1)
tj = len(sequence_2)

while(ti > 0 and tj > 0):

    # Comparing for diagonal match or mismatch
    if (ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj-1]+ match_checker_matrix[ti-1][tj-1]):

        aligned_1 = sequence_1[ti-1] + aligned_1
        aligned_2 = sequence_2[tj-1] + aligned_2

        ti = ti - 1
        tj = tj - 1

    # Comparing for gap from upper side
    elif(ti > 0 and main_matrix[ti][tj] == main_matrix[ti-1][tj] + gap_penalty):
        aligned_1 = sequence_1[ti-1] + aligned_1
        aligned_2 = "-" + aligned_2

        ti = ti -1

    # Comparing for gap from left side
    else:
        aligned_1 = "-" + aligned_1
        aligned_2 = sequence_2[tj-1] + aligned_2

        tj = tj - 1

print(aligned_1)
print(aligned_2)
