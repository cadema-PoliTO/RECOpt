## Levenshtein distance

import numpy as np

####################################################################################################################################################################################

# This module contains a method that uses the Levensthein distance between two string, in order
# to suggest a correction to a word provided from the user (test_word), drawn from a list of possible 
# words (words)

def Leven_dist_comparison(test_word, words):

    ## Levensthein algorithm
    # Given two words (test_word and word), being n_rows and n_cols, repsectively, the lengths of the two
    # strings increased by one, a matrix is created accounting for the distance between each two possible 
    # pair of substrings of the first two strings (including a null string). For each pair of substrings
    # (:j_test)|(:j_word) such distance is given by the minimum value between d((:j_test - 1),(:j_word)) + 1,
    # that accounts for the addition of a character to the first substring, d((:j_test1),(:j_word - 1)) + 1,
    # that accounts for the removal of a character from the second substring and d((:j_test1 - 1),(:j_word - 1)) + increase,
    # that accounts for the change of a character in the substring, where increase can either be 0 if the current characters
    # are equal, or 1 otherwise.

    # Creating a vector where to store the distances between the test word and the other ones (count runs through the vector)
    distances = np.zeros((len(words)))
    count = 0

    # Setting the number of rows (axis = 0) of the Levenshtein matrix 
    n_rows = len(test_word) + 1

    # Running through the list of words, in order to evaluate their distance from the test word
    for word in words:

        # Setting the number of columns (axis = 1) of the Levenhtein matrix
        n_cols = len(word) + 1

        # Initializing the Levensthein matrix
        leven_matrix = np.zeros((n_rows, n_cols))

        # Setting the first row and columns of the Levensthein matrix (distance of each sub-string from the null string)
        leven_matrix[:, 0] = np.arange(n_rows)
        leven_matrix[0, :] = np.arange(n_cols)
    
        # Running thorught the rows and the columns in order to evaluate the distance between each sub-string
        for i_rows in range(1, n_rows):

            # Index to be used for slicing the test_word string
            j_test = i_rows - 1

            for i_cols in range(1, n_cols):

                # Index to be used for slicing the word string
                j_word = i_cols - 1

                # Increase in the distance between the substrings (: jtest)|(:jword) and (: jtest - 1)|(:jword - 1),
                # that is equal to 0 if the current characters are equal, 1 otherwise
                increase = 0
                if test_word[j_test] != word[j_word]: increase = 1

                # Evaluating the distance between the substring (: jtest)|(:jword)
                leven_matrix[i_rows, i_cols] = min(leven_matrix[i_rows - 1, i_cols] + 1,
                                                   leven_matrix[i_rows, i_cols - 1] + 1,
                                                   leven_matrix[i_rows- 1, i_cols - 1] + increase)
        
        # Storing the distance evaluated and updating the counter
        distances[count] = leven_matrix[n_rows - 1, n_cols - 1]
        count += 1


    ## Post-processing

    # The closest word is evaluated as the word with the minimum distance from the test_word.
    # Since many minima may appear, they must be all found
    indices = np.nonzero(distances == np.min(distances))[0]

    
    return([words[index] for index in indices])
