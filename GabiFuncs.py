import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import operator as op
from functools import reduce


def ncr(n, r):
    # calculates {n \choose r}
    r = min(r, n - r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer / denom


def flip_probability(k, r, p):
    '''
    returns the probability of having at most r flips in k bits, where each flips occurs with probability p
    '''
    return np.sum([ncr(k, j) * p ** j * (1 - p) ** (k - j) for j in range(r + 1)])


# flip_probability(20,2,0.15)

def find_window_for_flip_probability(delta, r, p, epsilon=0.001, k_max=1000):
    '''
    returns the largest window size k for which the probability of having at most r flips in the window is at least delta
    '''
    k_min = r
    k = k_max
    while flip_probability(k, r, p) < delta + epsilon:
        print 'k_min=%d, k_max=%d' % (k_min, k_max)
        if flip_probability((k_min + k_max) / 2, r, p) > delta:
            k_min = (k_min + k_max) / 2
        else:
            k_max = (k_min + k_max) / 2
        if k_max - k_min < 2:
            k = k_min
            break
    return k


# find_window_for_flip_probability(0.5, 2, 0.15)

def expected_number_flips(k, p):
    '''
    return the expected number of flips in a window of length k
    '''
    return np.sum([j * ncr(k, j) * p ** j * (1 - p) ** (k - j) for j in range(k + 1)])


# expected_number_flips(17,0.15)

def hamming_dist(s1, s2):
    '''
    Calculate the Hamming distance between two bit strings
    '''
    if len(s1) == len(s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    else:
        return -1  # error return value

def iterative_levenshtein(s, t, costs=(1, 1, 1), print_match_length=0):
    """
        SLOW!!!!
        iterative_levenshtein(s, t) -> ldist
        ldist is the Levenshtein distance between the strings
        s and t.
        For all i and j, dist[i,j] will contain the Levenshtein
        distance between the first i characters of s and the
        first j characters of t

        costs: a tuple or a list with three integers (d, i, s)
               where d defines the costs for a deletion
                     i defines the costs for an insertion and
                     s defines the costs for a substitution
        print_match_length: block length for printing matche description.
                            if 0, then no print
        the conclusion is for s=source and t =res (important for the relation between s and t, deletion and insertions)
    """
    rows = len(s) + 1
    cols = len(t) + 1
    deletes, inserts, substitutes = costs
    conclusion = {"F": 0, "D": 0, "I": 0, "DIST":0} # F for flips, D for del, I for insertion
    pointer_dict = {}

    dist = np.array([[0 for x in range(cols)] for x in range(rows)] )#zeros matrix

    # source prefixes can be transformed into empty strings
    # by deletions:
    for row in range(1, rows): dist[row][0] = row * deletes

    # target prefixes can be created from an empty source string
    # by inserting the characters
    for col in range(1, cols): dist[0][col] = col * inserts

    for col in range(1, cols):
        for row in range(1, rows):
            if s[row - 1] == t[col - 1]:
                cost = 0
            else:
                cost = substitutes
            dist[row][col] = min([dist[row - 1][col] + deletes,
                                  dist[row][col - 1] + inserts,
                                  dist[row - 1][col - 1] + cost])  # substitution
            if print_match_length > 0: #the path throw the matrix, step by step
                if dist[row][col] == dist[row - 1][col - 1] + cost:
                    pointer_dict[(row, col)] = (row - 1, col - 1)
                elif dist[row][col] == dist[row][col - 1] + inserts:
                    pointer_dict[(row, col)] = (row, col - 1)
                else:
                    pointer_dict[(row, col)] = (row - 1, col)
    #    for r in range(rows):
    #        print(dist[r])
    print pointer_dict
    if print_match_length > 0:
        new_s = ''
        new_t = ''
        match = ''
        p_row = row
        p_col = col
        while (p_row != 0) and (p_col != 0):
            #            print (p_row,p_col)
            if (p_row == pointer_dict[(p_row, p_col)][0] + 1) and (p_col == pointer_dict[(p_row, p_col)][1] + 1):
                new_s += s[p_row - 1]
                new_t += t[p_col - 1]
                if new_s[-1] == new_t[-1]:
                    match += '*'
                else:
                    match += 'f'
                    conclusion['F']+=1
            elif (p_row == pointer_dict[(p_row, p_col)][0] + 1):  # only the column moves back
                new_s += s[p_row - 1]
                new_t += '-'
                match += '-'
                conclusion['D'] += 1
            else:
                new_s += '-'
                new_t += t[p_col - 1]
                match += '-'
                conclusion['I'] += 1
            p_row, p_col = pointer_dict[(p_row, p_col)]
        if (p_row == 0) and (p_col > 0):
            new_s += '-' * p_col
            new_t += t[:p_col]
            match += '-' * (p_col)
            conclusion['I'] += 1*p_col
        elif (p_row > 0) and (p_col == 0):
            new_s += s[:p_row]
            new_t += '-' * (p_row)
            match += '-' * (p_row)
            conclusion['D'] += 1*(p_row)
        # reverse the strings
        new_s = new_s[::-1]
        new_t = new_t[::-1]
        match = match[::-1]
        block = 0
        while (block + 1) * print_match_length < len(new_s):
            print new_s[block * print_match_length:(block + 1) * print_match_length]
            print new_t[block * print_match_length:(block + 1) * print_match_length]
            print match[block * print_match_length:(block + 1) * print_match_length] + '\n'
            block += 1
        if len(new_s) > block * print_match_length:
            print new_s[block * print_match_length:]
            print new_t[block * print_match_length:]
            print match[block * print_match_length:] + '\n'
    conclusion['DIST']=dist[row][col]
    return conclusion


def levenshtein_edit_dist(s1, s2, show_strings=False):
    '''
    FASTER THAN iterative_levenshtein()
    Calculate the Levenshtein edit distance between two bit strings
    also returns the s1 locations where an exact match was found
    '''
    #    if len(s1) > len(s2):
    #        s1, s2 = s2, s1
    #    s1 = '1001101'
    #    s2 = '001011'
    conclusion={'F':0, 'D':0, 'I':0,'DIST':0}
    pointers = {}
    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2 + 1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
                pointers[(i2, i1)] = (i2 - 1, i1 - 1)
            else:
                candidates = (distances[i1], distances[i1 + 1], distances_[-1])
                candidates_pointers = ((i2 - 1, i1 - 1), (i2 - 1, i1), (i2, i1 - 1))
                min_ind = np.argmin(candidates)
                distances_.append(1 + candidates[min_ind])
                pointers[(i2, i1)] = candidates_pointers[min_ind]
        distances = distances_


    # 9-3-18
    # if distances[-1] > len(s1)/8:
    #     dict = ({"DIST": 1, "I": 2, "D": 3, "F": 4},)
    #     dict[0]["DIST"]=len(s1)/8
    #     dict[0]["I"] = len(s1) / 8
    #     dict[0]["F"] = len(s1) / 8
    #     dict[0]["D"] = len(s1) / 8
    #     return dict

    s1_match_indices = []
    curr_pointer = (len(s2) - 1, len(s1) - 1)
    while (curr_pointer[0] > -1) and (curr_pointer[1] != -1):
        if (curr_pointer[0] == pointers[curr_pointer][0] + 1) and (curr_pointer[1] == pointers[curr_pointer][1] + 1):
            if s2[curr_pointer[0]] == s1[curr_pointer[1]]:
                s1_match_indices.append(curr_pointer[1])
            else: conclusion['F']+=1
        elif (curr_pointer[0] == pointers[curr_pointer][0]) and (curr_pointer[1] == pointers[curr_pointer][1] + 1):
            conclusion['D'] += 1
        else: conclusion['I'] += 1
        curr_pointer = pointers[curr_pointer]
    if(curr_pointer[0] == -1) and (curr_pointer[1] != -1):
        conclusion['D'] += 1*(curr_pointer[1]+1)
    if (curr_pointer[0] > -1) and (curr_pointer[1] == -1):
        conclusion['I'] += 1 * (curr_pointer[0] + 1)
    conclusion['DIST']=distances[-1]
    return conclusion, s1_match_indices

def bitwise_majority_string(string_array):
    '''
    returns the majority string corresponding to an array of strings of the same length
    '''
    return ''.join(np.array(
        np.array([np.mean(x) for x in np.array([list(x) for x in string_array]).transpose().astype(int)]) > 0.5).astype(
        int).astype(str))


def stringify(array):
    '''
    turns a binary numpy vector into a string of 0/1
    '''
    return ''.join(array.astype(str))


def intify(bin_string):
    '''
    turns a string of 0/1 into a binary numpy vector
    '''
    return np.array(list(bin_string))


def init_key(n, rand_seed=-1):
    '''
    returns a random key of n bits
    '''
    print 'initiating key...'
    if rand_seed > 0:
        np.random.seed(rand_seed)
    print 'DONE!'
    return ''.join((np.random.rand(n) > 0.5).astype(int).astype(str))


def build_samples(key, num_samples, sample_len, window_size, flip_probability, delete_probability, insert_probability):
    '''
    build snippets dataset, where each sample is noisified, and then sliced using a sliding window into snippets
    '''
    print 'building samples...'
    result_dict = {}
    for sample_idx in range(num_samples):
        if sample_idx % 1000 == 0:
            print sample_idx
        sample_start = np.random.randint(n - sample_len + 1)
        sample = np.array(list(key[sample_start:sample_start + sample_len])).astype(int)
        # flip random bits
        rand_flip = (np.random.rand(len(sample)) < flip_probability).astype(int)
        sample[np.where(rand_flip > 0)] = 1 - sample[np.where(rand_flip > 0)]
        # delete random bits
        rand_delete = (np.random.rand(len(sample)) < delete_probability).astype(int)
        sample = sample[np.where(rand_delete == 0)]
        # insert random bits
        rand_insert = (np.random.rand(len(sample) + 1) < insert_probability).astype(int)
        rand_bits_to_insert = (np.random.rand(sum(rand_insert)) > 0.5).astype(int)
        sample = np.insert(sample, np.where(rand_insert == 1)[0], rand_bits_to_insert)
        # scan windows of sample
        for window_start in range(len(sample) - window_size + 1):
            window = sample[window_start:window_start + window_size]
            window_key = ''.join(window.astype(str))
            if window_key not in result_dict:
                result_dict[window_key] = {'sample': window_key,
                                           'count': 1,
                                           'weight': sum(window),
                                           'similar_count': 0,
                                           'sample_start': [sample_start + window_start],
                                           'similar_samples': [],
                                           'closest_majority_sample': ''}
            else:
                result_dict[window_key]['count'] += 1
                if sample_start + window_start not in result_dict[window_key]['similar_samples']:
                    result_dict[window_key]['similar_samples'] = result_dict[window_key]['similar_samples'] + [
                        sample_start + window_start]
    result_df = pd.DataFrame.from_dict(result_dict, orient='index').sort_values(by='weight')
    print 'DONE!'
    return result_df


def prune_samples(result_df, min_count=-1):
    '''
    returns a subset of the snippets dataset which consists only of snippets that show high statistical significance
    '''
    print 'prunning samples...'
    if min_count < 0:
        min_count = result_df['count'].quantile(.5)
    common_samples_df = result_df[result_df['count'] > min_count]
    #    noisy_samples_df = result_df[result_df['count'] < 3].sort_values(by='weight')
    #    for idx, common_sample in enumerate(common_samples_df.index):
    #        print idx, common_sample
    #        candidate_noisy_samples_above_weight_df = noisy_samples_df[noisy_samples_df['weight'] > common_samples_df.loc[common_sample]['weight'] -2]
    #        candidate_noisy_samples_below_weight_df = candidate_noisy_samples_above_weight_df [candidate_noisy_samples_above_weight_df['weight'] < common_samples_df.loc[common_sample]['weight'] +2]
    #        candidate_noisy_samples_df = candidate_noisy_samples_below_weight_df[candidate_noisy_samples_below_weight_df['weight'] != common_samples_df.loc[common_sample]['weight']]
    #        for noisy_sample in candidate_noisy_samples_df.index:
    #            if hamming_dist(common_sample,noisy_sample) < 2:
    #                common_samples_df.set_value(common_sample, 'similar_count', common_samples_df.loc[common_sample]['similar_count'] + 1)
    print 'DONE!'
    return common_samples_df.sort_values(by='count', ascending=False)  # more reliable samples first


def build_shift_pointers(common_samples_array, stitch_shift_size):
    '''
    build DAG where snippets are connected if they can be stitched by a small shift
    only the highest-ranking snippet that can be stitched is used, where the snippets array is assumed to be sorted by popularity
    '''
    print 'building DAG...'
    shift_pointers = {'right_index': {}, 'left_index': {}}
    for idx1, left_sample in enumerate(common_samples_array):
        if idx1 % 100 == 0:
            print idx1
        for stitch_shift in range(1, stitch_shift_size + 1):
            #            print 'stitch_shift %d' % stitch_shift
            for idx2, right_sample in enumerate(common_samples_array):
                #                print idx2, sample2
                if hamming_dist(left_sample[stitch_shift:], right_sample[:-stitch_shift]) == 0:
                    if right_sample not in shift_pointers['right_index']:
                        shift_pointers['right_index'][right_sample] = {'right_sample': right_sample,
                                                                       'left_sample': left_sample,
                                                                       'shift': stitch_shift}
                    if left_sample not in shift_pointers['left_index']:
                        shift_pointers['left_index'][left_sample] = {'right_sample': right_sample,
                                                                     'left_sample': left_sample, 'shift': stitch_shift}
                    break_stitch_shift_loop = True
                    break
            if break_stitch_shift_loop:
                break
    print 'DONE!'
    return shift_pointers


def stitch(common_samples_array, shift_pointers):
    '''
    traverse the DAG, starting from the sinks, and generate as long sequences as possible
    the algorithm assumes each snippet (node) has at most one incoming link
    if it should support multiple incoming links, then the function should be adjusted
    '''
    #    shift_pointers_right_index_df = pd.DataFrame(shift_pointers['right_index']).transpose()
    #    shift_pointers_left_index_df = pd.DataFrame(shift_pointers['left_index']).transpose()
    start_samples = []
    for sample in common_samples_array:
        if sample not in shift_pointers['right_index']:
            start_samples += [sample]
    retrieved_key = []
    for start_sample in start_samples:
        print 'START SAMPLE: ' + start_sample
        curr_sample = start_sample
        curr_key = curr_sample
        total_shift = 0
        while curr_sample in shift_pointers['left_index']:
            #            print ' '*total_shift + curr_sample
            curr_key += shift_pointers['left_index'][curr_sample]['right_sample'][
                        -shift_pointers['left_index'][curr_sample]['shift']:]
            total_shift += shift_pointers['left_index'][curr_sample]['shift']
            curr_sample = shift_pointers['left_index'][curr_sample]['right_sample']
        print curr_key
        retrieved_key += [curr_key]
    return retrieved_key


def validate_sample(sample, near_sample_df, radius):
    '''
    sample is a Series representing a single sample in the original dataframe
    near_sample_df is a dataframe of samples that have a weight similar to that of sample
    radius is the threshold criteria for considering a near_sample as a noisy version of sample

    returns the number of near samples that satisfy the threshold for sufficiently many bits of sample
    '''
    sample_count = np.ones(len(sample['sample'])) * sample['count']
    for near_sample in near_sample_df.index:
        (dist_arr, match_array) = levenshtein_edit_dist(sample['sample'], near_sample)
        dist=dist_arr['DIST']
        if len(match_array) >= len(sample['sample']) - radius:
            sample_count[match_array] += 1
    #            print sample['sample']
    #            print near_sample, len(match_array), '\n'
    return min(sample_count)






def prune_samples_extended(result_df, min_count=-1, ignore_similar=True, min_count_radius=1, levenshtein_radius=2):
    '''
    Consider common_samples_df: contain samples that are common, but not very common (have count > min_count).
    Consider noisy_samples_df: contain samples that are even more rare, but not complete outliers (have count > min_count - min_count_radius).
    For each common_sample, look for noisy samples that can be near this common_sample. Such samples must have weight close to that of common_sample.
    Add the noisy samples that are close to common_sample to the similar_count of common_sample.

    returns a subset of the snippets dataset which consists only of snippets that show high statistical significance,
    when considering their similar_count, i.e., their actual count as well as their near noisy samples
    '''

    #    min_count=16
    #    ignore_similar=False
    #    min_count_radius=2
    #    levenshtein_radius=2

    print 'prunning samples...'
    if min_count < 0:
        min_count = result_df['count'].quantile(.5)
    if ignore_similar:
        common_samples_df = result_df[result_df['count'] > min_count]
        print 'DONE!'
        return common_samples_df.sort_values(by='count', ascending=False)  # more reliable samples first
    else:
        weight_radius = levenshtein_radius  # *2 # for slicing potential close samples
        common_samples_df = result_df[result_df['count'] > min_count].sort_values(by='count', ascending=False)
        noisy_samples_df = result_df[result_df['count'] > min_count - min_count_radius].sort_values(by='weight')
        for idx, common_sample in enumerate(common_samples_df.index):
            if idx % 10 == 0:
                print idx
            near_weight_sample_df = noisy_samples_df[noisy_samples_df['weight'].isin(
                np.arange(common_samples_df.iloc[idx]['weight'] - weight_radius,
                          common_samples_df.iloc[idx]['weight'] + weight_radius))]
            #    common_samples_df['similar_count'] = common_samples_df.apply(lambda row: validate_sample(row, near_weight_sample, levenshtein_radius))
            #    print idx, common_sample
            #    common_sample = common_samples_df.head(2).tail(1)
            #    candidate_noisy_samples_in_radius_df = noisy_samples_df[noisy_samples_df['weight'].isin(
            #                    np.arange(common_sample['weight'] - weight_radius , common_sample['weight'] + weight_radius)) ]
            similar_count = validate_sample(common_samples_df.iloc[idx], near_weight_sample_df, levenshtein_radius)
            common_samples_df.at[common_sample, 'similar_count'] = similar_count
        print 'DONE!'
        common_samples_df.sort_values(by='similar_count', ascending=False)['count'].hist(bins=100)
        return common_samples_df.sort_values(by='similar_count', ascending=False)  # more reliable samples first


