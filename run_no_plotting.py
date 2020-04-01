#!/usr/bin/python
from __future__ import division
import numpy as np
import pandas as pd
import random
import csv
import sys
import os

def mq(folder, protein_list, filter_con = True):
    """ Reads from evidence and peptide MQ results files and returns a
    dictionaryof the form:
    total_data[sample][protein] = [intensity, pre, start, seq,
                                   end, post, sample, mods]
    """
    print "Reading file from MQ"

    # Evidence and peptide folders
    evi_file = os.path.join(folder, "evidence.txt")
    pep_file = os.path.join(folder, "peptides.txt")

    # To store evidence and peptide data
    e_data = []
    p_data = {}

    # Read relevant columns from evidence file
    with open(evi_file) as csvfile:
        e_reader = csv.reader(csvfile, delimiter='\t')
        e_headers = next(e_reader)
        e_protein = e_headers.index("Leading razor protein")
        e_seq = e_headers.index("Sequence")
        e_i = e_headers.index("Intensity")
        e_pep_id = e_headers.index("Peptide ID")
        e_mod_seq = e_headers.index("Modified sequence")
        # Get experiment, if this was run with only 1, make a default value
        try:
            e_exp = e_headers.index("Experiment")
        except ValueError:
            e_exp = "Default"

        for row in e_reader:
            if e_exp != "Default":
                experiment = row[e_exp]
            else: experiment = e_exp

            protein = row[e_protein]
            # Filter out contaminants if active
            if (filter_con and "CON_" not in protein) or (not filter_con):
                if row[e_protein] in protein_list or len(protein_list) == 0:

                    # Experiment, seq, mods, intensity, ID, protein
                    e_data.append([experiment, row[e_seq], row[e_mod_seq], row[e_i],
                               row[e_pep_id], row[e_protein]])

    # Read positional info from peptide file
    with open(pep_file) as csvfile2:
        p_reader = csv.reader(csvfile2, delimiter='\t')
        p_headers = next(p_reader)
        p_pre = p_headers.index("Amino acid before")
        p_post  = p_headers.index("Amino acid after")
        p_start = p_headers.index("Start position")
        p_end = p_headers.index("End position")
        p_ID = p_headers.index("id")

        for row in p_reader:
            # ID : Pre, post, start, end
            p_data[row[p_ID]] = [row[p_pre], row[p_post],
                                 row[p_start], row[p_end]]

    # Combine p_data and e_data to total_data
    total_data = {}
    for row in e_data:
        sample, seq, mods, intensity, ID, protein = row
        if sample not in total_data:
            total_data[sample] = {}
        # No REVERSE hits
        if intensity != "" and "REV" not in protein:
            if ID in p_data:
                pre, post, start, end = p_data[ID]
                if pre != "" and start != "":
                    mods_list = []
                    possible_mods = ["de"] #TODO: Add more mods
                    # If mods found, add to total_data in a readable format
                    for m in possible_mods:
                        if m in mods:
                            actual_pos = 1 # For counting position in sequence
                            for i in xrange(1, len(mods)):
                                if mods[i] == "(":
                                    if mods[i+1] == m[0] and mods[i+2] == m[1]:
                                        mod_label = "%s-%i-%s" % (mods[i-1],
                                                     (int(start)+actual_pos-2),
                                                     m)
                                        mods_list.append(mod_label)
                                # Do not advance unless actual sequence
                                if mods[i].isupper(): actual_pos += 1
                    info = [float(intensity), pre, int(start), seq, int(end),
                            post, sample, mods_list]
                    if protein not in total_data[sample]:
                        total_data[sample][protein] = [info]
                    else:
                        total_data[sample][protein].append(info)
    print "Data loaded"
    return total_data


def get_robinson(aa):
    """ Takes an amino acid, N or Q, and calls read_rate
    Returns a Pandas dataframe
    """
    if aa == "Q": return read_rate("Info/gln.csv")
    if aa == "N": return read_rate("Info/asn.csv")


def read_rate(filename):
    """ Reads the RR rates into a Pandas DataFrame
    If format of RR csvs changes, edit this
    """
    data = []
    row_names, col_names = [], []
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        row = reader.next()
        col_names = row[1:] # X amino acid names
        for line in reader:
            line = filter(None, line)
            row_names.append(line[0]) # Y amino acid names
            data.append(line[1:]) # Half times
    data = np.array(data)
    df = pd.DataFrame(data, index=row_names, columns=col_names)
    return df


def get_mid(total_data):
    """ Forms two MIDs from total data
    mid[sample][protein][label] = [hl, [mod_i, total_i], [mod_i, total_i]...]
    Has a dataset of mod intensity vs total intensity and half time for each
    * mid_ss: unique X {N/Q} Y combo
    * mid_classic: position, amino acid combo
    """
    robs = {}
    robs["N"] = get_robinson("N")
    robs["Q"] = get_robinson("Q")
    mid_classic, mid_ss = {}, {}
    for sample in total_data:
        if sample not in mid_classic: mid_classic[sample] = {}
        if sample not in mid_ss: mid_ss[sample] = {}
        for protein in total_data[sample]:
            if protein not in mid_classic[sample]: mid_classic[sample][protein] = {}
            if protein not in mid_ss[sample]: mid_ss[sample][protein] = {}
            for aa in ["N", "Q"]: # For each deamiting amino acid
                rob = robs[aa]
                for row in total_data[sample][protein]:
                    intensity, pre, start, seq, end, post, sample, mods = row
                    cur_pos = start
                    for aa_i in range(len(seq)): # Iterating through sequence
                        aa_in_seq = seq[aa_i]
                        if aa == aa_in_seq: # If we stumble upon the current deaminating amino acid in the sequence
                            mod_i, total_i = 0, 0
                            # Default half time
                            hf = -1

                            total_i += intensity
                            # If this position was modified, its mod label would be this
                            q_mod = "%s-%i-de" % (aa, cur_pos)
                            # If that label is in the mod list for this pos
                            if q_mod in mods:
                                # Increase mod intensity
                                mod_i += intensity

                            if 0 < aa_i < len(seq)-1: # As long as the AA has two neighbours
                                x = seq[aa_i+1] # Look right
                                y = seq[aa_i-1] # Look left
                                combo = "%s %s %s" % (y, aa, x) # SS label

                                # Get half time from Pandas dataframe
                                if (y in rob.index) and (x in rob.columns):
                                    hf = rob.ix[y, x]

                                # Add to SS MID
                                if combo in mid_ss[sample][protein]:
                                    mid_ss[sample][protein][combo].append([mod_i, total_i])
                                else: mid_ss[sample][protein][combo] = [hf, [mod_i, total_i]]

                            # Add to classic MID
                            label = "%s %s" % (cur_pos, aa) # Classic label
                            if label in mid_classic[sample][protein]:
                                [old_mod_i, old_total_i] = mid_classic[sample][protein][label]
                                mod_i += old_mod_i
                                total_i += old_total_i
                            mid_classic[sample][protein][label] = [mod_i, total_i]
                        # On to the next AA
                        cur_pos += 1
    return mid_classic, mid_ss


def calc_deam(mid, to_print = True):
    """ Calculates bulk deamidation per sample, per protein
    Returns a list of [sample, protein, rel_asn, rel_gln]
    Where rel_asn and _gln are mod_intensity/total_intensity
    """
    relative = []
    print_results = []
    for sample in mid:
        for protein in mid[sample]:
            asn_m, asn_t, gln_m, gln_t = 0, 0, 0, 0
            rel_asn, rel_gln = -1, -1
            for label, val in mid[sample][protein].items():
                mod, total = val
                cur_pos, aa = label.split(" ")
                print_results.append([sample, protein, aa, 1-(mod/total)])
                if "N" == aa: # Asn
                    asn_m += mod
                    asn_t += total
                if "Q" == aa: # Gln
                    gln_m += mod
                    gln_t += total
            # Calculate relative amounts
            if asn_t > 0:
                rel_asn = 1 - (asn_m/asn_t)
            if gln_t > 0:
                rel_gln = 1 - (gln_m/gln_t)
            # Each sample relative amounts
            relative.append([sample, protein, rel_asn, rel_gln])
    if to_print: save_fine_bulk(print_results)
    if to_print: save_csv_results(relative, "Bulk")
    print "Bulk deamidation calculated"
    return relative


def bulk_deam(mid, show = False, to_print = True):
    """ Plots bulk deamidation
    """
    relative = np.array(calc_deam(mid))
    index = np.arange(len(relative[:,0]))
    relative = np.array(sorted(relative, key=lambda row:row[0])) # Sort by sample

    if to_print: save_csv_results(relative, "Bulk")


def ss_wrangle(mid):
    """ Takes a mid and returns a NP friendly format for plotting
    """
    data = []
    for sample in mid:
        for protein in mid[sample]:
            for label, val in mid[sample][protein].items():
                hf = val[0]
                after = val[1:]
                mod = [item[0] for item in after]
                total = [item[1] for item in after]
                info = [sample, protein, label, hf, (np.mean(mod)/np.mean(total)), np.mean(total)]
                data.append(info)
    data = np.array(data)
    return data


def get_relative_size(ti, mid, sample, protein, single_sample):
    """ Takes a total intensity of a point combo and returns it as a size
    relative to all total intensities in that sample
    (or, if single sample, relative to that protein)
    """
    # Get total intensities
    totals = []
    if single_sample:
        for label in mid[sample][protein]:
            for row in mid[sample][protein][label][1:]:
                totals.append(row[1])
    else:
        for p in mid[sample]:
            for label in mid[sample][p]:
                for row in mid[sample][p][label][1:]:
                    totals.append(row[1])
    # Get relative size
    ti = float(ti)
    o_max = np.amax(totals)
    o_min = np.amin(totals)
    n_max = 1000
    n_min = n_max/10
    o_range = (o_max - o_min)
    n_range = (n_max - n_min)
    new_size = (((ti - o_min)*n_range)/o_range)+n_min
    if np.isnan(new_size): new_size = n_min
    return new_size


def site_spef(mid, show = False, to_print = True):
    """ Plots site-specific deamidation plot
    """
    # Get data and sort by sample
    data = np.array(sorted(ss_wrangle(mid), key=lambda row:row[0]))

    # None for now
    prev_sample, prev_protein = data[0][0:2]

    # Data to print later
    data_to_print = []

    single_sample = False
    if len(mid.keys()) == 1:
        single_sample = True

    for line in data:

        sample, protein, label, hf, rmi, ti = line

        # Get a relative size for the points
        size = get_relative_size(ti, mid, sample, protein, single_sample)

        # We want intact intensity, not deamidated intensity
        rmi = 1 - float(rmi)

        if to_print:
            data_to_print.append([hf, rmi, size, sample, protein])

        prev_sample = sample
        prev_protein = protein

    # Send to printing?
    if to_print:
        save_csv_results(data_to_print, "Site-Specific")


def save_fine_bulk(data):
    """ Saves bulk deamidation data pre-averaging
    """
    results_dir = os.path.join(data_folder, "Results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    title = "Bulk_fine_grain_results.csv"
    path = os.path.join(results_dir, title)
    with open(path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sample", "Protein", "AA", "RelNonDeam"])
        writer.writerows(data)
    csvfile.close()
    print "%s saved in %s" % (title, results_dir)


def save_csv_results(data, method):
    """Saves the results needed to recreate the plots
    """
    results_dir = os.path.join(data_folder, "Results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    title = "%s_results.csv" % method
    path = os.path.join(results_dir, title)
    with open(path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        if method == "Site-Specific":
            writer.writerow(["Half-time", "RelNonDeam", "Size", "Sample", "Protein"])
        else:
            writer.writerow(["Sample", "Protein", "NNonDeam", "QNonDeam"])
        writer.writerows(data)
    csvfile.close()
    print "%s saved in %s" % (title, results_dir)


def read_protein_list(protein_list_file):
    """ Reads a list of relevant proteins in order to filter data later
    """
    f = open(protein_list_file, "r")
    protein_list = f.read()
    return protein_list


data_folder = ""
def main():
    global data_folder
    try:
        data_folder = sys.argv[1]
    except IndexError as e:
        print "Specify path to data"
    protein_list = []
    if len(sys.argv) > 2:
        protein_list_file = sys.argv[2]
        protein_list = read_protein_list(protein_list_file)
    total_data = mq(data_folder, protein_list, filter_con = True)
    mid_classic, mid_ss = get_mid(total_data)
    bulk_deam(mid_classic, show = False, to_print = True)
    site_spef(mid_ss, show = False, to_print = True)


main()
