#!/usr/bin/python
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import sys
import os

def mq(folder):
    """ Reads from evidence and peptide MQ results files and returns a
    dictionaryof the form:
    total_data[sample][protein] = [intensity, pre, start, seq,
                                   end, post, sample, mods]
    """
    print "Reading file from MQ"

    # Evidence and peptide folders
    evi_file = folder + "/evidence.txt"
    pep_file = folder +"/peptides.txt"

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
                    possible_mods = ["de", "hy", "ox"] #TODO: Add more mods
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


def peaks(folder):
    #TODO
    place = True


def gpm(folder):
    #TODO
    place = True


def mascot(mascot):
    #TODO
    place = True


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
                            if 0 < aa_i < len(seq)-1: # As long as the AA has two neighbours
                                x = seq[aa_i+1] # Look right
                                y = seq[aa_i-1] # Look left
                                combo = "%s %s %s" % (y, aa, x) # SS label
                                label = "%s %s" % (cur_pos, aa) # Classic label

                                # Get half time from Pandas dataframe
                                if (y in rob.index) and (x in rob.columns):
                                    hf = rob.ix[y, x]
                                total_i += intensity
                                # If this position was modified, its mod label would be this
                                q_mod = "%s-%i-de" % (aa, cur_pos)
                                # If that label is in the mod list for this pos
                                if q_mod in mods:
                                    # Increase mod intensity
                                    mod_i += intensity

                                # Add to SS MID
                                if combo in mid_ss[sample][protein]:
                                    mid_ss[sample][protein][combo].append([mod_i, total_i])
                                else: mid_ss[sample][protein][combo] = [hf, [mod_i, total_i]]

                                # Add to classic MID
                                if label in mid_classic[sample][protein]:
                                    old_hf, [old_mod_i, old_total_i] = mid_classic[sample][protein][label]
                                    # Double check we're not trying to update the half time
                                    try:
                                        assert old_hf == hf
                                    except AssertionError:
                                        print "Half times not equal"
                                    mod_i += old_mod_i
                                    total_i += old_total_i
                                mid_classic[sample][protein][label] = [hf,[mod_i, total_i]]
                        # On to the next AA
                        cur_pos += 1

    return mid_classic, mid_ss


def calc_deam(mid):
    """ Calculates bulk deamidation per sample, per protein
    Returns a list of [sample, protein, rel_asn, rel_gln]
    Where rel_asn and _gln are mod_intensity/total_intensity
    """
    relative = []
    for sample in mid:
        for protein in mid[sample]:
            asn_m, asn_t, gln_m, gln_t = 0, 0, 0, 0
            rel_asn, rel_gln = -1, -1
            for label, val in mid[sample][protein].items():
                hl = val[0]
                mod, total = val[1]
                cur_pos, aa = label.split(" ")
                if "N" == aa: # Asn
                    asn_m += mod
                    asn_t += total
                if "Q" == aa: # Gln
                    gln_m += mod
                    gln_t += total
            # Calculate relative amounts
            if asn_t > 0:
                rel_asn = 1 - (asn_m/asn_t)
            else: print "No Asn %s %s" % (sample, protein)
            if gln_t > 0:
                rel_gln = 1 - (gln_m/gln_t)
            else: print "No Gln %s %s" % (sample, protein)
            # Each sample relative amounts
            relative.append([sample, protein, rel_asn, rel_gln])
    print "Bulk deamidation calculated"
    return relative


def bulk_deam(mid, show = False, debug = False):
    """ Plots bulk deamidation
    """
    relative = np.array(calc_deam(mid))
    index = np.arange(len(relative[:,0]))
    relative = np.array(sorted(relative, key=lambda row:row[0])) # Sort by sample

    width = .35
    fig, ax = plt.subplots()

    # Asn and Gln bars
    asn_bars = ax.bar(index, relative[:,2], width, color='blue')
    gln_bars = ax.bar(index + width, relative[:,3], width, color='red')

    # Make the bar names protein and sample
    names = ["%s %s" % (x, y) for x, y in zip(relative[:,0], relative[:,1])]
    # Lables
    ax.set_ylabel('% (Asn, Gln)')
    ax.set_xticks(index + width / 2)
    ax.set_xticklabels(names)
    fig.autofmt_xdate()

    # Legend
    ax.legend((asn_bars[0], gln_bars[0]), ('% Asn', '% Gln'), loc="lower right")

    # Limits
    xmin, xmax, ymin, ymax = plt.axis()
    plt.axis((xmin, xmax, ymin, 1.0))

    # Title
    plot_title = "Deamidation"
    plt.title(plot_title)
    save_plots("Bulk")
    if debug: print relative
    if show: plt.show()


def save_plots(method):
    """ Save plots to a Results dir, which is created or located inside the
    data directory given as args
    """
    results_dir = "%s/Results" % data_folder
    if not os.path.exists(results_dir):
		os.makedirs(results_dir)
    title = "%s_plot.png" % method
    path = "%s/%s" % (results_dir, title)
    plt.savefig(path)
    print "%s saved in %s" % (title, results_dir)



data_folder = ""
def main():
    global data_folder
    IMPLEMENTED_SOFTWARE = {"MQ": mq, "PEAKS": peaks,
                        "GPM": gpm, "MASCOT": mascot}
    try:
        software = sys.argv[1].upper()
        data_folder = sys.argv[2]
    except IndexError as e:
        print "Specify software type, followed by path to data"

    total_data = IMPLEMENTED_SOFTWARE[software](data_folder)
    mid_classic, mid_ss = get_mid(total_data)
    bulk_deam(mid_classic, show = False, debug = False)


main()
