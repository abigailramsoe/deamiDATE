#!/usr/bin/python
import csv


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
        e_exp = e_headers.index("Experiment")

        for row in e_reader:
            # Experiment, seq, mods, intensity, ID, protein
            e_data.append([row[e_exp], row[e_seq], row[e_mod_seq], row[e_i],
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
                                        mod_label = "%s%i-%s" % (mods[i-1],
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

mq("/home/abby/Documents/1. DeamiDATE/Test files")
