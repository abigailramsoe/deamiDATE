#!/usr/bin/python
import csv
import numpy as np

def mq(folder):

    print "Reading file from MQ"
    evi_file = folder + "/evidence.txt"
    pep_file = folder +"/peptides.txt"
    first = True
    e_data = []
    p_data = {}

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
            e_data.append([row[e_exp], row[e_seq], row[e_mod_seq], row[e_i], row[e_pep_id], row[e_protein]])

    with open(pep_file) as csvfile2:
        p_reader = csv.reader(csvfile2, delimiter='\t')
        p_headers = next(p_reader)
        p_pre = p_headers.index("Amino acid before")
        p_post  = p_headers.index("Amino acid after")
        p_start = p_headers.index("Start position")
        p_end = p_headers.index("End position")
        p_ID = p_headers.index("id")

        for row in p_reader:
            p_data[row[p_ID]] = [row[p_pre], row[p_post], row[p_start], row[p_end]]

    c = 1

    total_data = {}


    for row in e_data:
        sample, seq, mods, intensity, ID, protein = row
        print mods
        if sample not in total_data:
            total_data[sample] = {}
        if intensity != "" and "REV" not in protein:
            if ID in p_data:
                pre, post, start, end = p_data[ID]
                if pre != "" and start != "":
                    mods_list = []
                    possible_mods = ["de", "hy", "ox"]
                    for m in possible_mods:
                        if m in mods:
                            actual_pos = 1
                            for i in xrange(1, len(mods)):
                                if mods[i] == "(":
                                    if mods[i+1] == m[0] and mods[i+2] == m[1]:
                                        mod_label = "%s%i#0.98402" % (mods[i-1], (int(start)+actual_pos-2))
                                        mods_list.append(mod_label)
                                if mods[i].isupper(): actual_pos += 1


                    info = [float(intensity), pre, int(start),seq,int(end),post, sample, mods_list]


                    #protein = xp #PARCHMENT
                    if protein not in total_data[sample]:
                        total_data[sample][protein] = [info]
                    else:
                        total_data[sample][protein].append(info)



        #print p_data[row[-1]], c
        c += 1
        #print ""

                    #print "reader2", [float(intensity), pre, int(start),seq,int(end),post,mods_list],sample,e_protein

    return total_data

mq("/home/abby/Documents/1. DeamiDATE/Test files")
