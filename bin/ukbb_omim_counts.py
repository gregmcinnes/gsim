"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from __future__ import print_function
import argparse
import sys

class UKBBOmimCounts(object):
    def __init__(self, disease_ontology, ukbb_counts, debug):
        self.debug = debug
        self.run(disease_ontology, ukbb_counts)

    def run(self, do_file, ukbb_file):
        # load disease ontology into dictionary
        terms = {}
        do_icds = {}
        with open(do_file) as f:
            for line in f:
                term = Term()
                term.parse_line(line)
                for i in term.clean_icds():
                    if not i in do_icds:
                        do_icds[i] = set()
                    do_icds[i].add(term.id)
                terms[term.id] = term
        f.close()
        # For each line in ukbb counts try to map to dictionary and if matches, add count
        with open(ukbb_file) as f:
            for line in f:
                icd, count = line.rstrip().split(",")
                if icd in do_icds:
                    for t in do_icds[icd]:
                        terms[t].ukbb_count += int(count)
                        terms[t].print_summary()

class Term:
    def __init__(self):
        self.id = None
        self.omim = set()
        self.mesh = set()
        self.icd10 = set()
        self.name = None
        self.ukbb_count = 0

    def print_summary(self):
        #print("%s\t%s\t%s\t%s\t%s" % (self.id, self.name, ",".join(self.icd10), ",".join(self.omim), ",".join(self.mesh)))
        print("%s\t%s\t%s\t%s" % (self.name, ",".join(self.icd10), ",".join(self.omim), self.ukbb_count))

    def unique_omim(self):
        if len(self.omim) == 1:
            return True
        return False

    def has_omim(self):
        if len(self.omim) > 0:
            return True
        return False

    def has_icd(self):
        if len(self.icd10) > 0:
            return True
        return False

    def parse_line(self, line):

        fields = line.rstrip().split("\t")
        self.id = fields[0]
        self.name = fields[1]
        #self.mesh = fields[4].split(",")
        self.omim = fields[3].split(",")
        self.icd10 = fields[2].split(',')

    def clean_icds(self):
        icds = set()
        for i in self.icd10:
            icds.add(i.replace(".",""))
        return icds
"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("-d", "--disease_ontology", help="Input")
    parser.add_argument("-u", "--ukbb_counts", help="Input")
    parser.add_argument("--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    UKBBOmimCounts(options.disease_ontology, options.ukbb_counts, options.debug)

