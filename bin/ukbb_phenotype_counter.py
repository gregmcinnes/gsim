"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from __future__ import print_function
import argparse

import sys

class Counter(object):
    def __init__(self, file, debug):
        self.debug = debug
        self.run(file)

    def run(self, file):
        code_counts = {}
        with open(file) as f:
            f.readline()
            for line in f:
                codes = line.rstrip().replace('"', '').split(",")
                codes = list(filter(None, codes))
                for c in codes[1:]:
                    if not c in code_counts:
                        code_counts[c] = 0
                    code_counts[c] += 1
        for c in code_counts:
            print("%s,%s" % (c, code_counts[c]))

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("-f", "--file", help="Input")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    Counter(options.file, options.debug)

