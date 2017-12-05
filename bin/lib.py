import sys
import json
import numpy as np
from random import randint
from random import uniform

class SampleStore:
    def __init__(self, debug=False):
        self.debug = debug
        self.samples = {}
        self.controls = set()
        self.cases = set()
        self.patterns = {}

    def add_sample(self, s):
        if not self.isSampleObject(s):
            print("Invalid class type of s!", sys.stderr())
            exit()
        self.samples[s.id] = s
        if s.is_control():
            self.controls.add(s.id)
        else:
            self.cases.add(s.id)

        if not s.pattern in self.patterns:
            self.patterns[s.pattern] = 0
        self.patterns[s.pattern] += 1

    def get_controls(self):
        print('get the controls')

    def get_cases(self):
        print("get the cases")

       # print(self.patterns)



    def print_summary(self):
        print("---------------- Sample store summary ----------------")
        print("Total number of samples: %s" % len(self.samples.keys()))
        print("Case count: %s" % len(self.cases))
        print("Control count: %s" % len(self.controls))
        print("Pattern summary:")
        for p, v in sorted(self.patterns.items()):
            print("%s: %s" % (p, v))


    def isSampleObject(self, object):
        """Return true if the object is a class.

        Class objects provide these attributes:
            __doc__         documentation string
            __module__      name of module in which this class was defined"""
        return isinstance(object, Sample)


class Sample:
    def __init__(self, id, phenotype, pattern):
        self.id = id
        self.phenotype = phenotype
        self.pattern = pattern
        self.genotype = {}

    def print_sample(self):
        print("ID: %s Phenotype: %s Pattern: %s Genotype: %s" % (self.id, self.phenotype, self.pattern, self.genotype))

    def is_control(self):
        if self.phenotype == "control":
            return True
        return False

    def is_case(self):
        if self.is_control():
            return False
        return True

    def add_genotype(self, genotype):
        for v in genotype:
            gt = genotype[v].split("/")
            self.genotype[v] = gt



class Manifest:
    def __init__(self, file, debug=False):
        self.file = file
        self.debug = debug
        self.manifest = self.read_manifest(file)

        self.variants = {}
        self.regions = {}
        self.patterns = {}

        self.get_variants()
        self.get_regions()

        # Need to be able to return regions
        # Need to be able to return variants
        # Need to be able to return patterns


    def read_manifest(self, manifest_file):
        if self.debug:
            print("Reading in manifest", file=sys.stderr)
        manifest = json.load(open(manifest_file))
        return manifest

    def get_variants(self):
        if self.debug:
            print("Parsing manifest variants", file=sys.stderr)
        variants = set()
        phenotype_patterns = self.manifest['phenotype_patterns']
        for p in phenotype_patterns:
            new_pattern = Pattern(p, debug=self.debug)
            new_pattern.case_proportion = self.manifest['phenotype_patterns'][p]['case_proportion']
            new_pattern.control_proportion = self.manifest['phenotype_patterns'][p]['control_proportion']
            new_pattern.protective_variants = self.manifest['phenotype_patterns'][p]['protective_variants']
            new_pattern.type = self.manifest['phenotype_patterns'][p]['type']

            for v in self.manifest['phenotype_patterns'][p]['variants']:
                chrom = self.manifest['phenotype_patterns'][p]['variants'][v]['chromosome']
                pos = self.manifest['phenotype_patterns'][p]['variants'][v]['position']
                ref = self.manifest['phenotype_patterns'][p]['variants'][v]['ref']
                alt = self.manifest['phenotype_patterns'][p]['variants'][v]['alt']
                hom_ref_freq = self.manifest['phenotype_patterns'][p]['variants'][v]['hom_ref_freq']
                het_freq = self.manifest['phenotype_patterns'][p]['variants'][v]['het_freq']
                hom_alt_freq = self.manifest['phenotype_patterns'][p]['variants'][v]['hom_alt_freq']
                new_variant = Variant(chrom, pos, ref, alt, hom_ref_freq, het_freq, hom_alt_freq, debug=self.debug)
                if not new_variant.id in self.variants:
                    self.variants[new_variant.id] = set()
                self.variants[new_variant.id].add(p)
                new_pattern.add_variant(new_variant)

            self.patterns[p] = new_pattern
        return variants

    def get_regions(self):
        if self.debug:
            print("Parsing manifest regions", file=sys.stderr)
        manifest_regions = self.manifest['regions']
        for r in manifest_regions:
            chrom = self.manifest['regions'][r]['chromosome']
            start = self.manifest['regions'][r]['start']
            end = self.manifest['regions'][r]['end']
            new_region = Region(chrom, start, end)
            self.regions[new_region.id] = new_region


    def pattern_proportions(self):
        case_total = 0
        control_total = 0
        case_patterns = {}
        control_patterns = {}

        for p in self.manifest.patterns.keys():  # this will break.   Should be the manifests
            case = self.patterns[p].case_proportion
            control = self.patterns[p].control_proportion

            case_total += case
            control_total += control

            if case > 0:
                case_patterns[p] = case
            if control > 0:
                control_patterns[p] = control

        if case_total > 1:
            print("Case proportion values invalid", file=sys.stderr)
        if control_total > 1:
            print("Control proportion values invalid", file=sys.stderr)

        if case_total < 1:
            case_patterns['random'] = 1 - case_total

        if control_total < 1:
            control_patterns['random'] = 1 - control_total

class Variant:
    def __init__(self, chrom, pos, ref, alt, hom_ref_freq, het_freq, hom_alt_freq, debug=False):
        self.debug = False
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.pattern_frequency = None
        self.case_frequency = None
        self.type = None
        self.id = "%s-%s-%s-%s" % (chrom, pos, ref, alt)
        self.hom_ref_freq = hom_ref_freq
        self.het_freq = het_freq
        self.hom_alt_freq = hom_alt_freq


class Region:
    def __init__(self, chrom, start, end, debug=False):
        self.debug = False
        self.chrom = chrom
        self.start = start
        self.end = end
        self.id = "%s-%s-%s" % (chrom, start, end)

class Pattern:
    def __init__(self, id, case_proportion=None, control_protportion=None, protective_variants=None, debug=False):
        self.id = id
        self.debug = debug
        self.case_proportion = case_proportion
        self.control_proportion = control_protportion
        self.protective_variants = protective_variants
        self.type = None # Could add another object for these.
        self.variants = {}

    def add_variant(self, v):
        self.variants[v.id] = v

    def print_summary(self):
        print("PATTERN %s - CASE %s - CONTROL %s -  TYPE %s" % (self.id, self.case_proportion, self.control_proportion, self.type))
        #print(self.protective_variants)
        #print(self.variants)


    # Returns a dictionary of variants and their genotype.
    def generate_variants(self, phenotype="case"):
        # Generate genotype based on the stored pattern.  Should be random each time.
        # Check the stored variants.  Create a list of the variants and a list of the case frequencies
        if self.type == "phenotype_linked_variant":
            return self.phenotype_linked_variant()
        elif self.type == "phenotype_variant_with_protective_SNV":
            return self.phenotype_variant_with_protective_SNV(phenotype)
        elif self.type == "synergistic_phenotype_linked":
            return self.synergistic_phenotype_linked(phenotype)

    def phenotype_linked_variant(self):
        #print("phenotype_linked_variant")
        # Create a list of the variants and a list of the case frequencies - make a function
        variants, freqs = self.variant_pattern_freq_lists()
        # sample a single variant from this list
        variant = np.random.choice(variants, 1, p=freqs).tolist()[0]
        # Check het-hom proportions
        genotype = self.generate_variant_genotype(variant)
        result = {variant: genotype}
        return result

    def generate_variant_genotype(self, variant):
        #print("Generating genotype")
        genotypes = ["0/0", "0/1", "1/1"]
        freqs = [self.variants[variant].hom_ref_freq, self.variants[variant].het_freq, self.variants[variant].hom_alt_freq]
        if None in freqs:
            freqs[0] = 0
            freqs[1] = uniform(0.0, 1.0)
            freqs[2] = 1 - freqs[1]
        genotype = np.random.choice(genotypes, 1, p=freqs).tolist()[0]
        return genotype


    def phenotype_variant_with_protective_SNV(self, phenotype):
        #print("phenotype_variant_with_protective_SNV")
        ## If a control genotype is being generated, make sure the protective variant is set
        # If a case genotype is being generated, no protective variant, but make sure a phenotype variant is set
        # Identify all non-protective variants
        variants, freqs = self.variant_pattern_freq_lists()
        phenotype_variant = np.random.choice(variants, 1, p=freqs).tolist()[0]
        phenotype_genotype = self.generate_variant_genotype(phenotype_variant)
        if phenotype == "control":
            protective_genotype = self.generate_variant_genotype(self.protective_variants[0]) # assuming only one protective variant
        else:
            protective_genotype = "0/0"

        result = {
            phenotype_variant: phenotype_genotype,
            self.protective_variants[0]: protective_genotype
        }
        return result

    def synergistic_phenotype_linked(self, phenotype):
        #print('synergistic_phenotype_linked')
        # in this pattern all variants in the pattern must be present for the phenotype to occur
        # maybe randomly select one variant to be missing for controls then randomly set genotypes for the others
        # for cases just randomly set the genotypes
        # Basically just protective module flipped.   I should rename "protective" to "linked"
        variants, freqs = self.variant_pattern_freq_lists()
        phenotype_variant = np.random.choice(variants, 1, p=freqs).tolist()[0]
        phenotype_genotype = self.generate_variant_genotype(phenotype_variant)
        if phenotype == "case":
            linked_genotype = self.generate_variant_genotype(self.protective_variants[0]) # assuming only one protective variant
        else:
            linked_genotype = "0/0"

        result = {
            phenotype_variant: phenotype_genotype,
            self.protective_variants[0]: linked_genotype
        }
        return result



    def variant_pattern_freq_lists(self, randomize_null=True):
        variants = []
        freqs = []
        for v in self.variants:
            if self.protective_variants is not None and v in self.protective_variants:
                continue
            variants.append(self.variants[v].id)
            freqs.append(self.variants[v].pattern_frequency)


        if randomize_null:
            # The assumption is that either all variants have a specified frequency of none of them do
            if None in freqs:
                new_freqs = []
                sum = 0
                for v in variants:
                    r = randint(1, 100)
                    new_freqs.append(r)
                    sum += r
                for i in range(len(new_freqs)):
                    new_value = new_freqs[i] / sum
                    freqs[i] = new_value
        return variants, freqs







