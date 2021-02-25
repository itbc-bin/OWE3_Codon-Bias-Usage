import regex as re


def extract_data(file):
    """Reads Fasta file, saves the header and sequence
    for every gene in a list

    return: genes - 2d list
    """
    genes = []
    seq = ""
    with open(file) as data:
        for line in data:
            if re.search("^>", line):
                if seq != "":
                    genes.append([header, seq])
                    seq = ""
                header = line.replace("\n", "")
            else:
                seq += line.replace("\n", "")
        genes.append([header, seq])

    return genes


def define_groups():
    """Divides the genes into two groups, one for surface proteins
    and the other for internal proteins

    return: surface_gene - list
    return: internal_gene - list
    """


def codons():
    """Makes a dictionary which can be used for codon analysis

    return: codon_dict - dict
    """
    codon_dict = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                  'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                  'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
                  'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
                  'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
                  'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
                  'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
                  'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
                  'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
                  'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
                  'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
                  'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
                  'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
                  'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
                  'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
                  'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
                  }

    return codon_dict


def calc_bias():
    """

    return:
    """


def calc_usage(seq, codons):
    """Calculates What codons are used in the sequence and what

    return: Values
    """

    for header, sequence in seq:
        for codon in range(0, len(sequence), 3):
            print(sequence[codon:codon + 3])


    return ""

def make_graph():
    """

    return:
    """


if __name__ == '__main__':
    Aconitate_genes = "sequenties.txt"
    genes_four_organisms = extract_data(Aconitate_genes)
    codons = codons()
    usage_value = calc_usage(genes_four_organisms, codons)
