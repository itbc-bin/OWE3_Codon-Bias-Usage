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
                    seq = seq.lower()
                    genes.append([header, seq])
                    seq = ""
                header = line.replace("\n", "")
            else:
                seq += line.replace("\n", "")
        seq = seq.lower()
        genes.append([header, seq])

    return genes


def process(virus):
    """

    :return:
    """
    virus_data = extract_data(virus)
    virus_groups = define_groups(virus_data)
    virus_bias = calc_bias(virus_groups)

    return virus_bias

def define_groups(virus):
    """Divides the genes into two groups, one for surface proteins
    and the other for internal proteins

    return: surface_gene - list
    return: internal_gene - list
    """
    group_internal = []
    group_envelop = []
    for gene in virus:
        for line in gene:
            if re.search(">", line):
                header = line.split(" ")
                for attribute in header:
                    attribute = attribute.lower()
                    if "protein=" in attribute \
                            and "env" in attribute:
                        group_envelop.append(gene)
                    elif "protein=" in attribute \
                            and "env" not in attribute:
                        group_internal.append(gene)

    return group_envelop, group_internal


def codons():
    """Makes a dictionary which can be used for codon analysis

    return: codon_dict - dict
    """
    codon_dict = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                  'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                  'tta': 'L', 'tca': 'S', 'taa': 'X', 'tga': 'X',
                  'ttg': 'L', 'tcg': 'S', 'tag': 'X', 'tgg': 'W',
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
                  'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G',
                  }
    return codon_dict


def count():
    """

    :return:
    """
    codon_code_count = dict(ttt=0, tct=0, tat=0, tgt=0, ttc=0, tcc=0,
                            tac=0, tgc=0, tta=0, tca=0, taa=0, tga=0,
                            ttg=0, tcg=0, tag=0, tgg=0, ctt=0, cct=0,
                            cat=0, cgt=0, ctc=0, ccc=0, cac=0, cgc=0,
                            cta=0, cca=0, caa=0, cga=0, ctg=0, ccg=0,
                            cag=0, cgg=0, att=0, act=0, aat=0, agt=0,
                            atc=0, acc=0, aac=0, ata=0, agc=0, aca=0,
                            aaa=0, aga=0, atg=0, acg=0, aag=0, agg=0,
                            gtt=0, gct=0, gat=0, ggt=0, gtc=0, gcc=0,
                            gac=0, ggc=0, gta=0, gca=0, gaa=0, gga=0,
                            gtg=0, gcg=0, gag=0, ggg=0,
                            )

    amino_acid_count = dict(F=0, L=0, I=0, M=0, V=0, S=0, P=0, T=0, A=0,
                            Y=0, X=0, H=0, Q=0, N=0, K=0, D=0, E=0, C=0,
                            W=0, R=0, G=0,
                            )

    return codon_code_count, amino_acid_count


def calc_bias(groups):
    """

    return:
    """
    codon_bias_percent = []
    codon_dict = codons()
    all_percent = []
    codon_code_count, amino_acid_count = count()
    for group in groups:
        for header, sequence in group:
            for position in range(0, len(sequence), 3):
                codon = sequence[position:position + 3]
                codon_code_count[codon] += 1
                amino_acid_count[codon_dict[codon]] += 1
    for codon in codon_dict:
        try:
            percent = round(codon_code_count[codon] /
                            amino_acid_count[codon_dict[codon]]
                            * 100, 2)
        except ZeroDivisionError:
            percent = "amino acid " + codon_dict[codon] +\
                    " is not used in sequence"
        all_percent += [codon, percent]
    codon_bias_percent += all_percent

    return codon_bias_percent

def calc_usage(seq):
    """Calculates What codons are used in the sequence and what

    return: Values
    """
    codon_usage_percent = []
    codon_dict = codons()
    codon_percentages = []
    for header, sequence in seq:
        header = header[:header.find(" " or "_")]
        header = header.replace(">", "")
        header = header.replace(":", "_")
        codon_code_count, amino_acid_count = count()
        for position in range(0, len(sequence), 3):
            codon = sequence[position:position + 3]
            codon_code_count[codon] += 1
            amino_acid_count[codon_dict[codon]] += 1
        globals()[header] = []
        for codon in codon_dict:
            try:
                percent = round(codon_code_count[codon] /
                                amino_acid_count[codon_dict[codon]]
                                * 100, 2)
            except ZeroDivisionError:
                percent = "amino acid", codon_dict[codon], \
                          "is not used in sequence"
            codon_percentages += [[codon, percent]]
        codon_usage_percent += header, globals()[header]

    return codon_usage_percent


def make_graph():
    """

    return:
    """


if __name__ == '__main__':
    Aconitate_genes = "sequenties.txt"
    hiv1 = "nucleo_hiv_1.txt"
    hiv2 = "nucleo_hiv_2.txt"
    siv1 = "nucleo_siv.txt"
    siv2 = "nucleo_sivmnd2.txt"

    hiv1_bias = process(hiv1)
    hiv2_bias = process(hiv2)
    siv1_bias = process(siv1)
    siv2_bias = process(siv2)

    genes_four_organisms = extract_data(Aconitate_genes)
    usage_value = calc_usage(genes_four_organisms)
