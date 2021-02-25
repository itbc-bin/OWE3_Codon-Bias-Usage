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


def calc_bias():
    """

    return:
    """


def calc_usage(seq):
    """Calculates What codons are used in the sequence and what

    return: Values
    """

    codon_dict = codons()
    for header, sequence in seq:
        header = header[:header.find(" ")]
        header = header.replace(">", "")
        header = header.replace(":", "_")
        codon_code_count, amino_acid_count = count()
        for position in range(0, len(sequence), 3):
            codon = sequence[position:position + 3]
            codon_code_count[codon] += 1
            amino_acid_count[codon_dict[codon]] += 1

        globals()[header] = amino_acid_count, codon_code_count
    hsa_48_percent = []
    rcn_112186046_percent = []
    aasc_A4S02_13295_percent = []
    acij_JS278_01766_percent = []
    for codon in codon_dict:
        percent = round(hsa_48[1][codon] / hsa_48[0][codon_dict[codon]] * 100, 2)
        hsa_48_percent += [[codon, percent]]
        percent = round(rcn_112186046[1][codon] / rcn_112186046[0][codon_dict[codon]] * 100, 2)
        rcn_112186046_percent += [[codon, percent]]
        percent = round(aasc_A4S02_13295[1][codon] / aasc_A4S02_13295[0][codon_dict[codon]] * 100, 2)
        aasc_A4S02_13295_percent += [[codon, percent]]
        percent = round(acij_JS278_01766[1][codon] / acij_JS278_01766[0][codon_dict[codon]] * 100, 2)
        acij_JS278_01766_percent += [[codon, percent]]
    print(hsa_48_percent)
    print(rcn_112186046_percent)
    print(aasc_A4S02_13295_percent)
    print(acij_JS278_01766_percent)

    return hsa_48_percent, rcn_112186046_percent, \
           aasc_A4S02_13295_percent, acij_JS278_01766_percent


def make_graph():
    """

    return:
    """


if __name__ == '__main__':
    Aconitate_genes = "sequenties.txt"
    genes_four_organisms = extract_data(Aconitate_genes)
    usage_value = calc_usage(genes_four_organisms)
