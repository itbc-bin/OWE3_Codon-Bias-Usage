import regex as re
import matplotlib.pyplot as plt


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
    virus_env_usage = calc_usage(virus_groups[0])
    virus_internal_usage = calc_usage(virus_groups[1])
    virus_usage = [virus_env_usage, virus_internal_usage]
    #print(virus_internal_usage)

    return virus_bias, virus_usage


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
                      " is not used in this sequence"
        all_percent += [codon, percent]
    codon_bias_percent += all_percent

    return codon_bias_percent


def calc_usage(seq_group):
    """Calculates What codons are used in the sequence and what

    return: Values
    """
    codon_usage_percent = []
    codon_dict = codons()

    for header, sequence in seq_group:
        codon_percentages = []
        header = header[:header.find(" " or "_")]
        header = header.replace(">", "")
        header = header.replace(":", "_")
        header = header[:15]
        codon_code_count, amino_acid_count = count()
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
                percent = "amino acid", codon_dict[codon], \
                          "is not used in sequence"
            codon_percentages += [[codon, percent]]
        codon_usage_percent += header, codon_percentages
    #print(codon_usage_percent)
    return codon_usage_percent


def lijsten_organismen(usage_value):
    #print(usage_value)
    hsa = usage_value[1]
    #print(hsa)
    codon = []
    hsa_percentage = []
    rcn_percentage = []
    aasc_percentage = []
    acij_percentage = []
    for line in hsa:
        hsa_percentage.append(line[1])
        codon.append(line[0])
    print(hsa_percentage)
    rcn = usage_value[3]
    for line in rcn:
        rcn_percentage.append(line[1])

    aasc = usage_value[5]
    for line in aasc:
        aasc_percentage.append(line[1])

    acij = usage_value[7]
    for line in acij:
        acij_percentage.append(line[1])

    return codon, hsa_percentage, rcn_percentage, \
           aasc_percentage, acij_percentage

def lijsten_hiv1(hiv1_usage, hiv1_bias):
    hiv1_usage_1 = hiv1_usage[0]
    hiv1_usage_1 = hiv1_usage_1[1]

    hiv1_usage_2 = hiv1_usage[1]
    hiv1_usage_2 = hiv1_usage_2[1]

    hiv1_usage_1_percentage = []
    hiv1_usage_1_codon = []
    hiv1_usage_2_percentage = []
    hiv1_usage_2_codon = []

    for line in hiv1_usage_1:
        hiv1_usage_1_percentage.append(line[1])
        hiv1_usage_1_codon.append(line[0])

    for line in hiv1_usage_2:
        hiv1_usage_2_percentage.append(line[1])
        hiv1_usage_2_codon.append(line[0])

    return hiv1_usage_1_codon, hiv1_usage_1_percentage, hiv1_usage_2_codon, hiv1_usage_2_percentage


def make_graph_hsa(codon, hsa_percentage, codon_dict):
    """

    return:
    """
    plt.figure(num=1, figsize=(15, 10))
    plt.bar(codon, hsa_percentage, color="darkgreen")
    plt.ylabel("percentage(%)")
    plt.xlabel("codons")
    plt.xticks(rotation=90, fontsize=12)
    plt.title("hsa")

    #labels = list(codon_dict.values())
    #handles = [plt.Rectangle((0, 0), 1, 1)]
    #plt.legend(handles, labels)
    plt.show()


def make_graph_rcn(codon, rcn_percentage,codon_dict):
    """

    return:
    """
    plt.figure(num=1, figsize=(15, 10))
    plt.bar(codon, rcn_percentage, color="maroon")
    plt.ylabel("percentage(%)")
    plt.xlabel("codons")
    plt.xticks(rotation=90, fontsize=12)
    plt.title("rcn")

    labels = list(codon_dict.values())
    handles = [plt.Rectangle((0, 0), 1, 1)]
    plt.legend(handles, labels)
    plt.show()


def make_graph_aasc(codon, aasc_percentage, codon_dict):
    """

    return:
    """
    plt.figure(num=1, figsize=(15, 10))
    plt.bar(codon, aasc_percentage, color="midnightblue")
    plt.ylabel("percentage(%)")
    plt.xlabel("codons")
    plt.xticks(rotation=90, fontsize=12)
    plt.title("aasc")

    labels = list(codon_dict.values())
    handles = [plt.Rectangle((0, 0), 1, 1)]
    plt.legend(handles, labels)
    plt.show()


def make_graph_acij(codon, acij_percentage,codon_dict):
    """

    return:
    """
    plt.figure(num=1, figsize=(15, 10))
    plt.bar(codon, acij_percentage, color="darkslategrey")
    plt.ylabel("percentage(%)")
    plt.xlabel("codons")
    plt.xticks(rotation=90, fontsize=12)
    plt.title("acij")

    labels = list(codon_dict.values())
    handles = [plt.Rectangle((0, 0), 1, 1, color="darkslategrey")]
    plt.legend(handles, labels)
    plt.show()

def make_graph_hiv1u1(hiv1_usage_1_codon, hiv1_usage_1_percentage, codon_dict):
    """

    return:
    """
    plt.figure(num=1, figsize=(15, 10))
    plt.bar(hiv1_usage_1_codon, hiv1_usage_1_percentage, color = "pink")
    plt.ylabel("percentage(%)")
    plt.xlabel("codons")
    plt.xticks(rotation = 90, fontsize = 12)
    plt.title("hiv1 usage 1")

    labels = list(codon_dict.values())
    handles = [plt.Rectangle((0, 0), 1, 1, color = "darkslategrey")]
    plt.legend(handles, labels)
    plt.show()

if __name__ == '__main__':
    Aconitate_genes = "sequenties.txt"
    hiv1 = "nucleo_hiv_1.txt"
    hiv2 = "nucleo_hiv_2.txt"
    siv1 = "nucleo_siv.txt"
    siv2 = "nucleo_sivmnd2.txt"
    codon_dict = codons()

    hiv1_bias, hiv1_usage = process(hiv1)
    hiv2_bias, hiv2_usage = process(hiv2)
    siv1_bias, siv1_usage = process(siv1)
    siv2_bias, siv2_usage = process(siv2)

    genes_four_organisms = extract_data(Aconitate_genes)
    usage_value = calc_usage(genes_four_organisms)
    print(usage_value)
    codon, hsa_percentage, rcn_percentage, \
           aasc_percentage, acij_percentage = lijsten_organismen(usage_value)

    hiv1_usage_1_codon, hiv1_usage_1_percentage, hiv1_usage_2_codon, hiv1_usage_2_percentage = lijsten_hiv1(hiv1_usage,
                                                                                                            hiv1_bias)

    make_graph_hsa(codon, hsa_percentage, codon_dict)
    make_graph_rcn(codon, rcn_percentage, codon_dict)
    make_graph_aasc(codon, aasc_percentage, codon_dict)
    make_graph_acij(codon, acij_percentage, codon_dict)
    make_graph_hiv1u1(hiv1_usage_1_codon, hiv1_usage_1_percentage, codon_dict)

