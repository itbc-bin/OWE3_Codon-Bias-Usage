import random as r


def sequence():
    """Generates random sequences of amino acid, which have
    a length between 600 and 800.

    :return: list - with 100 random sequences
    """
    sequences = []
    amino_acid = "ACDEFGHIKLMNPQRSTVWY"
    for number in range(400):
        seq = "M"
        for position in range(r.randint(600, 800)):
            seq += r.choice(amino_acid)
        sequences.append(seq)

    return sequences


def fasta_gen(sequences):
    """creates a Fasta formatted file from the inputted sequences

    :return: fasta file
    """
    n = 1
    try:
        with open("Random_seq.txt", "x") as file:
            for sequence in sequences:
                file.write(">Random_Seq" + str(n) + "\n")
                n += 1
                for part in range(1, len(sequence), 61):
                    file.write((sequence[part:part + 61]) + "\n")
                file.write("\n")
    except FileExistsError:
        print("There is already a random sequence file in the directory,"
              " if you would like to make a new file please delete"
              " the old one. \nIf you want a second file with random"
              " sequences rename the first file or"
              " move it to another directory")
        quit()


if __name__ == '__main__':
    random = sequence()
    fasta_gen(random)

