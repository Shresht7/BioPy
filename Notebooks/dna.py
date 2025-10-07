# Library
import random

# The nucleotide bases present in the DNA
BASES = ['A', 'T', 'G', 'C']

# Maps base short-code to full-names
BASE_NAMES = {
    "A": "Adenine",
    "T": "Thymine",
    "G": "Guanine",
    "C": "Cytosine"
}

# Maps nucleotide bases to their complements
BASE_COMPLEMENTS = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}

def generate_random_sequence(length: int = 12) -> str:
    """Generates a random DNA sequence of given length"""
    seq = [random.choice(BASES) for _ in range(length)]
    return "".join(seq)

def count_nucleotides(seq: str) -> dict[str, int]:
    """Count the number of nucleotides in the given DNA sequence"""
    counts = { base: seq.count(base) for base in BASES }
    return counts

def is_valid(seq: str) -> bool:
    """Checks whether the given sequence is a valid DNA sequence or not"""
    return all(base in BASES for base in seq)

def complement(seq: str) -> str:
    """Generates the complement for the given DNA sequence"""
    complement = [BASE_COMPLEMENTS[base] for base in seq]
    return "".join(complement)

def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of the given DNA sequence"""
    return complement(seq)[::-1]

def gc_content(seq: str, as_percentage: bool = False) -> float:
    """Calculate the ratio of GC content in the DNA sequence"""
    gc_count = seq.count('G') + seq.count('C')
    gc_content = gc_count / len(seq)
    return gc_content * 100 if as_percentage else gc_content

def visualize(seq: str, reverse: bool = False):
    """Prints a simple text-based visualization of the DNA strand"""
    top = seq
    connector = "".join("|" for _ in range(len(top)))
    bottom = reverse_complement(seq) if reverse else complement(seq)
    print(f"5' {top} 3'")
    print(f"   {connector}   ")
    print(f"3' {bottom} 5'")

# DNA Sequence Class
class DNASequence:
    def __init__(self, sequence: str):
        self.sequence = DNASequence.clean(sequence)
        self.validate() # Ensure that the DNA sequence is valid

    @staticmethod
    def clean(seq: str) -> str:
        seq = seq.replace("\n", "")
        seq = seq.replace(" ", "")
        seq = seq.upper()
        return seq

    @staticmethod
    def random(length = 12):
        """Generates a random DNA sequence of given length"""
        seq = [random.choice(BASES) for _ in range(length)]
        return DNASequence("".join(seq))

    def validate(self) -> bool:
        """Ensure that the DNA sequence is valid"""
        for base in self.sequence:
            if base not in BASES:
                raise ValueError(f"Invalid base: {base}")
        return True

    def count_nucleotides(self) -> dict[str, int]:
        """Counts the number of nucleotides in the DNA sequence"""
        return {base: self.sequence.count(base) for base in BASES}

    def count_bases(self) -> dict[str, int]:
        """Alias for `count_nucleotides`"""
        return self.count_nucleotides()

    def complement(self, reverse: bool = False):
        """Returns the complement or reverse-complement of the DNA sequence"""
        seq = self.sequence.translate(str.maketrans("ATGC", "TACG"))
        if reverse:
            seq = seq[::-1]
        return DNASequence(seq)

    def count(self, base: str) -> int:
        """Returns the count of the given base in the DNA sequence"""
        return self.sequence.count(base)

    def gc_content(self, as_percentage: bool = False) -> float:
        """Calculate the ratio of GC content in the DNA sequence"""
        gc_count = self.count('G') + self.count('C')
        gc_content = gc_count / len(self)
        return gc_content * 100 if as_percentage else gc_content

    def visualize(self, reverse: bool = False):
        """Prints a simple text-based visualization of the DNA strand"""
        top = self.sequence
        connector = "".join("|" for _ in range(len(top)))
        bottom = self.complement(reverse)
        print(f"5' {top} 3'")
        print(f"   {connector}   ")
        print(f"3' {bottom} 5'")

    def hamming_distance(self, other) -> int:
        """Compare two sequences of equal length"""
        if len(self.sequence) != len(other.sequence):
            raise ValueError("Sequences must be of equal length")
        return sum(a != b for a, b in zip(self.sequence, other.sequence))

    def transcribe(self):
        """Transcribes this DNA Sequence into RNA"""
        return self.sequence.replace('T', 'U')

    def translate(self):
        """Translates an mRNA sequence into a list of amino acids"""
        amino_acids = []

        mrna = self.transcribe()
        mrna = mrna.upper()

        # Find the start codon
        start_index = mrna.find("AUG")
        if start_index == -1:
            return [] # No start codon found.
        
        # Read codons in triplets from the start
        for i in range(start_index, len(mrna), 3):
            codon = mrna[i:i+3]
            if len(codon) < 3:
                break # Incomplete codon at the end?
            amino_acid = CODON_TABLE.get(codon, '???')
            if amino_acid == "*":
                break

            amino_acids.append(amino_acid)
        
        return "".join(amino_acids)

    def __len__(self):
        """Returns the length of the DNA Sequence"""
        return len(self.sequence)

    def __repr__(self):
        return f"<DNASequence(length={len(self.sequence)})?"

    def __str__(self):
        return self.sequence


CODON_TABLE = {
    # Phenylalanine (F)
    "UUU": "F", "UUC": "F",
    # Leucine (L)
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    # Isoleucine (I)
    "AUU": "I", "AUC": "I", "AUA": "I",
    # Methionine (M) - Start codon
    "AUG": "M",
    # Valine (V)
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    # Serine (S)
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    # Proline (P)
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine (T)
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine (A)
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine (Y)
    "UAU": "Y", "UAC": "Y",
    # Histidine (H)
    "CAU": "H", "CAC": "H",
    # Glutamine (Q)
    "CAA": "Q", "CAG": "Q",
    # Asparagine (N)
    "AAU": "N", "AAC": "N",
    # Lysine (K)
    "AAA": "K", "AAG": "K",
    # Aspartic Acid (D)
    "GAU": "D", "GAC": "D",
    # Glutamic Acid (E)
    "GAA": "E", "GAG": "E",
    # Cysteine (C)
    "UGU": "C", "UGC": "C",
    # Tryptophan (W)
    "UGG": "W",
    # Arginine (R)
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine (G)
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stop codons
    "UAA": "*", "UAG": "*", "UGA": "*"
}

def transcribe(sequence: str):
    """Transcribes this DNA Sequence into RNA"""
    return sequence.replace('T', 'U')

def translate(sequence: str):
    """Translates an mRNA sequence into a list of amino acids"""
    amino_acids = []

    mrna = transcribe(sequence)
    mrna = mrna.upper()

    # Find the start codon
    start_index = mrna.find("AUG")
    if start_index == -1:
        return [] # No start codon found.
    
    # Read codons in triplets from the start
    for i in range(start_index, len(mrna), 3):
        codon = mrna[i:i+3]
        if len(codon) < 3:
            break # Incomplete codon at the end?
        amino_acid = CODON_TABLE.get(codon, '???')
        if amino_acid == "*":
            break

        amino_acids.append(amino_acid)
    
    return "".join(amino_acids)
