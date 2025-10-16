from typing import Dict, List

class DNASequence:
    """Represents a DNA sequence and provides basic operations that are commonly needed in bioinformatics"""

    # Valid DNA nucleotides
    NUCLEOTIDES = ['A', 'T', 'G', 'C']

    # Complement mapping for reverse-complement
    COMPLEMENT_MAP = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }

    # Genetic Codon Table
    GENETIC_CODE = {
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

    def __init__(self, sequence: str, name: str  = ""):
        """
        Initialize a DNA sequence

        Args:
            sequence: DNA sequence string (A, T, G, C only)
            name: Optional name for the sequence
        """
        self.sequence = sequence.upper().strip()
        self.name = name
        self._validate_sequence()

    def _validate_sequence(self):
        """Validate that the sequence contains only valid DNA nucleotides"""
        invalid_chars = set(self.sequence) - set(self.NUCLEOTIDES)
        if invalid_chars:
            raise ValueError(f"Invalid nucleotides found: {invalid_chars}."
                             f"Only {self.NUCLEOTIDES} are allowed.")

    def __str__(self) -> str:
        """String representation of the sequence"""
        if self.name:
            return f"{self.name}: {self.sequence}"
        return self.sequence
    
    def __len__(self) -> int:
        """Returns the length of the sequence"""
        return len(self.sequence)
    
    def __getitem__(self, key) -> 'DNASequence | str':
        """Allow slicing and indexing like a string"""
        if isinstance(key, slice):
            return DNASequence(self.sequence[key], f"{self.name}_slice")
        return self.sequence[key]

    def __eq__(self, other) -> bool:
        """Compare two DNA sequences"""
        if isinstance(other, DNASequence):
            return self.sequence == other.sequence
        return self.sequence == str(other)

    def complement(self) -> 'DNASequence':
        """Returns the complement of the given DNA sequence"""
        complement = ''.join(self.COMPLEMENT_MAP[base] for base in self.sequence)
        return DNASequence(complement, f"{self.name}_complement")

    def reverse_complement(self) -> 'DNASequence':
        """
        Returns the reverse complement of the DNA sequence
        
        This is fundamental in biology - the complementary strand runs in the opposite direction.
        """
        rev_complement = self.complement().sequence[::-1]
        return DNASequence(rev_complement, f"{self.name}_reverse_complement")

    def transcribe(self) -> str:
        """Transcribe DNA to RNA (replace T with U). Returns the RNA sequence as a str"""
        return self.sequence.replace('T', 'U')

    def translate(self, start_pos: int = 0) -> str:
        """
        Translate DNA sequence to amino-acid sequence.

        Args:
            start_pos: Position to start the translation (0, 1, or 2 for reading frames)

        Returns:
            Amino acid sequence as string
        """
        # First transcribe to RNA
        rna = self.transcribe()

        # Start from the specified position
        rna = rna[start_pos:]

        # Translate in triplets
        amino_acids = []
        for i in range(0, len(rna) - 2, 3):
            codon = rna[i:i+3]
            amino_acid = self.GENETIC_CODE.get(codon, '?')
            amino_acids.append(amino_acid)

            # Stop at stop codons
            if amino_acid == "*":
                break

        return "".join(amino_acids)

    def gc_content(self, as_percentage: bool = False) -> float:
        """
        Calculate GC content
        
        GC content is important for DNA stability. Higher GC content means more stable DNA due to stronger hydrogen bonding
        """
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        ratio =  gc_count / len(self.sequence)
        return ratio * 100 if as_percentage else ratio

    def nucleotide_composition(self) -> Dict[str, int]:
        """Returns count of each nucleotide"""
        return {
            'A': self.sequence.count('A'),
            'T': self.sequence.count('T'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C'),
        }

    def find_motif(self, motif: str) -> List[int]:
        """
        Find all occurrences of a motif in the sequence

        Args:
            motif: The motif to search for (e.g., 'ATG' for start codon)

        Returns:
            List of starting positions where motif is found
        """
        motif = motif.upper()
        positions = []
        start = 0

        while True:
            pos = self.sequence.find(motif, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1

        return positions
