from Bio.Seq import Seq

nucleotide_seq = """___""".replace("\n", "").replace(" ", "")

dna_seq = Seq(nucleotide_seq)
protein_seq = dna_seq.translate(to_stop=False)

def extract_region(codon_position, nucleotides_after):
    
    dna_start = (codon_position - 1) * 3  
    dna_end = dna_start + 3 + nucleotides_after
    return str(dna_seq[dna_start:dna_end])

region_25 = extract_region(25, 22)  #start position (by codons) here it is stated hat it is 25th amino-acid
region_100 = extract_region(100, 22)  #end position (by codons) here it is stated hat it is 100th amino-acid

pos_25_start = (25 - 1) * 3 + 1  
pos_100_start = (100 - 1) * 3 + 1

print("=== 25th Codon + 22nt Downstream ===")
print(f"DNA Position: {pos_25_start}-{pos_25_start + 3 + 22 - 1}")
print(f"Sequence (25nt total): {region_25}")
print(f"Starts with codon: {region_25[:3]} → {Seq(region_25[:3]).translate()}\n")

print("=== 100th Codon + 22nt Downstream ===")
print(f"DNA Position: {pos_100_start}-{pos_100_start + 3 + 22 - 1}")
print(f"Sequence (25nt total): {region_100}")
print(f"Starts with codon: {region_100[:3]} → {Seq(region_100[:3]).translate()}\n")

print("=== Suggested Forward/Reverse Primers ===")
print(f"Forward (25th codon region): 5'-{region_25[:20]}-3'")
print(f"Reverse (100th codon region): 5'-{Seq(region_100).reverse_complement()[:20]}-3'")

aa_start = 25
aa_end = 100
protein_fragment = protein_seq[aa_start-1:aa_end] 

dna_start = (aa_start - 1) * 3
dna_end = aa_end * 3  
dna_fragment = str(dna_seq[dna_start:dna_end])

print("\n=== 25th to 100th Amino Acid Fragment ===")
print(f"Amino Acid Range: {aa_start}-{aa_end}")
print(f"Protein Length: {len(protein_fragment)} residues")
print(f"Protein Sequence:\n{protein_fragment}")

print(f"\nCorresponding DNA Fragment ({dna_start+1}-{dna_end} bp):")
print(dna_fragment)

print(len(dna_fragment))
