import os
from Designer import Designer

os.system('clear')

# must include where gRNA cut site is with '|' symbol
gRNA = 'TAATTTCAATAGCATCC|GAA'
stop_codon = 'TGTGTAG'
file_name = 'ZFHX4_plus_1000.txt'
forward_primer_one = 'AGCACTTAGCCAACACCTCC'
reverse_primer_one = 'CCATCTAGGCCAGAAGCAGG'
forward_primer_two = 'GTGTGTAGGAGTGAAGACAGGA'
reverse_primer_two = 'TCCAGATTCCCTCTGCGGTA'

# ~750 HOMOLOGY ARM DESIGN
sequenceDesigner = Designer(gRNA=gRNA, stop_codon=stop_codon, forward_primer=forward_primer_one, second_forward_primer=forward_primer_two, rev_primer_one=reverse_primer_one, rev_primer_two=reverse_primer_two, file_name=file_name, hom_arm_length=750)
sequenceDesigner.writeHomArmToFile('hom_arm_left.txt', 'hom_arm_right.txt')
sequenceDesigner.print_sequence()
sequenceDesigner.printPrimersPlusGibson()
sequenceDesigner.printgRNAWithOverhang()
