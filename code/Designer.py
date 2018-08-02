from termcolor import colored


class Designer():

    gRNA = ''
    cutsite_index = ''
    stop_codon = ''
    three_letter_stop = ''
    sequence = ''
    hom_arm_left = ''
    hom_arm_right = ''
    forward_primer = ''
    second_forward_primer = ''

    reverse_primer_one = ''
    reverse_primer_two = ''
    rev_primer_one_complement = ''
    second_primer_complement = ''

    def __init__(self, gRNA='', stop_codon='', forward_primer='', second_forward_primer='', rev_primer_one='', rev_primer_two='', file_name='', hom_arm_length=0):

        # minus one because you are going to be deleting the '|'
        self.cutsite_index = gRNA.index('|') - 1
        self.gRNA = gRNA.replace('|', '')
        self.stop_codon = stop_codon
        self.three_letter_stop = stop_codon[-3:]

        with open(file_name, 'r') as f:
            self.sequence = f.read().replace('\n', '')

        self.sequence = self.sequence.rstrip()

        self.forward_primer = forward_primer
        self.second_forward_primer = second_forward_primer
        self.reverse_primer_one = rev_primer_one
        self.reverse_primer_two = rev_primer_two
        self.rev_primer_one_complement = self.complement(rev_primer_one)
        self.second_primer_complement = self.complement(rev_primer_two)

        # need to add 1 because looking at point in between the cutsite
        gRNA_index = self.sequence.index(self.gRNA) + self.cutsite_index + 1

        self.stop_codon_index = self.sequence.index(self.stop_codon) + self.stop_codon.index(self.three_letter_stop)

        self.hom_arm_left = self.sequence[(gRNA_index - hom_arm_length): gRNA_index]
        self.hom_arm_right = self.sequence[gRNA_index: (gRNA_index + hom_arm_length)]

        self.validateSequences()

    def validateSequences(self):

        if self.hom_arm_left in self.sequence and self.hom_arm_right in self.sequence:
            print('homology arm validation passed')
        else:
            self.printError('homology arm validation failed')

        if self.stop_codon in self.sequence:
            print('stop codon validation passed')
        else:
            self.printError('stop codon validation failed')

        if self.forward_primer in self.sequence:
            print('forward primer validation passed')
        else:
            self.printError('forward primer validation failed')

        if self.second_forward_primer in self.sequence:
            print('second forward primer validation passed')
        else:
            self.printError('second forward primer validation failed')

        if self.rev_primer_one_complement in self.sequence:
            print('reverse primer one complement validation passed')
        else:
            self.printError('reverse primer one validation failed')

        if self.second_primer_complement in self.sequence:
            print('reverse primer two complement validation passed')
        else:
            self.printError('reverse primer two validation failed')

        if self.sequence[(self.sequence.index(self.gRNA) + len(self.gRNA) + 1):self.sequence.index(self.gRNA) + len(self.gRNA) + 3] != 'GG':

            print(colored('Warning: ', 'red') + 'gRNA site is not immediately upstream of canonical Cas9 PAM Site(NGG). Ignore this if you are using not using spCas9.')

        # testing the output string
        self.print_sequence(colored=lambda x, y: x, test=True)

    def writeHomArmToFile(self, save_file_left, save_file_right):
        with open(save_file_left, 'w') as f:
            f.write(self.hom_arm_left)

        with open(save_file_right, 'w') as f:
            f.write(self.hom_arm_right)

    def print_sequence(self, colored=colored, test=False, upstream_margin=100, downstream_margin=50):

        output_string = ''

        # creating indexes for important features
        fp_index = self.sequence.index(self.forward_primer)
        rp_index = self.sequence.index(self.rev_primer_one_complement)
        second_primer_index = self.sequence.index(self.second_forward_primer)
        second_primer_complement_index = self.sequence.index(self.second_primer_complement)
        gRNA_print_index = self.sequence.index(self.gRNA)

        # forward primer one printing
        output_string += self.sequence[fp_index - upstream_margin: fp_index]
        output_string += colored(self.forward_primer, 'blue')
        output_string += self.sequence[fp_index + len(self.forward_primer): rp_index]

        # reverse primer one printing
        output_string += colored(self.rev_primer_one_complement, 'yellow')
        output_string += self.sequence[rp_index + len(self.rev_primer_one_complement): gRNA_print_index]

        # gRNA printing
        output_string += colored(self.gRNA, 'green')

        # stop codon and second forward primer printing

        # If forward primer for downstream homology arm is overlapping with the stop codon
        if self.stop_codon_index > second_primer_index and self.stop_codon_index < second_primer_index + len(self.second_forward_primer):
            output_string += self.sequence[gRNA_print_index + len(self.gRNA): second_primer_index]
            three_stop_index = self.second_forward_primer.index(self.three_letter_stop)
            output_string += colored(self.second_forward_primer[: three_stop_index], 'blue')
            output_string += colored(self.three_letter_stop, 'red')
            output_string += colored(self.second_forward_primer[three_stop_index + len(self.three_letter_stop):], 'blue')
            output_string += self.sequence[second_primer_index + len(self.second_forward_primer): second_primer_complement_index]

        else:
            output_string += self.sequence[gRNA_print_index + len(self.gRNA):self.stop_codon_index]
            output_string += colored(self.sequence[self.stop_codon_index:self.stop_codon_index + 3], 'red')

            # printing second primer
            output_string += self.sequence[self.stop_codon_index + 3:second_primer_index]
            output_string += colored(self.sequence[second_primer_index: second_primer_index + len(self.second_forward_primer)], 'blue')
            output_string += self.sequence[second_primer_index + len(self.second_forward_primer): second_primer_complement_index]

        output_string += colored(self.second_primer_complement, 'yellow')
        output_string += self.sequence[second_primer_complement_index + len(self.second_primer_complement): second_primer_complement_index + len(self.second_primer_complement) + downstream_margin]

        if test:
            if output_string in self.sequence:
                print('output string validation passed')
            else:
                print('output string validation failed')
            return

        # getting homology arm lengths

        gRNA_index = self.sequence.index(self.gRNA) + self.cutsite_index + 1
        hom_arm_length_left = gRNA_index - fp_index

        hom_arm_length_right = second_primer_complement_index + len(self.second_primer_complement) - gRNA_index

        # gettign difference from cut site
        left_difference = gRNA_index - (rp_index + len(self.rev_primer_one_complement))
        right_difference = second_primer_index - gRNA_index

        # printing
        print()
        print()
        print('gRNA cut site: ' + self.gRNA[:self.cutsite_index + 1] + colored('|', 'cyan') + self.gRNA[self.cutsite_index + 1:])
        print('Left Homology Arm: {} bp -> {} bp from gRNA cut site'.format(hom_arm_length_left, left_difference))
        print('Right Homology Arm: {} bp -> {} bp from gRNA cut site'.format(hom_arm_length_right, right_difference))
        print('-----------------------------------------------------')
        print('forward primer = ' + colored('blue', 'blue'))
        print('reverse primer = ' + colored('yellow', 'yellow'))
        print('gRNA = ' + colored('green', 'green'))
        print('stop codon = ' + colored('red', 'red'))
        print(output_string)

    def complement(self, sequence, three_to_five=False):
        complements = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        # have to flip bases AND reverse since its on opposite strand
        if three_to_five:
            # don't flip it want it 3' -> 5'
            return ''.join([complements[x] for x in sequence])
        else:
            return ''.join([complements[x] for x in sequence])[:: -1]

    def printPrimersPlusGibson(self):
        forward_gibson_one = colored('TCCCCGACCTGCAGCCCAGCT', 'magenta')
        forward_gibson_two = colored('AGTTCTTCTGATTCGAACATC', 'magenta')
        rev_gibson_two = colored('TGGAGAGGACTTTCCAAG', 'magenta')

        # this is the codon immediately upstream of stop codon(have to complement since this is obtained from positive strand)
        NNN = self.complement(self.sequence[(self.stop_codon_index - 3):self.stop_codon_index])
        rev_gibson_one = colored('CCGGAACCTCCTCCGCTCCC' + NNN, 'magenta')

        print('\n------Primers With Gibson Sequences: ')
        print('Forward Primer One: ' + forward_gibson_one + self.forward_primer + "   (5' -> 3')")
        print('Reverse Primer One: ' + rev_gibson_one + self.reverse_primer_one + " (5' -> 3')")
        print('Forward Primer Two: ' + forward_gibson_two + self.second_forward_primer + " (5' -> 3')")
        print('Reverse Primer Two: ' + rev_gibson_two + self.reverse_primer_two + " \t(5' -> 3')")

    def printgRNAWithOverhang(self):
        top_overhang_five_prime_end = 'CACC'
        bottom_overhang_five_prime_end = 'CAAA'

        # have to add G to 5' end of top-strand and C to 3' end of bottom strand b/c that G-C pair is necessary for the U6 promoter(according to paper)
        print('\n------gRNAs with overhang and added G-C base pair: ')
        print('Top Strand: ' + colored(top_overhang_five_prime_end, 'magenta') + colored('G', 'cyan') + self.gRNA + "\t  (5' -> 3')")
        print('Bottom Strand:  ' + colored('C', 'cyan') + self.complement(self.gRNA, three_to_five=True) + colored(bottom_overhang_five_prime_end, 'magenta') + " (3' -> 5')")

    def printError(self, s):
        print(colored('Error: ', 'red') + s)


if __name__ == '__main__':
    print(colored('Error: ', 'red') + 'Please run main.py')

