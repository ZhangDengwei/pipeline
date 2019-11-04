#/usr/bin/env python3

import argparse
import pyfaidx
import regex
import HTSeq
import sys
import time

def reverseComplement(seq):
	compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
	out_list = [compl[bp] for bp in seq]
	return ''.join(out_list[::-1])

""" Get sequences from some reference genome
"""
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    # In case start position is samller than one or end position is larger than total length of chromosome
    start_position  = int(start) if int(start) > 1 else 1
    end_position = int(end) if int(end) < int(chrom_length[chromosome]) else int(chrom_length[chromosome])
    if strand == "+":
        seq = reference_genome[chromosome][start_position:end_position]
    elif strand == "-":
        seq = reference_genome[chromosome][start_position:end_position].reverse.complement
    return str(seq)

def regexFromSequence(seq, indels=1, errors=7):
    seq = seq.upper()
    """
    Given a sequence with ambiguous base characters, returns a regex that matches for
    the explicit (unambiguous) base characters
    """
    IUPAC_notation_regex = {'N': '[ATCGN]',
                            'Y': '[CTY]',
                            'R': '[AGR]',
                            'W': '[ATW]',
                            'S': '[CGS]',
                            'A': 'A',
                            'T': 'T',
                            'C': 'C',
                            'G': 'G'}
    pattern = ''
    for c in seq:
        pattern += IUPAC_notation_regex[c]

    pattern_standard = '(' + '?b:' + pattern + ')' + '{{s<={0}}}'.format(errors)
    pattern_gap = '(' + '?b:' + pattern + ')' + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
    return pattern_standard, pattern_gap

def alignDeletion(targetsite, deletion_seq, indels=1, errors=7):
	pattern_standard, pattern_gap = regexFromSequence(targetsite, indels, errors)
	m = regex.search(pattern_gap, deletion_seq)
	seq = deletion_seq
	for x in m.fuzzy_changes[2]:
		seq = seq[:x+1] + '-' + seq[x+1:]
	return seq

def alignSequences(targetsite, window_sequence, max_score=7):
	window_sequence = window_sequence.upper()
	query_regex_standard, query_regex_gap = regexFromSequence(targetsite, errors=max_score)

	# Try both strands
	alignments_mm, alignments_bulge = [], []
	alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
	alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
	alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
	alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

	# iterate through alignments_mm at first, then check alignments_bulge to ensure that a bulge is worse than a mismatch.
	lowest_distance, lowest_mismatch = 10, max_score
	final_alignment_seq, final_alignment_strand, final_alignment_start, final_alignment_end, mutation_type, distance_score = "", "", ".", ".", ".", "."
	for aln_m in alignments_mm:
		strand_m, alignment_type_m, match_m = aln_m
		if match_m != None:
			substitution, insertion, deletion = match_m.fuzzy_counts
			if substitution <= lowest_mismatch:
				mutation_type = "S" + str(substitution) + "I" + str(insertion) + "D" + str(deletion)
				distance_score = substitution + (insertion + deletion) * 3
				final_alignment_seq = match_m.group()
				final_alignment_strand = strand_m
				final_alignment_start = match_m.start()
				final_alignment_end = match_m.end()
				lowest_mismatch = substitution

	if not final_alignment_seq:
		for aln_b in alignments_bulge:
			strand_b, alignment_type_b, match_b = aln_b
			if match_b != None:
				substitution, insertion, deletion = match_b.fuzzy_counts
				if insertion or deletion:
					distance = substitution + (insertion + deletion) * 3
					edistance = substitution + insertion + deletion
					if distance < lowest_distance and edistance <= lowest_mismatch:
						mutation_type = "S" + str(substitution) + "I" + str(insertion) + "D" + str(deletion)
						distance_score = substitution + (insertion + deletion) * 3
						seq = match_b.group()
						final_alignment_seq = alignDeletion(targetsite, seq, errors=max_score) if deletion else seq
						final_alignment_strand = strand_b
						final_alignment_start = match_b.start()
						final_alignment_end = match_b.end()

	return [final_alignment_seq, final_alignment_strand, final_alignment_start, final_alignment_end, mutation_type, distance_score]

def read_GTF(GTF_file):
	gtf_file = HTSeq.GFF_Reader(GTF_file)
	ga_gene = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	for feature in gtf_file:
		if feature.type == 'gene':
			start = feature.iv.start
			end = feature.iv.end
			chrom = feature.iv.chrom
			# getting gene name
			ga_gene[HTSeq.GenomicInterval(chrom, start, end)] = sorted(feature.attr.items())[1][1]

	return ga_gene

def parase_cigar(Cigar, start_position):
	cigar_letter = [x for x in regex.split('\d', Cigar) if x != '']
	cigar_numer = [x for x in regex.split('\D', Cigar) if x != '']
	end_position = int(start_position)
	for i, k in zip(cigar_letter, cigar_numer):
		if i == "M" or i == "D":
			end_position += int(k)
		else:
			pass
	# end_position is the end position on the reference genome for a particular read
	return end_position

def analysis_sam(SAMfile, gap_threshold, overlap_threshold, output_file):
	ga = HTSeq.GenomicArray("auto", stranded=False)
	linear_contaminate_read, circle_read, right_orientation_wrong_location = 0, 0, 0

	with open(SAMfile, "r") as fin, open(output_file, "w") as fo:
		reverse_strand_flag = ['81', '83', '145', '147']

		header = ["Read_name", "Chromosome", "first_read_end_position", "Second_read_start_position", "Gap"]
		print(*header, sep="\t", file=fo)

		for line in fin:
			if not regex.search('^@', line):
				items = line.rstrip("\n").split("\t")
				read_name = items[0]
				flag = items[1]
				chromosome = items[2]
				leftmost_pos = items[3]
				map_quality = items[4]
				cigar = items[5]
				mate_chromosome = items[6]
				mate_read_leftmost_pos = items[7]
				insert_size = items[8]

				# For circle read, the leftmost position of reverse strand must precede mate read.
				if flag in reverse_strand_flag:
					if mate_chromosome == "=" and insert_size != '0' and int(map_quality) >= 50:
						if int(leftmost_pos) > int(mate_read_leftmost_pos):
							linear_contaminate_read += 2
						else:
							end_position = parase_cigar(cigar, leftmost_pos) - 1
							gap = int(mate_read_leftmost_pos) - end_position
							if gap <= 0 and abs(gap) > overlap_threshold:
								right_orientation_wrong_location += 2
							elif gap <= 0 and abs(gap) < overlap_threshold:
								iv = HTSeq.GenomicInterval(chromosome, end_position - overlap_threshold, int(mate_read_leftmost_pos), "+")
								ga[iv] += 2
								print(read_name, chromosome, end_position, mate_read_leftmost_pos, gap, sep="\t", file=fo)
								circle_read += 2
							elif 0 < gap <= int(gap_threshold):
								iv = HTSeq.GenomicInterval(chromosome, end_position, int(mate_read_leftmost_pos), "+")
								ga[iv] += 2
								print(read_name, chromosome, end_position, mate_read_leftmost_pos, gap, sep="\t", file=fo)
								circle_read += 2
							else:
								circle_read += 2

	print("Pattern 1:", file=sys.stderr)
	print("-------->\n		<---------", file=sys.stderr)
	print("The counts of liner contaminate reads (count twice for each paired read): ", linear_contaminate_read, file=sys.stderr)
	print("Pattern 2:", file=sys.stderr)
	print("<--------\n		-------->", file=sys.stderr)
	print("The counts of ds-circle reads (count twice for each paired read): ", circle_read, file=sys.stderr)
	print("Pattern 3:", file=sys.stderr)
	print("<--------", file=sys.stderr)
	print("     -------->", file=sys.stderr)
	print("The counts of mapped reads with right orientation but wrong location (overlap between read 1 and read 2): ", right_orientation_wrong_location, file=sys.stderr)

	return ga

def output_offtarget(SAMfile, search_radius, reference_genome, GTF_file, targetsite, gap_threshold, overlap_threshold, mismatch_threshold, 
					  circle_reads_coordinates, offtarget_output, unmatch_output):
	ga = analysis_sam(SAMfile, gap_threshold, overlap_threshold, circle_reads_coordinates)
	ga_gene = read_GTF(GTF_file)

	match_dic = {}
	off_target_sites_number, total_broken_points = 0, 0
	for iv, value in ga.steps():
		if value:
			gene_name = [value for iv, value in ga_gene[iv].steps()][0]
			gene_name = gene_name if gene_name else '*'
			
			window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - search_radius, iv.end + search_radius)

			final_alignment_seq, final_alignment_strand, final_alignment_start, final_alignment_end, mutation_type, distance_score = \
				alignSequences(targetsite, window_sequence, max_score=mismatch_threshold)
			absolute_start, absolute_end, absolute_strand = "", "", ""
			if final_alignment_seq and final_alignment_strand == "+":
				absolute_start = iv.start - search_radius + int(final_alignment_start)
				absolute_end = iv.start - search_radius + int(final_alignment_end)
				absolute_strand = final_alignment_strand
				absolute_window_sequence = window_sequence
			elif final_alignment_seq and final_alignment_strand == "-":
				absolute_start = iv.end + search_radius - int(final_alignment_end)
				absolute_end = iv.end + search_radius - int(final_alignment_start)
				absolute_strand = final_alignment_strand
				absolute_window_sequence = reverseComplement(window_sequence)
			else:
				absolute_start = iv.start
				absolute_end = iv.end
				absolute_strand = "*"
				absolute_window_sequence = window_sequence

			name = iv.chrom + ":" + str(absolute_start) + "-" + str(absolute_end)

			if name not in match_dic.keys():
				match_dic[name] = [iv.chrom, final_alignment_seq, mutation_type, str(distance_score), absolute_strand, 
										gene_name, value, absolute_window_sequence]
			else:
				match_dic[name][-2] += value

	match_dic = sorted(match_dic.items(), key=lambda item:item[1][-2], reverse=True)
	with open(offtarget_output, "w") as fo:
		with open(unmatch_output, "w") as fu:
			header = ["Chromosome", "strand", "Coordinate", "Gene", "off-target site", "Mutation type", "Distance score", "Reads", "Window sequence"]
			print(*header, sep="\t", file=fo)
			print(*header, sep="\t", file=fu)

			for key, value in match_dic:
				if value[4] != "*":
					print(value[0], value[4], key, value[5], value[1], value[2], value[3], str(int(value[6])), value[7], sep="\t", file=fo)
					off_target_sites_number += 1
					total_broken_points += 1
				else:
					print(value[0], value[4], key, value[5], value[1], value[2], value[3], str(int(value[6])), value[7], sep="\t", file=fu)
					total_broken_points += 1

	print("Total broken points: ", total_broken_points, file=sys.stderr)
	print("Total off-target sites (exclude on-target site): ", off_target_sites_number - 1, file=sys.stderr)
	print("*" * 100, "\n", file=sys.stderr)

def main():
	parse = argparse.ArgumentParser(description="find off-target site from aligned SAM file.")

	parse.add_argument("--reference", "-r", help="Reference Genome Fasta", required=True)
	parse.add_argument("--GTF_file", "-g", help="GTF file", required=True)
	parse.add_argument("--SAM", "-S", help="aligned SAM file", required=True)
	parse.add_argument("--targetsite", "-t", help="Target site sequene", required=True)
	parse.add_argument("--search_radius", help="Search radius around the position window, default: 20", default=20, type=int)
	parse.add_argument("--gap", help="Gap threshold, default: 3", default=3, type=int)
	parse.add_argument("--overlap_threshold", help="overlap threshold between pair reads, default: 5", default=5, type=int)
	parse.add_argument("--mismatch_threshold", help="Maximun score threshold, default: 6", default=6, type=int)
	parse.add_argument("--coordinate_output", "-c", help="the file stored broken points", required=True)
	parse.add_argument("--offtarget_output", "-o", help="the file stored off-target sites' info", required=True)
	parse.add_argument("--unmatch_output", "-u", help="the file stored broken points without potential off-target sites", required=True)

	args = parse.parse_args()

	start_time = time.time()
	start_time_out = time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(start_time))

	print("Strating at: ", start_time_out, file=sys.stderr)
	print("*" * 100, "\n", file=sys.stderr)

	# store the length of each chromosome
	global chrom_length
	chrom_length = {}
	reference_genome = pyfaidx.Fasta(args.reference)
	for chrom in reference_genome.keys():
		chrom_length[chrom] = reference_genome[chrom][:].end

	output_offtarget(args.SAM, args.search_radius, reference_genome, args.GTF_file, args.targetsite, args.gap, args.overlap_threshold, 
						args.mismatch_threshold, args.coordinate_output, args.offtarget_output, args.unmatch_output)

	end_time = time.time()
	end_time_out = time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(end_time))
	print("Ending at: ", end_time_out, file=sys.stderr)

	interval_time = end_time - start_time
	m, s = divmod(interval_time, 60)
	h, m = divmod(m, 60)
	print("Running time consumption: %02d:%02d:%02d" % (h, m, s), file=sys.stderr)


if __name__ == "__main__": 
	main()

