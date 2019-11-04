#/usr/bin/env python3

import svgwrite
import argparse
import regex

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', '-': '#B3B3B3'}

def parseSiteFile(file, top_number=50):
	offtargets = []
	total_seq = 0
	with open(file, "r") as f:
		f.readline()
		for line in f:
			line = line.rstrip("\n")
			line_items = line.split("\t")
			offtargets_site = line_items[4]
			reads = line_items[7]
			coordinate = line_items[2]
			strand = line_items[1]
			gene = line_items[3]
			locus = coordinate + "(" + strand + ")"
			offtargets.append([offtargets_site, reads, locus, gene])
			total_seq += 1
		if total_seq <= top_number:
			pass
		else:
			offtargets = offtargets[:top_number]
		return offtargets

def realign(seq, target):
	regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
	pattern = ''
	for letter in seq:
		pattern += regex_extended[letter]
	regex_pattern = '(?b:' + pattern + '){s<=6,d<=1}'
	m = regex.search(regex_pattern, target)
	realign_target = target
	for k in m.fuzzy_changes[2]:
		realign_target = realign_target[:k+1] + "-" + realign_target[k+1:]
	return seq, realign_target

def visualization(infile, target_seq, outfile, title, top_number=50):
	offtargets = parseSiteFile(infile, top_number)

	dwg = svgwrite.Drawing(outfile + ".svg", profile="full", size=(u'100%', 400 + top_number*(box_size + 1)))

	if title is not None:
		# Define top and left margins
		x_offset = 50
		y_offset = 50
		dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Calibri"))
	else:
		x_offset = 50
		y_offset = 20

	# Draw ticks
	tick_locations = [1, len(target_seq)]
	tick_locations += range(len(target_seq) + 1)[::10][1:]
	tick_locations.sort()
	for x in tick_locations:
		dwg.add(dwg.text(str(x), insert=(x_offset + (x -1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Calibri"))

	# Draw reference sequence row
	for i, c in enumerate(target_seq):
		y = y_offset
		x = x_offset + i * box_size
		dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
		dwg.add(dwg.text(c, insert=(x + 3, y + box_size -3), fill='black', style="font-size:15px; font-family:Calibri"))
	
	dwg.add(dwg.text('Reads', insert=(box_size * (len(target_seq) + 1) + 50, y_offset + box_size - 3), style="font-size:15px; font-family:Calibri; font-weight:bold"))
	dwg.add(dwg.text('Gene', insert=(box_size * (len(target_seq) + 1) + 100, y_offset + box_size - 3), style="font-size:15px; font-family:Calibri; font-weight:bold"))
	dwg.add(dwg.text('Locus', insert=(box_size * (len(target_seq) + 1) + 200, y_offset + box_size - 3), style="font-size:15px; font-family:Calibri; font-weight:bold"))


	# Draw aligned sequence rows
	y_offset += 1 # leave some extra space after the reference row
	line_number = 0 # keep track of ploted sequences
	
	for j, (sequences, reads, locus, gene) in enumerate(offtargets):
		k = 0
		line_number += 1
		y = y_offset + line_number * box_size
		dwg.add(dwg.text(str(j+1), insert=(x_offset - 30, 2 * box_size + y - 3), fill='black', style="font-size:10px; font-family:Calibri"))
		if len(sequences) == len(target_seq):
			for i, (c, r) in enumerate(zip(sequences, target_seq)):
				x = x_offset + k * box_size
				if c == "-":
					if 0 < k < len(target_seq):
						dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
						dwg.add(dwg.text(c, insert=(x + 4.5, 2 * box_size + y - 4), fill="black", style="font-size:15px; font-family:Calibri"))
						k += 1
				elif c == r:
					dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Calibri"))
					k += 1
				elif r == "N":
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Calibri"))
					k += 1
				else:
					dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Calibri"))
					k += 1
		else:
			sequence, realign_target = realign(sequences, target_seq)
			for i, (c, r) in enumerate(zip(sequence, realign_target)):
				x = x_offset + k * box_size
				if r == "-":
					if 0 < k < len(realign_target):
						x = x_offset + (k - 0.25) * box_size
						dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
						dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
				elif c == r:
					dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Calibri"))
					k += 1
				elif r == "N":
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Calibri"))
					k += 1
				else:
					dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
					dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Calibri"))
					k += 1

		dwg.add(dwg.text(str(reads), insert=(box_size * (len(target_seq) + 1) + 50, y_offset + box_size * (line_number + 2) - 2), 
					fill='black', style="font-size:15px; font-family:Calibri"))
		dwg.add(dwg.text(gene, insert=(box_size * (len(target_seq) + 1) + 100, y_offset + box_size * (line_number + 2) - 2), 
					fill='black', style="font-size:15px; font-family:Calibri"))
		dwg.add(dwg.text(locus, insert=(box_size * (len(target_seq) + 1) + 200, y_offset + box_size * (line_number + 2) - 2), 
					fill='black', style="font-size:15px; font-family:Calibri"))


	dwg.save()

def main():
	parser = argparse.ArgumentParser(description='Plot visualization plots for re-aligned reads.')
	parser.add_argument("--off_target_file", "-i", help="file stored all off-target sites", required=True)
	parser.add_argument("--top_number", "-n", help="the number of display off-target sites, default: 50", default=50, type=int)
	parser.add_argument("--target_seq", "-s", help="target sequence, which could contain 'N' letter", required=True)
	parser.add_argument("--outfile", "-o", help="output file", required=True)
	parser.add_argument("--title", "-t", help="Plot title", required=True)

	args = parser.parse_args()

	visualization(args.off_target_file, args.target_seq, args.outfile, args.title, args.top_number)

if __name__ == "__main__":
	main()

