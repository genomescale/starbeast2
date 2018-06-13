import sys
import os

def read_fasta(path):
	fasta_file = open(path)
	label = ""
	sequence = ""
	fasta_sequences = {}
	l = fasta_file.readline()
	while l != "":
		if l[0] == ">":
			if label != "":
				fasta_sequences[label] = sequence

			label = l[1:].rstrip()
			sequence = ""
		else:
			sequence += l.strip()

		l = fasta_file.readline()

	if label != "":
		fasta_sequences[label] = sequence

	return fasta_sequences

def write_fasta(path, fasta_sequences):
	output_file = open(path, "w")
	for label in sorted(fasta_sequences):
		sequence_block = ">" + label + "\n" + fasta_sequences[label] + "\n"
		output_file.write(sequence_block)
	output_file.close()

caninae_genera = {
	"Canis",
	"Cuon",
	"Lycaon"
}

caninae_folder = sys.argv[1]
caninae_filenames = os.listdir(caninae_folder)

for fn in caninae_filenames:
	if fn.endswith(".fasta"):
		fasta_path = os.path.join(caninae_folder, fn)
		caninae_sequences = read_fasta(fasta_path)
		canis_sequences = {}
		for label in caninae_sequences:
			genus, species, strand = label.split()[0].split("_")
			if genus in caninae_genera:
				canis_sequences[label] = caninae_sequences[label]

		write_fasta(fn, canis_sequences)
