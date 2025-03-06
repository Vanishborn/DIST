#!/usr/bin/env python3

import sys
import os
import argparse
import math

def read_file(file):
	data = {}
	with open(file) as fp:
		for line in fp:
			gene, val = line.split()
			data[gene] = float(val)
	return data

def manhattan(d1, d2):
	total = 0
	for gene in d1:
		if gene not in d2: continue
		total += abs(d1[gene] - d2[gene])
	return total

def cartesian(d1, d2):
	total = 0
	for gene in d1:
		if gene not in d2: continue
		total += (d1[gene] - d2[gene]) ** 2
	result = math.sqrt(total)
	return result

def normalize(d1, d2):
	shared_keys = set(d1.keys()).intersection(set(d2.keys()))
	shared_d1 = {gene: d1[gene] for gene in shared_keys}
	shared_d2 = {gene: d2[gene] for gene in shared_keys}
	total_d1 = sum(shared_d1.values())
	total_d2 = sum(shared_d2.values())
	p1 = {gene: val / total_d1 for gene, val in shared_d1.items()}
	p2 = {gene: val / total_d2 for gene, val in shared_d2.items()}
	return p1, p2

def kullback(d1, d2):
	p1, p2 = normalize(d1, d2)
	total = 0
	for gene in p1:
		ratio = (p1[gene]) / (p2[gene])
		divergence = p1[gene] * math.log(ratio)
		total += divergence
	return total

def get_files(input_dir):
	files = []
	for file in os.listdir(input_dir):
		if file.endswith(".tsv"):
			full_path = os.path.join(input_dir, file)
			files.append(full_path)
	return files

parser = argparse.ArgumentParser(
	description="Calculate distances between gene expression files.")
parser.add_argument("-i", "--input", required=True,
	help="Input directory containing TSV files.")
parser.add_argument("-m", "--manhattan", action="store_true",
	help="Use Manhattan distance (default).")
parser.add_argument("-c", "--cartesian", action="store_true",
	help="Use Cartesian (Euclidean) distance.")
parser.add_argument("-kl", "--kullback-leibler",
	action="store_true", help="Use Kullback-Leibler divergence.")
parser.add_argument("-s", "--sort", action="store_true",
	help="Sort distances from smallest to largest.")
args = parser.parse_args()

files = sorted(get_files(args.input))
if len(files) < 2:
	print("Error: At least two TSV files are required for distance calculation.", file=sys.stderr)
	sys.exit(1)

data = []
for file in files:
	file_data = read_file(file)
	data.append(file_data)

distances = []
for i in range(len(data)):
	for j in range(i + 1, len(data)):
		if args.cartesian:
			dist = cartesian(data[i], data[j])
		elif args.kullback_leibler:
			dist = kullback(data[i], data[j])
		else:
			dist = manhattan(data[i], data[j])
		distances.append((files[i], files[j], dist))

if args.sort:
	distances.sort(key=lambda x: x[2])

for file1, file2, dist in distances:
	print(file1, file2, dist)
