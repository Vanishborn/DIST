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
		if gene in d2:
			total += abs(d1[gene] - d2[gene])
	return total

def cartesian(d1, d2):
	total = 0
	for gene in d1:
		if gene in d2:
			total += (d1[gene] - d2[gene]) ** 2
	result = math.sqrt(total)
	return result

def normalize(data):
	total = sum(data.values())
	probabilities = {gene: count / total for gene, count in data.items()}
	correction_factor = 1.0 / sum(probabilities.values())
	probabilities = {gene: prob * correction_factor for gene, prob in probabilities.items()}
	return probabilities

def kullback(d1, d2):
	p1 = normalize(d1)
	p2 = normalize(d2)
	total = 0
	for gene in p1:
		if gene in p2:
			ratio = (p1[gene]) / (p2[gene])
			divergence = p1[gene] * math.log(ratio)
			total += divergence
	return abs(total)

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
