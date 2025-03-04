#!/usr/bin/env python3

import argparse
import json
import os
import sys

def load_gene_mappings(json_file):
	with open(json_file, 'r') as fp:
		return json.load(fp)

def read_tsv(file):
	with open(file, 'r') as fp:
		lines = fp.readlines()
	return [line.strip().split('\t') for line in lines]

def write_tsv(data, output_file):
	with open(output_file, 'w') as fp:
		for line in data:
			fp.write('\t'.join(map(str, line)) + '\n')

def convert_gene_ids(data, gene_mappings, found_only):
	converted_data = []
	for gene_id, count in data:
		found = False
		for wb_id, aliases in gene_mappings.items():
			if gene_id in aliases or gene_id == wb_id:
				converted_data.append([wb_id, count])
				found = True
				break
		if not found and not found_only:
			converted_data.append([gene_id, count])
	return converted_data

parser = argparse.ArgumentParser(description='Convert gene IDs to WBGene IDs')
parser.add_argument('file', help='Input TSV file')
parser.add_argument('json', help='JSON file containing gene ID mappings')
parser.add_argument('-s', '--sort', action='store_true',
	help='Sort output by WBGene IDs')
parser.add_argument('-o', '--output',
	help='Output file (default: <input_file_basename>_WB.tsv)')
parser.add_argument('--found-only', action='store_true',
	help='Include only those gene IDs that have a found match')
args = parser.parse_args()

gene_mappings = load_gene_mappings(args.json)

data = read_tsv(args.file)

converted_data = convert_gene_ids(data, gene_mappings, args.found_only)

if args.sort: converted_data.sort(key=lambda x: x[0])

output_file = args.output if args.output else os.path.splitext(args.file)[0] + "_WB.tsv"
write_tsv(converted_data, output_file)
