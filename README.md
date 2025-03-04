# DIST - Distance Integration and Similarity Tool

A tool that compares gene expression similarities in bulk.

Using components from [Project Intronomicon](https://github.com/KorfLab/intronomicon) at [The Korf Lab](https://github.com/KorfLab).

Currently serving C. elegans only.

## Manifest

- `nameloader.pl` queries ensembl by gene name for aliases
- `name-resolver.py` makes index of aliases in json
- `genecount-extractor` prints gene id and count columns to screen
- `id-converter.py` converts aliases from expression files to WBGene ID
- `DIST.py` calculates distances between expression files in bulk, prints to screen
- `fakermaker.py` generates testing expression files

## Usage

```
perl nameloader.pl <download dir> <gene names>
```

```
python3 name-resolver.py --build file index.json
```

```
./genecount-extractor <tabular-file> --tsv <gene id col> <count col> > outfile.tsv
```

```
python3 id-converter.py <outfile.tsv> <index.json> -s
```

```
python3 DIST.py -i <directory of tsv files> -m/-c/-kl -s
```

## License

This project is licensed under the MIT License.
