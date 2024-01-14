# ITERATIVE K MEANS LONGITUDINAL TRAJECTORY CALCULATIONS

### Example usage

```
usage: ./kml.R [-h] --input INPUT --metadata METADATA --nclusters NCLUSTERS
               --algorithm ALGORITHM --subset SUBSET --transformation
               TRANSFORMATION --level LEVEL

Command Line Argument Parsing in R

options:
  -h, --help            show this help message and exit
  --input INPUT         Input .rds with count matrices per taxa
  --metadata METADATA   Input metadata path
  --nclusters NCLUSTERS
                        Number of trajectories to generate per taxa.
  --algorithm ALGORITHM
                        Clustering algorithm to use. If you're unsure use
                        parWithEuclidean_rndm
  --subset SUBSET       Naming variable only. Is this data a particular
                        subset?
  --transformation TRANSFORMATION
                        Naming variable only. Is this data transformed?
  --level LEVEL         Naming variable only. What taxonomic level does the
                        data represent?
```
