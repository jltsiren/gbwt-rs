```
# GFA to GBZ; ensure the right order for the paths
gfa2gbwt --approx-jobs 1 example

# GBZ to GG/GBWT; GG in Simple-SDS format
gfa2gbwt --decompress-graph --simple-sds-graph example

# GFA to GBZ; chop to length 2
gfa2gbwt --max-node 2 translation

# GBZ to GG/GBWT; GG in Simple-SDS format
gfa2gbwt --decompress-graph --simple-sds-graph translation
```
