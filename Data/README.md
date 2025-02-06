For readability, please refer to the 'readable_' files.
Each raw data file contains multiple data blocks.

The first line of each block contains metadata describing that block.
Example: dv0 force10 site148 [AMP]30

1. 'dv0' corresponds to P_{close} = 0.28, representing wild-type AdK
2. 'force10' indicates an applied tensile force of 10 pN
3. 'site148' specifies force application at residues 148 and 177
4. '[AMP]30' denotes an AMP concentration of 30 μM

In 'other_site.dat':

1. 'site144' refers to residues 42 and 144
2. 'site177' refers to residues 55 and 177

In 'dna.dat':

'lod' represents DNA length in base pairs (bp)
Example: 'lod40' indicates a 40 bp dsDNA

Mapping between dv values and P_{close}:
dv0 0.28
dv1 0.58
dv2 0.78
dv10 0.12
dv14 0.99
dv15 0.02
dv16 0.01