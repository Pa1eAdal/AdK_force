In each file, data consists of multiple blocks.
The first line of each block contains information about that block.
For example: dv0 force10 site148 [AMP]30,
where 'dv0' indicates P_{close} = 0.28, representing wild-type AdK;
'force10' indicates an applied tensile force of 10pN;
'site148' indicates that force is applied to the 148th and 177th residues;
'[AMP]30' indicates that the concentration of AMP is 30 μM.

In 'other_site.dat', 'site144' refers to the 42nd and 144th residues, while 'site177' refers to the 55th and 177th residues.
In 'dna.dat', 'lod' represents the length of DNA in base pairs (bp), e.g., lod40 means the dsDNA has 40 bps.

Below is a list showing the correspondence between dv values and P_{close}:
dv0 0.28
dv1 0.58
dv2 0.78
dv10 0.12
dv14 0.99
dv15 0.02
dv16 0.01
