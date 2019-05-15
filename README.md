# Pysam Iterators

This module contains functions and iterators which help in using Pysam

## Functionality

* Iteration over reads belonging to the same fragment (R1, R2)
* Caching of pysam.FastaFile reference sequences to memory for faster access
* get_aligned_pairs using a supplied fasta file when the MD tag doesn't match.
* get_aligned_pairs with additional sequencer cycle

### Prerequisites

The only prerequisite is [PySAM](https://github.com/pysam-developers/pysam)

### Examples

Iteration over R1 and R2:
```
with  pysam.AlignmentFile('test.bam') as bam:
    for R1,R2 in MatePairIterator( bam ):
        pass
```
Iteration over R1 and R2, on chromosome 2: (Arguments are passed to pysam.AlignmentFile fetch)
```
with  pysam.AlignmentFile('test.bam') as bam:
    for R1,R2 in MatePairIterator( bam, contig='chr2' ):
        print(R1,R2)
```

Iteration over query and reference base, where the reference base is extracted from a fasta file.
Make sure the fasta file has been indexed (samtools faidx).

```
reference = pysamiterators.CachedFasta( pysam.FastaFile('test.fasta') )
with  pysam.AlignmentFile('test.bam') as bam:
    for R1,R2 in MatePairIterator( bam, contig='chr2' ):
        break
    for read_index, reference_pos, reference_base in ReferenceBackedGetAlignedPairs(R1,
        reference=reference,
        matches_only=True,
        with_seq=True ):
        print(read_index, reference_pos, reference_base)
```
