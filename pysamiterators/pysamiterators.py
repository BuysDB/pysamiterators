#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import random
import itertools
import numpy as np
import collections
import functools

# Global:
complement = str.maketrans('ATCGN', 'TAGCN')

class ReferenceBackedGetAlignedPairs(object):
    """
    This is a function which works similar to pysam get_aligned_pairs but
    the difference is that the reference base is fetched from a supplied fasta
    file. This can be usefull when mapping to a masked genome or when the mapper
    generated an invalid MD tag.
    """
    def __init__(self, read, reference, matches_only=False,with_seq=True):
        """Initialise  The iterator.
        Argument(s):
        handle to pysam.AlignmentFile(),
        handle to pysam.FastaFile(), (Tip: wrap in CachedFasta() )
        """
        self.read = read
        self.matches_only= matches_only
        self.with_seq = with_seq
        self.reference = reference

    def __repr__(self):
        return(f'{self.read}, iterator' )

    def __iter__(self):
        self.iterator = iter(self.read.get_aligned_pairs(matches_only=self.matches_only, with_seq=self.with_seq))
        return self

    def __next__(self):
        if self.with_seq:
            readIndex, referencePos, referenceBaseByPysam = next(self.iterator)
            if referencePos is not None:
                referenceBase = self.reference.fetch(self.read.reference_name,referencePos,referencePos+1)
            return readIndex, referencePos, referenceBase
        else:
            readIndex, referencePos = next(self.iterator)
            return readIndex, referencePos


class MatePairIterator():
    """Fast iteration over R1, R2"""

    def __init__(self, handle,performProperPairCheck=True, **kwargs):
        """Initialise  The iterator.

        Args:
            handle: handle to pysam.AlignmentFile(),
            **kwargs, arguments passed to fetch()

        Example:
        for R1,R2 in MatePairIterator( pysam.AlignmentFile('test.bam'), contig='chr1' )
        """
        self.handle = handle
        self.iterator = self.handle.fetch(**kwargs)
        self.cachedR1s = {}
        self.cachedR2s = {}
        self.performProperPairCheck = performProperPairCheck
        self.prevChromosome = None

    def __iter__(self):
        """Exectuted upon generator initiation."""
        return(self)

    def clearCache(self):
        self.cachedR1s = {}
        self.cachedR2s = {}


    def __next__(self):
        foundR1R2 = False
        while not foundR1R2:
            rec = next(self.iterator)
            if not rec.is_secondary and not  rec.is_supplementary:
                if rec.is_proper_pair:
                    haveR1 = False
                    haveR2 = False
                    if rec.is_read1:
                        if rec.query_name in self.cachedR1s:
                            raise( ValueError("Collision"))
                        self.cachedR1s[rec.query_name] = rec
                        haveR1 = True
                    else:
                        if rec.query_name in self.cachedR2s:
                            raise( ValueError("Collision"))
                        self.cachedR2s[rec.query_name] = rec
                        haveR2 = True
                    if (haveR1 or rec.query_name in self.cachedR1s) and (haveR2 or rec.query_name in self.cachedR2s):
                        return(self.cachedR1s.pop(rec.query_name), self.cachedR2s.pop(rec.query_name))
                elif not self.performProperPairCheck:
                    if rec.is_read1:
                        return (rec, None)
                    else:
                        return (None, rec)




# For a mapped read pair it is very important to figure out which bases are actual genomic signal
# Al sequence behind and aligned to the random primer(s) cannot be trusted and should be masked out
# secondly all signal before the starting location of R1 cannot be trusted
# This function returns a lower and higher bound of the locations within the fragment that can be trusted
# ASCII art: (H is primer sequence)
#           R1 H------------------------>
#   <------------------HH R2

# Output: (E is emitted)
#           R1 HEEEEEEE----------------->
#   <------------------HH R2

def getPairGenomicLocations(R1,R2, R1PrimerLength=4, R2PrimerLength=6):


    if  R1 is None or R2 is None or R1.is_unmapped:
        # This is an annoying situation, we cannot determine what bases can be trusted
        raise ValueError('Genomic locations cannot be determined')
    if R2.is_unmapped:
        # This is an annoying situation, we cannot determine what bases can be trusted
        raise ValueError('Genomic locations cannot be determined')

    if R1.is_reverse==R2.is_reverse:
        raise ValueError('Fragment incorrectly mapped') # The fragment is not correctly mapped

    if not R1.is_reverse: # R1 is forward, R2 is reverse
        #           R1 H------------------------>
        #   <------------------HH R2
        start = R1.reference_start+R1PrimerLength
        end = R2.reference_end -  R2PrimerLength
    else:
        #           R2 HH------------------------>
        #   <------------------HH R1
        start = R2.reference_start+R2PrimerLength
        end = R1.reference_end -  R1PrimerLength

    if start>=end:
        raise ValueError('Fragment has no size')

    return start,end


class JumpyMatePairIterator:
    """Fast iteration over R1, R2"""

    def __init__(self, handle, **kwargs):
        """Initialise  The iterator.

        Argument(s):
        handle to pysam.AlignmentFile(),
        Keyword arguments: fetch() arguments
        Example:
        for R1,R2 in MatePairIterator( pysam.AlignmentFile('test.bam'), contig='chr1' )
        """
        self.handle = handle
        self.iterator = self.handle.fetch(**kwargs)
        self.cachedR1s = {}
        self.cachedR2s = {}

    def __iter__(self):
        """Exectuted upon generator initiation."""
        return(self)

    def __next__(self):
        foundR1R2 = False
        while not foundR1R2:
            rec = next(self.iterator)
            if rec.is_paired and rec.is_proper_pair and not rec.is_secondary and not  rec.is_supplementary: # is_proper_pair
                haveR1 = False
                haveR2 = False
                if rec.is_read1:
                    if rec.query_name in self.cachedR1s:
                        raise( ValueError("Collision"))
                    self.cachedR1s[rec.query_name] = rec
                    haveR1 = True

                    if rec.next_reference_name != rec.reference_name:
                        # multi chromosome jump
                        pass

                else:
                    if rec.query_name in self.cachedR2s:
                        raise( ValueError("Collision"))
                    self.cachedR2s[rec.query_name] = rec
                    haveR2 = True
                if (haveR1 or rec.query_name in self.cachedR1s) and (haveR2 or rec.query_name in self.cachedR2s):
                    return(self.cachedR1s.pop(rec.query_name), self.cachedR2s.pop(rec.query_name))


class ReadSource():

    def __init__(self, handle):
        self.fetchingChromsome = None
        self.handle = handle
        self.currRecord = None
        self.ejectedRecord=False # If  the record was already yielded
        self.depleted =False
        self.currPos = 0


    def nextRead(self, chrom, pos):
        if chrom!=self.fetchingChromsome:
            self.fetchingChromsome = chrom
            try:
                self.readFetchIterator = self.handle.fetch(self.fetchingChromsome)
            except StopIteration:
                self.depleted=True
                self.currPos=None
                raise
        try:
            self.currRecord = next(self.readFetchIterator)
        except StopIteration:
            self.depleted=True
            raise
        self.ejectedRecord=False
        self.currPos = self.currRecord.reference_start
        return( self.currRecord )


    def eject(self):
        self.ejectedRecord = True

    def nextReadBuffered(self,chrom, pos):

        # There is a record waiting:
        if self.ejectedRecord==False and self.currRecord is not None:
            if self.currRecord.reference_start<=pos :
                self.eject()
                return self.currRecord
            else: # The record has to wait more
                return None
        # We need to get a new record:
        record = self.nextRead(chrom ,pos)
        if record.reference_start<=pos:
            self.eject()
            return self.currRecord




class SyncedIterator():
    """
    Yield reads from multiple bam files at the same location
    """
    def __init__(self, handles, contig):
        self.handles = [ReadSource(handle) for handle in handles]
        self.curpos = 0 # Current position on the contig
        self.contig=contig

    def __iter__(self):
        """Exectuted upon generator initiation."""
        return(self)

    def __next__(self):
        while True:

            for handle in self.handles:
                record  = handle.nextReadBuffered(self.contig, self.curpos)
                if record is not None:
                    return record

            # Obtain minimum next site:
            positions = [handle.currPos for handle in self.handles if handle.currPos is not None]
            if len(positions)==0:
                raise StopIteration
            self.curpos = min(positions)

class MixedIterator():
    """interleave reads  from the supplied bam files"""
    def __init__(self, handles, **kwargs):
        """Initialise  The iterator.

        Argument(s):
        handles to pysam.AlignmentFile():
        Kwargs: Keyword arguments: fetch() arguments
        """
        self.handles = handles
        # Yields: handleIndex, read
        self.iterator = ((i%len(self.handles),v) for i,v in enumerate(itertools.chain(*itertools.zip_longest( *self.handles ) )))
        self.depleted=False # if any of the iterators is depleted
        self.stopAtAnyDeplete = True

    def __iter__(self):
        """Exectuted upon generator initiation."""
        return(self)

    def __next__(self):
        while True:
            handleIndex, record = next(self.iterator)
            if record is not None:
                break
            else:
                self.depleted=True
                if self.stopAtAnyDeplete:
                    raise StopIteration
        return( handleIndex, record )

# Obtain the range where the reads map to
# emits contig, range start, range end
def getListSpanningCoordinates(l):
    surfaceStart = None
    surfaceEnd = None
    contig = None
    for read in l:
        if contig is None and read.reference_name is not None:
            contig = read.reference_name
        if surfaceStart is None or read.reference_start<surfaceStart:
            surfaceStart=read.reference_start
        if surfaceEnd is None or read.reference_end>surfaceEnd:
            surfaceEnd=read.reference_end
    return contig,surfaceStart,surfaceEnd


class CachedFasta():
    """
    Wrapper around pysam.FastaFile, stores the content of one or more chromosomes
    into ram to allow for faster retrieval
    """
    def __init__(self, handle, maximumCached=1):
        self.handle = handle
        self.maximumCached=maximumCached
        self.references = handle.references
        self.cachedSequences = {}
    def flush(self):
        self.cachedSequences = {}
    def fetch(self, chromosome=None, start=None, end=None):

        if not chromosome in self.cachedSequences:
            if len(self.cachedSequences)>=(self.maximumCached-1):
                self.cachedSequences={}

            self.cachedSequences[chromosome] = self.handle.fetch(chromosome)
        return(self.cachedSequences[chromosome][start:end])

class CachedFastaRange():
    def __init__(self, handle):
        self.handle = handle

    @functools.lru_cache(maxsize=10, typed=False)
    def fetch(self, *args):
        return self.handle.fetch(*args)
