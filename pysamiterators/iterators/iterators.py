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

def PileIterator(pile):
    """Obtain all read pairs overlapping with the supplied pileupcolumn

    Args:
        pile(pysam.pileupcolumn) : Pileup column to extract read pairs from.

    Yields:
        R1,R2 (list) : Read pairs overlapping the column

    Example:
        alignments =pysam.AlignmentFile(bamPath)
        for pc in alignments.pileup(chromosome,pos,pos+1 ):
            for R1,R2 in PileIterator(pc.pileups):
                print(R1,R2)
    """
    pair_dict = collections.defaultdict(dict) # name->True/False->read
    for pileupRead in pile:
        rec = pileupRead.alignment
        if not rec.is_secondary and not rec.is_supplementary:
            if not rec.is_proper_pair:
                if rec.is_read1:
                    yield (rec, None)
                else:
                    yield (None, rec)
                continue

            pair_dict[rec.query_name][rec.is_read1] = rec
            if True in pair_dict[rec.query_name] and False in pair_dict[rec.query_name]:
                R1,R2 = pair_dict[rec.query_name][True], pair_dict[rec.query_name][False]
                yield R1,R2
                del pair_dict[rec.query_name]
    # Obtain all non-matched reads:
    for query_name, singleton_read in pair_dict.items():
        is_read1 =  list(singleton_read.keys())[0]
        rec = list(singleton_read.values())[0]
        mate = alignments.mate( x.alignment )
        pair_dict[query_name][mate.is_read1] = mate

        R1,R2 = pair_dict[rec.query_name][True], pair_dict[rec.query_name][False]
        yield R1,R2

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
        self.iterator = iter(self.read.get_aligned_pairs(
            matches_only=self.matches_only, with_seq=self.with_seq))
        return self

    def __next__(self):
        if self.with_seq:
            readIndex, referencePos, referenceBaseByPysam = next(self.iterator)
            referenceBase=None
            if referencePos is not None:
                referenceBase = self.reference.fetch(self.read.reference_name,
                    referencePos,referencePos+1)
            return readIndex, referencePos, referenceBase
        else:
            readIndex, referencePos = next(self.iterator)
            return readIndex, referencePos


class MatePairIterator():
    """Fast iteration over R1, R2"""

    def __init__(self, handle,performProperPairCheck=True, ignore_collisions=False, **kwargs):
        """Initialise  The iterator.

        Args:
            handle: handle to pysam.AlignmentFile()

            ignore_collisions(bool) : when a read name is present multiple times for a single mate, do not throw an exception

            **kwargs, arguments passed to fetch()

        Example:
        for R1,R2 in MatePairIterator( pysam.AlignmentFile('test.bam'), contig='chr1' )
        """
        self.handle = handle
        if 'contig' in kwargs:
            self.iterator = self.handle.fetch(**kwargs)
        else:
            self.iterator = iter(self.handle)
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

            try:
                rec = next(self.iterator)
            except StopIteration:
                if not self.performProperPairCheck:
                    # yield remains in buffer..
                    if len(self.cachedR1s):
                        for query_name in self.cachedR1s:
                            return (self.cachedR1s.pop(query_name),None)
                    if len(self.cachedR2s):
                        for query_name in self.cachedR2s:
                            return (None,self.cachedR2s.pop(query_name))
                raise

            if not rec.is_secondary and not  rec.is_supplementary:
                # Check if mates are on same chromsome
                if (not rec.mate_is_unmapped and rec.reference_name ==rec.next_reference_name):
                    haveR1 = False
                    haveR2 = False
                    if rec.is_read1:
                        if rec.query_name in self.cachedR1s:
                            raise( ValueError(f"Collision {rec.query_name} is present multiple times as R1"))
                        self.cachedR1s[rec.query_name] = rec
                        haveR1 = True
                    else:
                        if rec.query_name in self.cachedR2s:
                            raise( ValueError(f"Collision  {rec.query_name} is present multiple times as R2"))
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

def getPairGenomicLocations(R1,R2, R1PrimerLength=4, R2PrimerLength=6, allow_unsafe=False):

    if R1 is None and allow_unsafe:
        if R2.is_reverse: # R1 is forward, R2 is reverse
            #           ?????????????????????
            #   <------------------HH R2
            start = R2.reference_start+R1PrimerLength
            end = R2.reference_end -  R2PrimerLength
        else:
            #           R2 HH------------------------>
            #   <------------------HH R1
            start = R2.reference_start+R2PrimerLength
            end = R2.reference_end -  R1PrimerLength
        return start,end

    if  (R1 is None and not allow_unsafe) or R2 is None or R1.is_unmapped:
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



def getCycleOffset(read, trimmed_begin_tag_R1='eB', trimmed_begin_tag_R2='EB'):
    if read is None:
        return None
    start=0
    if read.is_read1:
        start = read.get_tag(trimmed_begin_tag_R1) if read.has_tag(trimmed_begin_tag_R1) else 0
    elif read.is_read2:
        start = read.get_tag(trimmed_begin_tag_R2) if read.has_tag(trimmed_begin_tag_R2) else 0
    else:
        raise ValueError('Designed for single or mate pair only')
    return (start)

"""Obtain the total amount of cycles for the read,
    including bases which have been trimmed off from the start
"""
def getReadTotalCycles(read, cycleOffset=None,trimmed_begin_tag_R1='eB', trimmed_begin_tag_R2='EB'):
    # The obvious part:
    totalCycles = read.infer_read_length()
    # Add trimmed cycles:
    if cycleOffset is None:
        cycleOffset  = getCycleOffset(read, trimmed_begin_tag_R1, trimmed_begin_tag_R2)
    totalCycles += cycleOffset #@warn: end is not defined!
    return totalCycles



class ReadCycleIterator():
    """ This iterator is similar and a wrapper of the Pysam get_aligned_pairs function
    The difference is that the cycle of the sequencer is emitted (distToFirstCycle) (int or float)
    yields cycle, queryIndex, referencePos, (refbase)
    The cycle of the sequencer is obtained by the index, the read orientation and two tags which
    store how many bases have been trimmed of from the beginning of R1 and R2. (eB and EB by default)
    The second added feature is that a reference handle can be added which will yield
    reference bases from the supplied fasta file. This feature is neccesary when mapping to masked genomes
    """

    def __init__(self,
        read,
        matches_only=False, # transfered to pysam api
        with_seq=False, # emit reference bases
        emitFloats=False,  # Emit as percentage of total cycles instead of absolute cycles
        reference=None, # obtain reference base from a reference (Should be type pysam FastaFile)
        trimmed_begin_tag_R1='eB',
        trimmed_begin_tag_R2='EB'
        ):
        self.read = read
        self.trimmed_begin_tag_R1 = trimmed_begin_tag_R1
        self.trimmed_begin_tag_R2 = trimmed_begin_tag_R2
        self.with_seq = with_seq
        self.matches_only = matches_only
        self.emitFloats = emitFloats
        self.start = getCycleOffset(read,
                                    trimmed_begin_tag_R1=trimmed_begin_tag_R1,
                                    trimmed_begin_tag_R2=trimmed_begin_tag_R2)
        self.len = getReadTotalCycles(read,
                                      cycleOffset=self.start,
                                      trimmed_begin_tag_R1=trimmed_begin_tag_R1,
                                      trimmed_begin_tag_R2=trimmed_begin_tag_R2)
        self.reference = reference

    def __repr__(self):
        return(f'{self.read}, {self.len} cycles starting at {self.start}' )

    def __iter__(self):
        if self.reference is None:
            self.iterator = iter(self.read.get_aligned_pairs(matches_only=self.matches_only,
            with_seq=self.with_seq))
        else:
            self.iterator = iter(
            ReferenceBackedGetAlignedPairs(self.read,
                                           self.reference,
                                           matches_only=self.matches_only,
                                           with_seq=True))
        return self

    def __next__(self):
        if self.with_seq:
            readIndex, referencePos, referenceBase = next(self.iterator)
        else:
            readIndex, referencePos = next(self.iterator)
        # Obtain cycle:
        if readIndex is None:
            cycle=None
        else:
            if not self.read.is_reverse:
                cycle = readIndex+self.start
            else:
                #print(self.len, readIndex, self.start )

                cycle = self.len - readIndex - self.start -1 # minus one as the index starts at 0

            if self.emitFloats:
                if self.len==0:
                    cycle=0
                else:
                    cycle /= self.len
        if self.with_seq:
            return cycle, readIndex, referencePos, referenceBase
        else:
            return cycle, readIndex, referencePos



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
        self.iterator = ((i%len(self.handles),v) for i,v in enumerate(
            itertools.chain(*itertools.zip_longest( *self.handles ) )))
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
