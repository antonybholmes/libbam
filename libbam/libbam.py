#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 17:41:12 2018

@author: Antony Holmes
"""

import sys
import subprocess


class SamRead(object):
    def __init__(self,
                 qname,
                 flag,
                 rname,
                 pos,
                 mapq,
                 cigar,
                 rnext,
                 pnext,
                 tlen,
                 seq,
                 qual,
                 tags = []):
        self.__qname = qname
        self.__rname = rname
        self.__flag = flag
        self.__pos = pos
        self.__mapq = mapq
        self.__cigar = cigar
        self.__rnext = rnext
        self.__pnext = pnext
        self.__tlen = tlen
        self.__seq = seq
        self.__qual = qual
        self.__tags = tags
     
    @property
    def qname(self):
        return self.__qname

    @property
    def flag(self):
        return self.__flag
    
    @property
    def rname(self):
        """
        Return the reference name, usually the chromosome id
        
        Returns
        -------
        str
            Reference name
        """
        
        return self.__rname
    
    @property
    def pos(self):
        """
        Returns the start position. Assume 1-based unless user modified.
        
        Returns
        -------
        int
            Start position of alignment.
        """
        
        return self.__pos
    
    @property
    def mapq(self):
        """
        Returns the mapping quality
        
        Returns
        -------
        int
            Mapping quality
        """
        
        return self.__mapq
    
    @property
    def cigar(self):
        """
        Returns the CIGAR
        
        Returns
        -------
        str
            CIGAR alignment
        """
        
        return self.__cigar
    
    @property
    def rnext(self):
        return self.__rnext
    
    @property
    def pnext(self):
        return self.__pnext
    
    @property
    def tlen(self):
        return self.__tlen
    
    @property
    def seq(self):
        return self.__seq
    
    @property
    def qual(self):
        return self.__qual
    
    @property
    def is_paired(self):
        return self.flag & 2
    
    @property
    def is_proper_pair(self):
        return self.flag & 1
    
    @property
    def length(self):
        return len(self.seq)
    
    @property
    def tags(self):
        return self.__tags
    
    @tags.setter
    def tags(self, tags):
        self.__tags = tags
        
    def __str__(self):
        """
        Returns a tab delimited SAM string.
        """
        
        return '\t'.join([self.qname,
                          str(self.flag),
                          self.chr,
                          str(self.pos),
                          str(self.mapq),
                          self.cigar,
                          self.rnext,
                          str(self.pnext),
                          str(self.tlen),
                          self.seq,
                          self.qual,
                          '\t'.join(self.tags)])
        
    #
    # Alternative property names
    #
    
    @property
    def chr(self):
        """
        Alias for rname since this usually contains the chr
        
        Returns
        -------
        str
            Chromosome
        """
        
        return self.rname
    

def parse_sam_read(sam):
    """
    Parses a SAM alignment and returns a SamFile object.
    
    Parameters
    ----------
    sam : str or list
        Either a tab delimited SAM string, or an already tokenized 
        list of SAM fields.
    
    Returns
    -------
    SamRead
        a SamRead object representation of the SAM alignment.
    """
    
    if isinstance(sam, str):
        sam = sam.strip().split('\t')
        
    if not isinstance(sam, list):
        return None
    
    qname = sam[0]
    flag = int(sam[1])
    rname = sam[2]
    pos = int(sam[3])
    mapq = int(sam[4])
    cigar = sam[5]
    rnext = sam[6]
    pnext = int(sam[7])
    tlen = int(sam[8])
    seq = sam[9]
    qual = sam[10]
    
    tags = sam[11:]
    
    read = SamRead(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags)
    
    return read

    
class BamReader(object):
    def __init__(self, 
                 bam, 
                 paired=False, 
                 samtools='samtools',
                 args=[{'-F', '4'}]):
        """
        Create a new SAM reader
        
        Parameters
        ----------
        bam : str
            SAM/BAM file path
        samtools : str, optional
            Path to samtools executable. Default assumes it can be
            found on the sys path
        """
        
        self.__bam = bam
        self.__paired = paired
        self.__samtools = samtools
    
    def header(self):
        """
        Return the BAM/SAM header
        
        Returns
        -------
        generator
            Each line of the header
        """
        
        cmd = [self.__samtools, 'view', '-H', self.__bam]

        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

        for l in stdout:
            yield l.decode('utf-8').strip()

        stdout.close()
    
    def print_header(self):
        """
        Print the BAM/SAM header
        """
        
        for l in self.header():
            print(l)
        
    def __iter__(self):
        """
        Iterate over the reads in the bam file.
        """
        
        if self.__paired:
            cmd = [self.__samtools, 'view', '-f', '3', self.__bam]
        else:
            cmd = [self.__samtools, 'view', '-F', '4', self.__bam]
            
        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
        
        for l in stdout:
            tokens = l.decode('utf-8').strip().split('\t')
    
            read = parse_sam_read(tokens)
            
            yield read
            
        stdout.close()
        
    def reads(self, l):
        """
        Iterate over the reads on a particular genome in the bam file.
        
        Parameters
        ----------
        l : str
            A genomic location e.g. chr1:1-10.
        """
        
        if self.__paired:
            cmd = [self.__samtools, 'view', '-f', '3', self.__bam, l]
        else:
            cmd = [self.__samtools, 'view', '-F', '4', self.__bam, l]
            
        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
        
        for l in stdout:
            tokens = l.decode('utf-8').strip().split('\t')
    
            read = parse_sam_read(tokens)
            
            yield read
            
        stdout.close()
        
    def count_reads(self, l):
        """
        Iterate over the reads on a particular genome in the bam file.
        
        Parameters
        ----------
        l : str
            A location.
        """
        
        if self.__paired:
            cmd = [self.__samtools, 'view', '-c', '-f', '3', self.__bam, l]
        else:
            cmd = [self.__samtools, 'view', '-c', '-F', '4', self.__bam, l]
         
        #print(' '.join(cmd))
        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
        
        ret = int(stdout.readline().decode('utf-8').strip())
        
        stdout.close()
    
        return ret
            
        
        

class BamWriter(object):
    def __init__(self, bam, paired=False, samtools='samtools'):
        """
        Create a new SAM reader
        
        Parameters
        ----------
        bam : str
            SAM/BAM file path
        samtools : str, optional
            Path to samtools executable. Default assumes it can be
            found on the sys path
        """
        
        self.__bam = bam
        
        # bam file to write to
        self.__out = open(bam, 'wb')
        
        # Maintain a pipe to output a sam read and write as bam ('-F', '4',)
        # Note the last '-' which samtools uses to get input from stdin
        self.__stdin = subprocess.Popen([samtools, 'view', '-Sb', '-'], stdin=subprocess.PIPE, stdout=self.__out).stdin
    
    
    def _write(self, text):
        """
        Internal method for writing to BAM file. Do not call this
        method directly.
        
        Parameters
        ----------
        text : str
            Text to write to BAM.
        """
        
        self.__stdin.write('{}\n'.format(text).encode('utf-8'))
        
    
    def write_header(self, samreader):
        """
        Writes the header from a SamReader to the BAM file. Note that
        if you pipe the header as text from a SAM file into samtools
        view and output in BAM format, the header will be encoded as
        BAM and can be written to the beginning of the file. This does
        not work if you output SAM.
        
        Parameters
        ----------
        samreader : SamReader
            A SamReader object.
        """
        
        for line in samreader.header():
            self._write(line)
    
    def write(self, read):
        """
        Writes a SamRead to BAM.
        
        Parameters
        ----------
        read : SamRead
            A SamRead object.
        """
        
        self._write(str(read))
        
    def close(self):
        """
        Close the BAM file when finished writing.
        """
        
        self.__out.close()
        self.__stdin.close()
