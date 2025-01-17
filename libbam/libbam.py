#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 17:41:12 2018

@author: Antony Holmes
"""

import re
import subprocess


class SamRead:
    def __init__(
        self,
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
        tags=[],
    ):
        self._qname = qname
        self._rname = rname
        self._flag = flag
        self._pos = pos
        self._mapq = mapq
        self._cigar = cigar
        self._rnext = rnext
        self._pnext = pnext
        self._tlen = tlen
        self._seq = seq
        self._qual = qual
        self._tags = tags

    @property
    def qname(self):
        return self._qname

    @property
    def flag(self):
        return self._flag

    @property
    def rname(self):
        """
        Return the reference name, usually the chromosome id

        Returns
        -------
        str
            Reference name
        """

        return self._rname

    @property
    def pos(self):
        """
        Returns the start position. Assume 1-based unless user modified.

        Returns
        -------
        int
            Start position of alignment.
        """

        return self._pos

    @property
    def mapq(self):
        """
        Returns the mapping quality

        Returns
        -------
        int
            Mapping quality
        """

        return self._mapq

    @property
    def cigar(self):
        """
        Returns the CIGAR

        Returns
        -------
        str
            CIGAR alignment
        """

        return self._cigar

    @property
    def rnext(self):
        return self._rnext

    @property
    def pnext(self):
        return self._pnext

    @property
    def tlen(self):
        return self._tlen

    @property
    def seq(self):
        return self._seq

    @property
    def qual(self):
        return self._qual

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
        return self._tags

    @tags.setter
    def tags(self, tags):
        self._tags = tags

    def __str__(self):
        """
        Returns a tab delimited SAM string.
        """

        return "\t".join(
            [
                self.qname,
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
                "\t".join(self.tags),
            ]
        )

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
        sam = sam.strip().split("\t")

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

    read = SamRead(
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags
    )

    return read


class BamReader:
    def __init__(self, bam, paired=False, samtools="samtools" ):
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

        self._bam = bam
        self._paired = paired
        self._samtools = samtools

    def header(self):
        """
        Return the BAM/SAM header

        Returns
        -------
        generator
            Each line of the header
        """

        cmd = [self._samtools, "view", "-H", self._bam]

        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

        for l in stdout:
            yield l.decode("utf-8").strip()

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

        if self._paired:
            cmd = [self._samtools, "view", "-f", "3", self._bam]
        else:
            cmd = [self._samtools, "view", "-F", "4", self._bam]

        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

        for l in stdout:
            tokens = l.decode("utf-8").strip().split("\t")

            read = parse_sam_read(tokens)

            yield read

        stdout.close()

    def chrs(self):
        """
        List the chromosomes in the file.
        """

        cmd = [self._samtools, "idxstat", self._bam]
 
        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

        chrs = []

        for l in stdout:
            tokens = l.decode("utf-8").strip().split("\t")
            chr = tokens[0]

            #if not chr.startswith("chr"):

            if re.match(r'^(chr)?(\d+|[XYxyMm])$', chr):
                chrs.append(chr)

        stdout.close()

        return chrs #list(sorted(chrs))

    def reads(self, loc: str=''):
        """
        Iterate over the reads on a particular location, e.g. chromosome in the bam file.

        Parameters
        ----------
        l : str
            A genomic location e.g. chr1:1-10.
        """

        if self._paired:
            # properly mapped pair, first in pair and we use
            # the pnext and tlen to get an idea of the fragment
            # length. We use the min of the start and pnext as
            # the read start and the abs of tlen to get the fragment
            # length
            cmd = [self._samtools, "view", "-f", "67", self._bam]
        else:
            # mapped reads
            cmd = [self._samtools, "view", "-F", "4", self._bam]

        if loc:
            cmd.append(loc)

        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

        for loc in stdout:
            tokens = loc.decode("utf-8").strip().split("\t")

            read = parse_sam_read(tokens)

            yield read

        stdout.close()

    def count_reads(self, loc:str):
        """
        Iterate over the reads on a particular genome in the bam file.

        Parameters
        ----------
        l : str
            A location.
        """

        if self._paired:
            cmd = [self._samtools, "view", "-c", "-f", "3", self._bam, loc]
        else:
            cmd = [self._samtools, "view", "-c", "-F", "4", self._bam, loc]

        # print(' '.join(cmd))
        stdout = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout

        ret = int(stdout.readline().decode("utf-8").strip())

        stdout.close()

        return ret


class BamWriter:
    def __init__(self, bam, paired=False, samtools="samtools"):
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

        self._bam = bam

        # bam file to write to
        self._out = open(bam, "wb")

        # Maintain a pipe to output a sam read and write as bam ('-F', '4',)
        # Note the last '-' which samtools uses to get input from stdin
        self._stdin = subprocess.Popen(
            [samtools, "view", "-Sb", "-"], stdin=subprocess.PIPE, stdout=self._out
        ).stdin

    def _write(self, text):
        """
        Internal method for writing to BAM file. Do not call this
        method directly.

        Parameters
        ----------
        text : str
            Text to write to BAM.
        """

        self._stdin.write("{}\n".format(text).encode("utf-8"))

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

        self._out.close()
        self._stdin.close()
