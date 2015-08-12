# -*- coding: utf-8 -*-


"""
@package
@brief
@copyright
@author
"""

# Standard library imports
from gzip import open as gopen
import os

# Local imports
from FastqSeq import FastqSeq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ReadPair(object):
    """Combine 2 FastqSeq object"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__(self, R1, R2):
        """
        @param R1, R2 : read1, read2
        """
        assert R1.sampleName = R2.sampleName, 'R1sampleName:{}\tR2sampleName:{} '.format(R1.sampleName, R2.sampleName)
        assert R1.sampleIndex = R2.sampleIndex, 'R1sampleIndex:{}\tR2sampleIndex:{} '.format(R1.sampleIndex, R2.sampleIndex)
        assert R1.molecularIndex = R2.molecularIndex, 'R1molecularIndex:{}\tR2molecularIndex:{} '.format(R1.molecularIndex, R2.molecularIndex)

        self.R1 = R1
        self.R2 = R2
        self.molecularIndex = R1.molecularIndex





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



    def __call__(self):
        """ Simple fastq reader returning a generator over a fastq file """
        try:

			# Open the file depending of the compression status
            fastq = gopen(self.fastq_file, "rb") if self.fastq_file[-2:] == "gz" else open(self.fastq_file, "rb")


			# Iterate on the file until the end
            while True:

				# Extract informations from the fastq file
                name, seq, sep, qual = next(fastq), next(fastq), next(fastq), next(fastq)

				# Try to generate a valid FastqSeq object
                try:
                    yield FastqSeq(
					sampleName =  ":".join(name.split(":")[0:-2]),
					seq = seq.rstrip(),
					qual = qual.rstrip(),
					sampleIndex = name.split(":")[-2],
					molecularIndex = name.split(":")[-1])

                #self.n_seq += 1

                except AssertionError as E:
                    print(E)
                    print ("Skipping the sequence")

        except IOError as E:
            print(E)
            print ("Error while reading {} file".format(self.fastq_file))
            exit()

        except StopIteration:
            raise StopIteration("\t{} sequences parsed".format(self.n_seq))
            fastq.close()


 #~~~~~~~PROPERTIES~~~~~~~#

    @property
    def score(self ):
        """Return a score compute with lenght and quality of sequence"""
        s = (len(R1.seq)*R1.qual.mean())+(len(R2.seq)*R2.qual.mean())
        return s

 #~~~~~~~MAGIC METHODS~~~~~~~#

    def __magicsort__ (self):
        """Support for len method"""
        return
