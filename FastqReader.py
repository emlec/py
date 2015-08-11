# -*- coding: utf-8 -*-


"""
@package    Sekator
@brief      Contain an iterator function to read fastq files and return FastqSeq objects
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
from gzip import open as gopen
import os

# Local imports
from FastqSeq import FastqSeq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqReader (object):
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~FUNDAMENTAL METHODS~~~~~~~#

    def __init__ (self, fastq_file):
        """
        @param fastq_file
        Check if fastq_file is readable
        """
        #if not os.path.exists(fastq_file):
         #   raise FileNotFoundError('%s: file not found' % fastq_file)
        #elif not os.access(fastq_file, os.R_OK):
        #    raise FileNotReadableError('%s: file not readable' % fastq_file)
		
        self.n_seq = 0
        self.fastq_file = fastq_file
			
	
        
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



    def __call__():
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
