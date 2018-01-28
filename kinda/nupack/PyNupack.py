# PyUtils.py
#
# For use with ChE/BE 163 problem set #2.
#
# Written by Justin Bois, May 2007.
#
# This file contains utility functions for to enable easy interfacing
# between Python and calls to the NUPACK core executables.  It additionally
# contains utility scripts for converting structures from dot-paren
# notation to pair lists.


# ############################################################### #
def getEnergy(sequence,structure,T,material,NUPACKHOME):

    # USAGE: energy = getEnergy(sequence,structure,T,material,NUPACKHOME)
    #
    # Runs the NUPACK executable energy to get the energy of a secondary
    # structure.  Parses all I/O to and from NUPACK.
    #
    # The inputs are:
    #   sequence: the base sequence of the strand
    #   structure: the structure in dot-parens format
    #   T: the temperature in degrees celsius
    #   material: a string, either 'dna', 'rna', or 'rna37'
    #    (see NUPACK User Guide)
    #   NUPACKHOME: a string with the name of the directory containing
    #     NUPACK.  There is NO trailing slash in this string.  E.g.,
    #     NUPACKHOME = '/home/justin/nupack2.0'
    #
    # The free energy of the secondary structure, in units of kcal/mol,
    # is returned.

    import os
    import re
    whiteSpaceSearch = re.compile('\s+')
    import random


    # Name of file for I/O
    fname = 'junk_file_delete_me%d' % random.randint(100000,999999)
    outputfile = '%s.out' % fname

    # Make the input file
    inputFileName = '%s.in' % fname
    f = open(inputFileName,'w')
    f.write('%s\n%s\n' % (sequence,structure))
    f.close()

    # Create the command to compute the energy
    cmd = '%s/bin/energy -multi -T %.1f -material %s %s > %s' % (NUPACKHOME,T,material,fname,outputfile)

    # Run the command in a subshell
    os.system(cmd)

    # Parse the output
    f = open(outputfile,'r')

    # Blow through comments
    line = f.readline()
    while line[0] == '%':
        line = f.readline()

    # Now we're at the energy, in kcal/mol
    lineData = whiteSpaceSearch.split(line)
    energy = float(lineData[0])
    
    f.close()

    # Remove junk files
    cmd = 'rm -f %s.in ; rm -f %s' % (fname,outputfile)
    os.system(cmd)
    
    return energy
# ############################################################### #



# ############################################################### #
def getProb(sequence,structure,T,material,NUPACKHOME):

    # USAGE: prob = getProb(sequence,structure,T,material,NUPACKHOME) 
    #
    # Runs the NUPACK executable prob to get the equilibrium probability
    # of a secondary structure.  Parses all I/O to and from NUPACK.
    #
    # The inputs are:
    #   sequence: the base sequence of the strand
    #   structure: the structure in dot-parens format
    #   T: the temperature in degrees celsius
    #   material: a string, either 'dna', 'rna', or 'rna37'
    #    (see NUPACK User Guide)
    #   NUPACKHOME: a string with the name of the directory containing
    #     NUPACK.  There is NO trailing slash in this string.  E.g.,
    #     NUPACKHOME = '/home/justin/nupack2.0'
    #
    # The equilibrium probability of the sequence having the secondary
    # structure is returned.

    import os
    import re
    whiteSpaceSearch = re.compile('\s+')
    import random

    # Name of file for I/O
    fname = 'junk_file_delete_me%d' % random.randint(100000,999999)
    outputfile = '%s.out' % fname

    # Make the input file
    inputFileName = '%s.in' % fname
    f = open(inputFileName,'w')
    f.write('%s\n%s\n' % (sequence,structure))
    f.close()

    # Create the command to compute the probability
    cmd = '%s/bin/prob -multi -T %.1f -material %s %s > %s' % (NUPACKHOME,T,material,fname,outputfile)

    # Run the command in a subshell
    os.system(cmd)

    # Parse the output
    f = open(outputfile,'r')

    # Blow through comments
    line = f.readline()
    while line[0] == '%':
        line = f.readline()

    # Now we're at the probability
    lineData = whiteSpaceSearch.split(line)
    print '"' + lineData[0] + '"'
    prob = float(lineData[0])
    
    f.close()

    # Remove junk files
    cmd = 'rm -f %s.in ; rm -f %s' % (fname,outputfile)
    os.system(cmd)
    
    return prob
# ############################################################### #


# ############################################################### #
def getMFEStruct(sequence,T,material,NUPACKHOME):

    # USAGE: mfeStruct = getMFEStruct(sequence,T,material,NUPACKHOME)
    #
    # Runs the NUPACK executable mfe to get the MFE structure for the
    # input sequence.  Parses all I/O to and from NUPACK.
    #
    # The inputs are:
    #   sequence: the base sequence of the strand
    #   T: the temperature in degrees celsius
    #   material: a string, either 'dna', 'rna', or 'rna37'
    #    (see NUPACK User Guide)
    #   NUPACKHOME: a string with the name of the directory containing
    #     NUPACK.  There is NO trailing slash in this string.  E.g.,
    #     NUPACKHOME = '/home/justin/nupack2.0'
    #
    # The minimal free energy structure in dot-paren notation is returned.


    import os
    import re
    whiteSpaceSearch = re.compile('\s+')
    import random


    # Name of file for I/O
    fname = 'junk_file_delete_me%d' % random.randint(100000,999999)
    outputfile = '%s.mfe' % fname

    # Make the input file
    inputFileName = '%s.in' % fname
    f = open(inputFileName,'w')
    f.write('%s\n' % sequence)
    f.close()

    # Create the command to compute the MFE struct
    cmd = '%s/bin/mfe -multi -T %.1f -material %s %s' % (NUPACKHOME,T,material,fname)

    # Run the command in a subshell
    os.system(cmd)

    # Parse the output
    f = open(outputfile,'r')

    # Blow through comments and blank lines
    line = f.readline()
    while line[0] == '%' or line[0] == '\n' or line[0] == '\0':
        line = f.readline()

    # Now we are at the entry containing number of bases
    line = f.readline()

    # Now we're at the free energy of the MFE structure
    line = f.readline()

    # Now we're at the dot paren structure
    lineData = whiteSpaceSearch.split(line)
    mfeStruct = lineData[0]
    
    f.close()
    
    # Remove junk files
    cmd = 'rm -f %s.in ; rm -f %s' % (fname,outputfile)
    os.system(cmd)

    return mfeStruct
# ############################################################### #

def getSubopt(sequences, T, material, dangles, energy_gap):
    # USAGE: getSubopt = getSubopt(sequence,T,material,energy_gap,NUPACKHOME)
    #
    # Runs the NUPACK executable subopt to get the MFE structures for the
    # input sequence.  Parses all I/O to and from NUPACK.
    #
    # The inputs are:
    #   sequence: the base sequence of the strand
    #   T: the temperature in degrees celsius
    #   material: a string, either 'dna', 'rna', or 'rna37'
    #    (see NUPACK User Guide)
    #   energy_gap: energy gap in which to get structures
    #   NUPACKHOME: a string with the name of the directory containing
    #     NUPACK.  There is NO trailing slash in this string.  E.g.,
    #     NUPACKHOME = '/home/justin/nupack2.0'
    #
    # The minimal free energy structure in dot-paren notation is returned.
    

    import os
    import re
    whiteSpaceSearch = re.compile('\s+')
    import random

    NUPACKHOME = os.environ['NUPACKHOME']

    # Name of file for I/O
    fname = 'junk_file_delete_me%d' % random.randint(100000,999999)
    outputfile = '%s.subopt' % fname

    # Make the input file
    inputFileName = '%s.in' % fname
    f = open(inputFileName,'w')
    inputContents = [str(len(sequences))]
    inputContents.extend(sequences)
    inputContents.append(" ".join([str(i+1) for i in range(len(sequences))]))
    inputContents.append(str(energy_gap))
    f.write("\n".join(inputContents))
    f.close()

    # Create the command to compute the MFE struct
    cmd = '%s/bin/subopt -multi -T %.1f -material %s -dangles %s %s' % (NUPACKHOME,T,material,dangles,fname)

    # Run the command in a subshell
    os.system(cmd)

    # Parse the output
    f = open(outputfile,'r')

    struct_list = []
    line = f.readline()
    while line[0] == '%' or line[0] == '\n' or line[0] == '\0':
        line = f.readline()
    while line != '':
      # Now we are at the entry containing number of bases
      line = f.readline()
      energy = float(line)

      # Now we're at the free energy of the MFE structure
      line = f.readline()

      # Now we're at the dot paren structure
      lineData = whiteSpaceSearch.split(line)
      mfeStruct = lineData[0]
      
      struct_list.append((mfeStruct, energy))
      
      # Blow through pair list notation
      while line != '' and line[0] != '%' and line[0] != '\n' and line[0] != '\0':
        line = f.readline()
      while line != '' and (line[0] == '%' or line[0] == '\n' or line[0] == '\0'):
          line = f.readline()
    
    f.close()
    
    # Remove junk files
    cmd = 'rm -f %s.in ; rm -f %s' % (fname,outputfile)
    os.system(cmd)

    return struct_list
# ############################################################### #

# ############################################################### #
def getPfunc(sequences, T, material, dangles):
    # USAGE: pfunc = getPfunc(sequence,T,material, dangles)
    #
    # Runs the NUPACK executable pfunc to get the partition function
    # input sequence.  Parses all I/O to and from NUPACK.
    #
    # The inputs are:
    #   sequence: the base sequence of the strand
    #   T: the temperature in degrees celsius
    #   material: a string, either 'dna', 'rna', or 'rna37'
    #    (see NUPACK User Guide)
    #
    

    import os
    import re
    whiteSpaceSearch = re.compile('\s+')
    import random

    NUPACKHOME = os.environ['NUPACKHOME']

    # Name of file for I/O
    fname = 'junk_file_delete_me%d' % random.randint(100000,999999)
    outputfile = '%s.pfunc' % fname

    # Make the input file
    inputFileName = '%s.in' % fname
    f = open(inputFileName,'w')
    inputContents = [str(len(sequences))]
    inputContents.extend(sequences)
    inputContents.append(" ".join([str(i+1) for i in range(len(sequences))]))
    f.write("\n".join(inputContents))
    f.close()

    # Create the command to compute the MFE struct
    cmd = '%s/bin/pfunc -multi -T %.1f -material %s -dangles %s %s' % (NUPACKHOME,T,material,dangles,fname)

    # Run the command in a subshell
    os.system(cmd)

    # Parse the output
    f = open(outputfile,'r')

    struct_list = []
    line = f.readline()
    while line[0] == '%' or line[0] == '\n' or line[0] == '\0':
      line = f.readline()
        
    line = f.readline()
    pfunc = f.readline()
    
    f.close()
    
    # Remove junk files
    cmd = 'rm -f %s.in ; rm -f %s' % (fname,outputfile)
    os.system(cmd)

    return pfunc
# ############################################################### #


# ############################################################### #
def getPairProbs(sequence,T,material,NUPACKHOME):

    # USAGE: (pairprobs,unpairprobs) = getPairProbs(sequence,T,material,NUPACKHOME)
    #
    # Runs the NUPACK executable pairs to get the pair probabilities
    # for the input sequence.  Parses all I/O to and from NUPACK.
    #
    # The inputs are:
    #   sequence: the base sequence of the strand
    #   T: the temperature in degrees celsius
    #   material: a string, either 'dna', 'rna', or 'rna37'
    #    (see NUPACK User Guide)
    #   NUPACKHOME: a string with the name of the directory containing
    #     NUPACK.  There is NO trailing slash in this string.  E.g.,
    #     NUPACKHOME = '/home/justin/nupack2.0'
    #
    # The call to pairs uses that -cutoff option, which is not
    # described in the NUPACK User Guide.  This enables calculation of
    # even small probabilities and therefore more accurate calculation
    # of n(s*), which is ultimately what we'll use the pair
    # probability matrix for.
    #
    # The pair probability matrix, pairprobs, is returned.
    # pairprobs[i][j] is the probability that base i and j are paired
    # at equlibrium.  Also returns unpairprobs, where unpairprobs[i]
    # is the probability that base is unpaired at equilibrium.  Note
    # that indexing begins at zero and only the upper triangle of the
    # pair probability matrix is nonzero.

    import os
    import re
    whiteSpaceSearch = re.compile('\s+')
    import random

    # Sequence length
    seqlen = len(sequence)

    # Initialize pairprobs and unpairprobs
    pairprobs = []
    for i in range(seqlen):
        pairprobs.append([0.0]*seqlen)
    unpairprobs = [0]*seqlen

    # Name of file for I/O
    fname = 'junk_file_delete_me%d' % random.randint(100000,999999)
    outputfile = '%s.ppairs' % fname

    # Make the input file
    inputFileName = '%s.in' % fname
    f = open(inputFileName,'w')
    f.write('%s\n' % sequence)
    f.close()

    # Create the command to compute the energy
    cmd = '%s/bin/pairs -cutoff 0.0 -T %.1f -material %s %s' % (NUPACKHOME,T,material,fname)

    # Run the command in a subshell
    os.system(cmd)

    # Parse the output
    f = open(outputfile,'r')

    # Blow through comments and blank lines
    line = f.readline()
    while not line[0].isdigit():
        line = f.readline()

    # Now we are at the entry containing number of bases
    line = f.readline()

    # Now we're at the pair probabilities, read until the end of the file
    while line != '':
        lineData = whiteSpaceSearch.split(line)
        i = int(lineData[0]) - 1
        j = int(lineData[1]) - 1
        if j == seqlen: # unpaired probability
            unpairprobs[i] = float(lineData[2])
        else:
            pairprobs[i][j] = float(lineData[2])
        line = f.readline()
    
    f.close()
    
    # Remove junk files
    cmd = 'rm -f %s.in ; rm -f %s' % (fname,outputfile)
    os.system(cmd)

    return (pairprobs,unpairprobs)
# ############################################################### #



# ############################################################### #
def dotParen2PairList(structure):

    # Converts the input dot-paren structure into a list of pairs
    # called pairlist.  Also returns plist, where plist[i] is the base
    # to which base i is paired and -1 if base is unpaired.
    # For example, an input of '.(((...)))' returns
    # pairlist = [[2,10],[3,9],[4,8]] and
    # plist = [-1,10,9,8,-1,-1,-1,4,3,2]
    # NOTE THAT INDEXING BEGINS AT ZERO!
    #
    #  This only works for single-stranded structures, and not
    #  complexes with multiple strands.

    # Length of the sequence
    seqlen = len(structure)
    
    pairlist = []
    leftlist = []
    ind = 0
    # While loop steps through list.  Each left bracket is stored.
    # Whenever we get a right bracket, it is necessarily pair with
    # the last left bracket in leftlist.  This pair is documented
    # and the first entry in leftlist is then deleted.
    while ind < seqlen:
        if structure[ind] == '(':
            leftlist.append(ind)
        elif structure[ind] == ')':
            pairlist.append([leftlist[-1],ind])
            leftlist.remove(leftlist[-1])
        ind = ind + 1

    pairlist.sort()

    # Get plist
    plist = [-1]*seqlen
    for x in pairlist:
        plist[x[0]] = x[1]
        plist[x[1]] = x[0]

    return (pairlist,plist)
# ############################################################### #



def sample(count, strands, T, material, dangles):
   
  import subprocess, tempfile, os
  
  tmp = tempfile.NamedTemporaryFile(delete=False,suffix=".sample")
  tmp.close()
    # we close it here as some OS's have issues opening the same file
    # simultaneous. Will reopen it later to get the data back out.

    # Notes:
    # 1) As of NUPACK 3.0.2, 'sample' is in the standard distribution
    #    and does the Boltzmann sampling we need. So we look in $NUPACKHOME/bin first.  
    # 2) Otherwise, 'sample' should be found using your user $PATH.
    # 3) Beware that there is a standard BSD tool with that name: if you are using OS X, 
    #    make sure that your path to the nupack 'sample' occurs before /usr/bin or it may not find it correctly.

  if 'NUPACKHOME' in os.environ:
    sample_exec = os.environ['NUPACKHOME']+'/bin/sample'
  else:
    sample_exec = 'sample' 

  material_out = ['-material', material]

  dangles_out = ['-dangles', dangles.lower()]   # Note that self._dangles should always be the string as long as the calls to the setter always invert the int back into string form. NUPACK appears to use the lowercase string name as the dangles names.

  temp_out = ['-T', '{0}'.format( T )]


    # editing here by EW 1/26/2014 to update from NUPACK 2.1 command line options & output to NUPACK 3.0.2
  popen_params = [sample_exec, "-multi"] + material_out + dangles_out + temp_out + ["-samples", str(count)]
  p = subprocess.Popen( popen_params,stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # new: NUPACK 3.0.2 takes the file name as user input
  input_str = "{0}\n{1}\n{2}\n{3}\n".format( tmp.name[:-7],
                                             len(strands),
                                             "\n".join( strands ),
                                             " ".join( [str(i+1) for i in range( len( strands ))])
                                           )
  result = p.communicate(input_str)[0]
    # note we toss the result as it's mostly just spam from the subprocess

  f = open(tmp.name, "rt")
  lines = f.readlines()

  f.close()
  os.remove(tmp.name) # was created by us [NamedTemporaryFile] and
                        # used by the sampler, thus we need to clean it up.
  if not "NUPACK 3.0" in lines[0]:
    raise IOError("Boltzmann sample function is not up to date. NUPACK 3.0.2 or greater needed.")

  sampled = [l.strip() for l in lines[14:]]
  return sampled
