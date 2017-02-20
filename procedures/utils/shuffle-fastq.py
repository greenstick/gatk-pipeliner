#! /home/users/cordier/.linuxbrew/bin/python3

if __name__ == "__main__":

    # Imports
    import argparse
    import random
    # Library Imports
    from Bio import SeqIO

    # 
    # Parse Arguments
    # 

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Input Paired-End FastQ File")
    parser.add_argument("-o", "--outputprefix", type = str, help = "Prefix for Output FastQ File (i.e. prefix.shuffled.fastq)")
    parser.add_argument("-s", "--chunksize", type = int, help = "Size of Chunk (i.e. n Pairs of Reads) to Shuffle at a Time. Default = 50000")
    argsDict = vars(parser.parse_args())

    fastq = argsDict["input"]
    prefix = argsDict["outputprefix"]
    chunksize = argsDict["chunksize"]

    # Assertions for Required Input
    assert (fastq is not None), "No FastQ input provided!"

    #  If No Prefix, Use Same as FASTQ
    if prefix is None: 
        prefix = ".".join(fastq.split(".")[0:-1])

    # If No Chunk Size, Set Default
    if chunksize is None:
        chunksize = 50000
    assert (chunksize % 2 == 0), "Chunksize must be even!"
    
    # 
    # Conversion
    #

    # Open FastQ, Chunk File, Shuffle, & Write Chunk to Output
    with open(fastq, "r") as handle, open(prefix + ".shuffled.fastq", "w") as output:
        records = SeqIO.parse(handle, "fastq")
        chunk = []
        i, nrecords = 0, 0
        for record in records:
            try:
                chunk.append([record, next(records)])
                if len(chunk) > chunksize:
                    random.shuffle(chunk) # Inplace
                    flattened = [item for sublist in chunk for item in sublist]
                    SeqIO.write(flattened, output, "fastq")
                    chunk = []
                    nrecords += len(flattened)
                    i += 1
                    print("Shuffle FastQ: %d Chunk(s) Written" % (i))
            except:
                random.shuffle(chunk) # Inplace
                flattened = [item for sublist in chunk for item in sublist]
                SeqIO.write(flattened, output, "fastq")
                chunk = []
                nrecords += len(flattened)
                i += 1
                print("Shuffle FastQ: %d Chunk(s) Written" % (i))
        # Write Tail
        try:
            random.shuffle(chunk) # Inplace
            flattened = [item for sublist in chunk for item in sublist]
            SeqIO.write(flattened, output, "fastq")
            nrecords += len(flattened)
            i += 1
            print("Shuffle FastQ: %d Chunk(s) Written (%d Records)" % (i, nrecords))
        except:
            pass

else:

    pass
