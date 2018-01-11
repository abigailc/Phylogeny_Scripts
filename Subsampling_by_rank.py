#!/usr/bin/python

# last update 1/11/18 from abigailc@Leviathan
# pulled base code out of FEAST (4/12/16 LED)


#USAGE:
#python Sumsampling_by_rank /path/to/dir FASTA Keep_Number Keep_Depth (-Show_Errors) (-NA_number [NUMBER]) (-Special_String)
#python Sumsampling_by_rank /myproject/ Example.fasta 2 3 -Show_Errors -NA_number 1 -Special_String "cyanobacteria"
#^ would open the file Example.fasta and keep 2 sequences per class, all things that have "cyanobacteria" in the seqID, and 1 per order of any sequences with class = "NA"



#Main Function

#keepn is how many to keep
#keepdepth is which rank to consider eg Kingdom|Phylum|Class phylum in depth = 2, class is depth = 3
#NAn should be True if you indicated you want subsampling at a lower level if given rank is "NA"
#showerrors should be True if you want to print a list of error sequences.
#special should be a string of things you want to keep (ignoring all other ss rules). eg if you put "Cyanobacteria" it will do subsampling properly for non-cyanos, but just keep every single cyano.
def Subsample_New(fasta, keepn, keepdepth, NAn = 0, showerrors = False, special = "$%$"):
    #make an output file name
    fa, ex = fasta.split(".")
    output = fa + "SubSamp" + keepn + "p" + keepdepth + ".fasta"
    #read in all data at given rank: keepdepth
    MLO = MakeLists(fasta, keepdepth)
    #prints error sequences if you ask it to.
    if showerrors == True:
        if NAnum == 0:
            ShowErrorsAll(fasta, output, MLO[2], MLO[1])
        else:
            ShowErrors(fasta, output, MLO[2])
    KE = SeqsToKeep(MLO, keepn, NAn)
    return GetSeqs(fasta, KE, output, special)


#this function is a mess and you should rewrite it because its ugly and half-unused.
def MakeLists(fasta, rank):
    #make sure depth is an int
    number = int(rank)
    #initiate dictionaries
    diTRL = {}
    diNA = {}
    #initiate error tracker
    listofERROR = []
    with open(fasta) as old:
        for line in old:
            #only look at seqids
            if ">" in line:
                if "gi#" in line:
                    gitax = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\5~\\2", line)

                    gi, tax = gitax.split("~")
                    tn = tax[:-1]
                    taxlist = tn.split("|")
        ##                print (taxlist)
                    if taxlist[0] == "NT":
                        listofERROR.append(gi)
                        continue
                    thetax = taxlist[number - 1]
                else:
                    taxlist = line.split("|")
        ##                print (taxlist)
                    if taxlist[0] == "NT":
                        listofERROR.append(gi)
                        continue
                    thetax = taxlist[number - 1]
                    gi = line
                if thetax == "NA":
                    try:
                        thetax = taxlist[number]
                    except:
                        pass
                    if thetax == "NA":
                        try:
                            thetax = taxlist[number + 1]
                        except:
                            pass
                        if thetax == "NA":
                            listofERROR.append(gi)
                            continue
                        if thetax in diNA:
                            diNA[thetax].append(gi)
                        else:
                            diNA[thetax] = [gi]
                elif thetax in diTRL:
                    diTRL[thetax].append(gi)
                else:
                    diTRL[thetax] = [gi]
    print("Created lists")
    return diTRL, diNA, listofERROR


def ShowErrorsAll(fasta, output, erlist, nadic):
    n, e = output.split(".")
    num = 0
    errorfa = n + "_SubErr.txt"
    with open(fasta) as old:
        with open(errorfa, "w") as new:
            for line in old:
                if ">" in line:
                    gi = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)(\n)", "\\5", line)
##                    print (gi)
                    if gi in erlist:
                        num += 1
                        new.write(line)
                    for thing in nadic:
                        if gi in nadic[thing]:
                            num += 1
                            new.write(line)
    print("There were " + str(num) + " sequences that were ignored for lack of taxonomic information. \n Their SeqIDs were saved in " +
          errorfa + " if you'd care to see.")

def ShowErrors(fasta, output, erlist):
    import re
    n, e = output.split(".")
    num = 0
    errorfa = n + "_SubErr.txt"
    with open(fasta) as old:
        with open(errorfa, "w") as new:
            for line in old:
                if ">" in line:
                    gi = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)(\n)", "\\5", line)
##                    print (gi)
                    if gi in erlist:
                        num += 1
                        new.write(line)
    print("There were " + str(num) + " sequences that couldn't be sorted due to lack of taxonomic information. \n Their SeqIDs were saved in " +
          errorfa + " if you'd care to see.")


def SeqsToKeep(makelistsout, number, NAnum):
    diTRL = makelistsout[0]
    number = int(number)
    diNA = makelistsout[1]
    import random
    keep = []
    for entry in diTRL:
        try:
            r = random.sample(diTRL[entry], number)
        except ValueError:
            le = len(diTRL[entry])
            r = random.sample(diTRL[entry], le)
        for thing in r:
            keep.append(thing)
##    print (keep)
    for entry in diNA:
        try:
            r = random.sample(diNA[entry], NAnum)
        except ValueError:
            le = len(diNA[entry])

            r = random.sample(diNA[entry], le)
        for thing in r:
            keep.append(thing)
    print("Picked a random subset of each group")
    return keep


def GetSeqs(fasta, keep, output, special):
    copy = "no"
    with open(fasta) as old:
        with open(output, "w") as new:
            for line in old:
                if ">" in line:
                    gi = re.sub(
                        "(>)(.*)(\|)(gi#?\|?)([0-9]*)(.*)(\n)", "\\5", line)
                    if gi in keep:
                        copy = "yes"
                        new.write(line)
                        # to stop it from copying multiple sequences in the
                        # case of shared gi numbers All(whole genome seq etc)
                        keep.remove(gi)
                    elif special in line:
                        copy = "yes"
                        new.write(line)
                    else:
                        copy = "no"
                else:
                    if copy == "yes":
                        new.write(line)
            new.flush()
    print("Created a subset file at " + output)
    return output


##PARSER
if __name__ == "__main__":

    print("Running subsampler")

    import sys
    import argparse
    import os
    import re

    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(
    ), type=str, help="type name of directory to run in (where .fasta(s) reside)")
    parser.add_argument("FASTA", default="no_fasta_specified", type=str,
                        help="type the name of your .fasta file or a space seperated list of files within quotes")

    parser.add_argument("Keep_Number", default="1", type=str,
                        help="type how many sequences you want per rank")

    parser.add_argument("Keep_Depth", default="1", type=str,
                        help="type the depth (by bar seperations >1|2|3|4) in your seqIDs of desired rank. eg '2'")
    parser.add_argument("-Show_Errors", "--se", action = "store_true", default=False, help="Prints a list of prolematic sequences / sequences with missing data")
    parser.add_argument("-NA_number", "--na", action = "store", default = 0, help="toggles 'if rank=NA, keep n sequences of one-rank-lower.' provide the n. default = 0.")
    parser.add_argument("-Special_String", "--ss", action = "store", default="$$$$$$$$", help="toggles 'ignore ss and keep ALL sequences contining the following string:'. provide string.")
    
#USAGE:
#python Sumsampling_by_rank /path/to/dir FASTA Keep_Number Keep_Depth -Show_Errors -NA_number, -Special_String
#python Sumsampling_by_rank /myproject/ Example.fasta 2 3 -Show_Errors -NA_number 1 -Special_String "cyanobacteria"
#^ would open the file Example.fasta and keep 2 sequences per class, all things that have "cyanobacteria" in the seqID, and 1 per order of any sequences with class = "NA"


#load args
args = parser.parse_args()

#change directory
directory = args.directory
os.chdir(directory)

#load fasta / list of fastas
fastalist = args.FASTA
try:
    a = fastalist.split()
except:
    a = [fastalist]

#one at a time run subsampling
for fasta in a:
    #avoid a dumb thing.
    if "./" in a:
        a = a[:-2]
    print("Beginning SubSampling on "+fasta)
    new_fasta = Subsample_New(fasta, args.Keep_Number, args.Keep_Depth, args.na, args.se, args.ss)
    print("Done")
    

print("Finished! Exiting.")







