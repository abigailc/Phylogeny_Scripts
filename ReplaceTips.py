#!/usr/bin/python



#@abigailc

#created 5/16/16 on actaeon

#update 11/11/17 on leviathan



#this correlates first ten protein id characters with bracketed species names

#for gregs thing

def MakeDictBlastOut(fasta):

    print("running with -s")

    CTdict = {}

    snlist = []

    iteration = 0

    with open(fasta) as old:

        for line in old:

            if ">" in line:

                iteration += 1

                line = line[:-1]

                pidsn = re.sub("(>gi\|)([0-9]*\|[a-z]*\|)(.{10})([^\[]*\[)([^\]]*)(].*)", "\\3~\\5", line)

                try:

                    pid, sn = pidsn.split("~")

                except:

##                    print("ERROR")

                    pid = re.sub("(>gi\|)([0-9]*\|[a-z]*\|)(.{10})(.*)", "\\3", line)

                    sn = re.sub("([^\[]*)(\[)([^\]]*)(.*)", "\\3", line)

                sn = re.sub("([\.;:\-\)\/\( ])", "_", sn)

                sn = re.sub("(___*)", "_", sn)

                print (sn)

                if sn in snlist:

                    sn=sn+"_copy"+str(iteration)

                snlist.append(sn)

                CTdict[pid] = sn



    return CTdict



#this correlates entire line to any combination of : depth-taxonomy(s) and gi | given a "Shorten" style format.

def MakeDictTaxo(fasta, replace):

    CTdict = {}

    with open(fasta) as old:

        for line in old:

            if ">" in line:

                line = line[:-1]

                if "|gi#" in line:

                    taxgi = re.sub("(>)([^#]*)(\|gi#\|?)([0-9]*)(.*)", "\\2~\\4", line)

                    tax, gi = taxgi.split("~")

                    taxlist = tax.split("|")

                    if replace == "gi":

                        CTdict[line[1:]] = gi

                    if type(replace) is int:

                        CTdict[line[1:]] = taxlist[replace-1]

                    if type(replace) is str:

                        listreplace = replace.split()

                        newid = ""

                        for item in listreplace:



                            if item == "gi":

                                newid = newid+"|"+gi

                            else:

                                newid = str(newid)+"|"+str(taxlist[int(item)-1])

                        newid = newid[1:]

                        CTdict[line[1:]] = newid

                        print(newid)

                else:

                    tax = re.sub("(>)([^#]*)(\|gi#\|?)([0-9]*)(.*)", "\\2", line)

                    taxlist = tax.split("|")

                    if replace == "gi":

                        pass

                    if type(replace) is int:

                        CTdict[line[1:]] = taxlist[replace-1]

                    if type(replace) is str:

                        listreplace = replace.split()

                        newid = ""

                        f = 1

                        for item in listreplace:

                            f += 1

                            if item == "gi":

                                newid = newid+"|NA"

                            else:

                                #SPECIFICALLY FOR CURRENT USE_CASE, REMOVE LATER

                                if f == 2:

                                    newid = str(newid)+"|"+str(taxlist[int(item)-1])

                                if f == 3:

                                    newid = str(newid)+"|"+str(taxlist[int(item)])

                        newid = newid[1:]

                        if ">" in newid:

                            newid = newid[1:]

                        CTdict[line[1:]] = newid

                        print(newid)

    return CTdict



def MakeDictTenChar(fasta):

    CTdict = {}

    iteration = 0

    try:

        name, extention = fasta.split(".")

        newfasta = name+"_TenChar.fasta"

    except:

        newfasta = fasta+"_TenChar.fasta"

    with open(fasta) as old:

        with open(newfasta, "w") as new:

            for line in old:

                if ">" in line:

                    iteration +=1

                    line = line.strip()

                    if iteration > 999999999:

                        print("ITERATION ERROR, shutting down. don't run with > 999999999 sequences?")

                    filler = 9-iteration

                    fill = ""

                    for i in range(filler):

                        fill = fill+"0"

                    newid = fill+str(iteration)+"x"

                    new.write(">"+newid+"\n")

                    CTdict[line] = newid

                else:

                    new.write(line)

    return CTdict, newfasta





def MakeDictTenCharSpecNam(fasta):

    CTdict = {}

    iteration = 0

    try:

        name, extention = fasta.split(".")

        newfasta = name+"_TenChar.fasta"

    except:

        newfasta = fasta+"_TenChar.fasta"

    with open(fasta) as old:

        with open(newfasta, "w") as new:

            for line in old:

                if ">" in line:

                    iteration +=1

                    line = line.strip()
                    #print(line)
                    if iteration > 99999:

                        print("ITERATION ERROR, shutting down. don't run with > 99999 sequences?")



##                    #i have something like

##                    >Methanococcoides_burtonii|gi|909890

##                    #i want

##                    MethBurt00

                    GenusSpecies = re.sub("(>)([A-Z][a-z]*)(_)([A-Z]*[a-z]*)(.*)", "\\2~\\4", line)
                   # print(GenusSpecies)
                    try:
                        Genus, Species = GenusSpecies.split("~")
                        g4 = Genus[:4]
                        try:
                            s4 = Species[:4]
                        except:
                            s4 = Species[:2]
                        if iteration < 10:
                            newid = g4+s4.capitalize()+"0"+str(iteration)
                        elif iteration < 100:
                            newid = g4+s4.capitalize()+str(iteration)
                        elif iteration < 1000:
                            newid = g4+s4.capitalize()[:-1]+str(iteration)
                        elif iteration < 10000:
                            newid = g4[:-1]+s4.capitalize()[:-1]+str(iteration)
                        else:
                            newid = g4+s4.capitalize()+str(iteration)
                     #   print("1"+newid)
                    except:

##                        print(GenusSpecies)

                        gs8 = GenusSpecies[1:9]

                        if iteration < 10:
                            newid = gs8+"0"+str(iteration)
                        elif iteration < 100:
                            newid = gs8+str(iteration)
                        elif iteration < 1000:
                            newid = gs8[:-1]+str(iteration)
                        elif iteration < 100000:
                            newid = gs8[:-6]+gs8[4:6]+str(iteration)
                        else:

                            print("too many sequences")

                            raise SystemExit

##                        print(newid)

                    new.write(">"+newid+"\n")
                    
                    #print(newid)

                    CTdict[line] = newid

                else:

                    new.write(line)

    return CTdict, newfasta



#this parses a one-line tree (newick format; will also work on newick trees embedded in nexus, and will also change thing is nexus taxa blocks (usually!)

#and converts things in CTdict(keys) to CTdict(values).

#so, full id to just to species or what have you.

def ConvertTreeTips(newick, CTdict, fasta):

    infoid = PrintDict(fasta, CTdict)

    print("In case of error, or desired changes, dictionary is being printed to "+infoid+" for future reference. Call with -fi")

    nnewick = ""

    if "." in newick:

        for item in newick:

            if item == ".":

                pass

            else:

                nnewick= nnewick+item

    newnewick = nnewick+"replaced.tree"

    with open (newick) as old:

        with open (newnewick, "w") as new:

            for line in old:

                for item in CTdict:

##                    print(item)

                    line = line.replace(item, CTdict[item])

##                print(line)

                new.write(line)

    print("finished, gi-newick file at: "+newnewick)



def PrintDict(fasta, ctdict):

    try:

        name, ext = fasta.split(".")

        info = name+"_info.txt"

    except:

        info = fasta+"_info.txt"

    with open(info, "w") as inf:

        for item in ctdict:

            inf.write(item.strip()+"\n"+ctdict[item.strip()].strip()+"\n\n")

    return info





def ParseInfo(infofile):

    kid = "no"

    vid = "no"

    CTdict = {}

    with open (infofile) as old:

        for line in old:

            if kid == "no":

                key = line.strip()

##                print("k:"+key)

                kid = "yes"

                continue

            elif kid == "yes":



                if vid == "yes":



                    if len(value) > 11:

                        print("Error: value too long! key:"+key)

                        print("value:"+value)

##                    print("dict")

                        ####want it backwards this time

                    if ">" in key:

                        key = key[1:]

                    if ">" in value:

                        value = value[1:]

                    CTdict[value]=key

                    vid = "no"

                    kid = "no"

                    continue

                if vid == "no":

                    value = line.strip()

##                    print("V:"+value)

                    vid = "yes"

        CTdict[key]=value

##    print(CTdict)

    return CTdict



def TenCharSwap(fasta, tree = "NA"):

    CTdict, newfasta = MakeDictTenChar(fasta)

    print("Please use the fasta generated at "+newfasta+" for further analyses. Once you have created your tree, enter the resulting filename for tip re-conversion.")

    tree = input("..... type new Tree File Name:")

    ConvertTreeTips(tree, CTdict, fasta)



def TenCharSwapSpecNam(fasta, tree = "NA"):

    CTdict, newfasta = MakeDictTenCharSpecNam(fasta)

    print("Please use the fasta generated at "+newfasta+" for further analyses. Once you have created your tree, enter the resulting filename for tip re-conversion.")

    tree = input("..... type new Tree File Name:")

    ConvertTreeTips(tree, CTdict, fasta)





def ConvertTipsFromInfo(info, tree, fasta):

    CTdict = ParseInfo(info)

    ConvertTreeTips(tree, CTdict, fasta)



if __name__ == "__main__":



    print("Running in terminal")

    import sys

    import argparse

    import os

    import re

    parser = argparse.ArgumentParser(description="All")

    parser.add_argument("directory", type=str, help="type name of directory to run in (where .nex resides)")

    parser.add_argument("FASTA", type=str, help="type the name of your fasta file w original seqids")

    parser.add_argument("-t", "--tree", dest = "TREE", action = "store", default = "NA", help="give filepath to tree for tip-replacing")

    parser.add_argument("-re", "--replace", dest = "RE", action = "store", default = "NA", help="type gi or bar-indent number to replace full line with, space-sep in quotes eg /'1 4 gi/'")

    parser.add_argument("-sp", "--species", action = "store_true", default = "NA", help="replaces first ten digits of pid with species ID if fasta in blastout brackets")

    parser.add_argument("-co", "--convert", action = "store_true", default = "NA", help="convert tipIDs to unique ten characters in fasta, with option to convert back after tree creation")

    parser.add_argument("-fi", "--fileinfo", action = "store", default = "NA", help="convert tipIDs in tree to original form as stored in infofile specified")

    parser.add_argument("-te", "--tencharspecnam", action = "store_true", default = "NA", help="forspeciesnameconversion")





    args = parser.parse_args()

#making input list

    import re

    try:

        os.chdir(args.directory)

        print("changed dir to "+args.directory)

    except:

        print ("didn't change dir")



    if args.RE == "NA":

        if args.species == False:

            if args.convert == False:

                if args.fileinfo == "NA":

                    print ("specify a job next time")

                    raise SystemExit

    if args.species == True:

        print("...")

        a = MakeDictBlastOut(args.FASTA)

    elif args.tencharspecnam == True:

        done = TenCharSwapSpecNam(args.FASTA, args.TREE)

    elif args.convert == True:

        done = TenCharSwap(args.FASTA, args.TREE)

        raise SystemExit

    elif args.fileinfo != "NA":

        ConvertTipsFromInfo(args.fileinfo, args.TREE, args.FASTA)

        raise SystemExit

    else:

        a = MakeDictTaxo(args.FASTA, args.RE)

    ConvertTreeTips(args.TREE, a, args.FASTA)

    print("Exiting")

