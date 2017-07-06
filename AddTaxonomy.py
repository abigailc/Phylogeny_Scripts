#!/usr/bin/python

#refresh after NCBI update on 11/9

#this requires internet access, uses eutils to get each bit of data.

#could also be built locally from the taxonomy database + crawling upwards parsing.
#implement that if it would save time for the larger dataset.

#testing 11/17


#last edit abigailc@ACTAEON
#purpose: adds taxonomy information to .fasta file sequence IDs
#purpose2: smaller and more customizable than implementation within FISH_2.py
from __future__ import print_function
#^for progress bar
###CLASS
class Fasta:
    def __init__(self, name):
        #this is the to-be-modified version of sequence IDs and sequence-Data
        # ALWAYS keep IDS and SEQS the same length. id[1] should ALWAYS correspond to seq[1].
        self.name = name
        self.ids = []
        self.seqs = []
        # these are the original SEQids and Sequences. They should never be modified after generation in gen_original_lists or blast_to_fasta
        self.original_ids = []
        self.original_seqs = []
        #obsolete
        self.species_names = []
        #list of gi / accession numbers
        self.numbers = []
        #list of taxid numbers
        self.taxid = []
        #list of dictionaries of taxonomy rank information from NCBI
        self.taxonomy = []
    #loads your file as IDS and SEQUENCES
    def gen_original_lists(self):
        fastaname = self.name
        with open(fastaname) as fastafile:
            for line in fastafile:
                if "\n" == line:
                    pass
                if ">" in line:
                    #write the previous AA seq
                    try:
                        AAseq=AAseq.strip()
                        self.seqs.append(AAseq)
                        self.original_seqs.append(AAseq)
                    except:
                        pass
                        #initialize a new AAseq
                    AAseq = ""
                    #format the seqID
                    newline = line.strip()
                    newline = line.strip(">")
                    #write the seqID
                    self.ids.append(newline.strip())
                    self.original_ids.append(newline.strip())
                else:
                    AAseq = AAseq+line
                    AAseq=AAseq.strip()
            #catch the last AAseq pass
            self.seqs.append(AAseq)
            self.original_seqs.append(AAseq)
        print("Initial sequence and ID lists created. Contains "+str(len(self.ids))+" sequences")
    #gets number (gi num or acc. num)
    def gen_numbers(self):
        print("Finding identification numbers...")
        for item in self.ids:
            number = re.sub("(.*)(\|)(.*)","\\3", item)
            try:
                n = int(number)
            except:
                do = "N"
                nlist = item.split("|")
                for t in nlist:
                    if do == "Y":
                        number = t
                        break
                    if t == "gi#":
                        do = "Y"
            self.numbers.append(number)

    def NumbersList(self, numfile):
        self.numbers = []
        with open(numfile) as num:
            for line in num:
                self.numbers.append(line.strip())
    #gets taxID
    def SetTaxID(self):
        from time import sleep
        print("acquiring taxids")
        self.taxid = []
        current = 0
        total = len(self.numbers)
        #print(self.numbers)
        for item in self.numbers:
         #   print(item)
            #because we all need a goddamn progress bar.
            current += 1
            printProgressBar(current, total, prefix = 'Progress:', suffix = 'Complete', length = 100)
           

            #####new bit added nov 16: no longer requires XMLLINT, does require minidom, urllib2 (probably are standard)
            #print(item)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+item+"&retmode=xml" # define XML location w. current taxid
            try:
                dom = minidom.parse(urllib2.urlopen(url)) # parse the data
            except:
                try:
                    time.sleep(10)
                    dom = minidom.parse(urllib2.urlopen(url)) 
                except:
                    print("couldnt parse url: "+url)
                    continue
            staffs = dom.getElementsByTagName("GBQualifier_value")
            #print (staffs)
            taxonid = "NA"
            for staff in staffs:
                s=staff.firstChild.data
                if "taxon:" in s:
                    #print (s)
                    try:
                        tax, num = s.split(":")
                    except:
                        print("error in parsing taxonid")
                        print(s)
                        blah, good = s.split("taxon:")
                        num = re.sub("([0-9]*)(.*)", "\\1", good)
                        print("tried to fix:"+num)
                    taxonid = num
            taxonid.strip()
            if taxonid == "NA":
                print("Error with "+item)
            self.taxid.append(taxonid)
            
    def GetTaxonomy(self):
        print("getting taxonomic information... this might take a while")
        total = len(self.taxid)
        current = 0
        self.taxonomy = []
        if self.taxid == []:
            print("You need to generate taxids first.. lets try")
            self.SetTaxID()
            total = len(self.taxid)
            
        for item in self.taxid:
            current += 1
            printProgressBar(current, total, prefix = 'Progress:', suffix = 'Complete', length = 100)
            taxid = item
            ranklist = "superkingdom kingdom phylum class order family genus"
            ranklist = ranklist.split()
            taxdict = {}
            #####new bit added nov 16: no longer requires XMLLINT, does require minidom, urllib2 (probably are standard)
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id='+item # define XML location w. current taxid
            try:
                dom = minidom.parse(urllib2.urlopen(url)) # parse the data
            except:
                print("error with url:"+ url)
                print(taxid+" is getting all NA ranks")
                for r in ranklist:
                    taxdict[r] = "NA"
            staffs = dom.getElementsByTagName("Taxon")
            #print (staffs)
            for staff in staffs:
                    tid = staff.getElementsByTagName("TaxId")[0]
                    tname = staff.getElementsByTagName("ScientificName")[0]
                    trank = staff.getElementsByTagName("Rank")[0]
                    taxid = tid.firstChild.data
                    taxname = tname.firstChild.data
                    taxrank = trank.firstChild.data
                    #print("taxid:%s, taxname:%s, taxrank:%s" %(taxid, taxname,taxrank))
                    if taxrank in ranklist:
                        #cleanup taxname here. going to enable.
                        taxname = taxname.replace(".", "")
                        taxname = taxname.replace(" ", "")
                        taxdict[taxrank] = taxname
            for r in ranklist:
                if r in taxdict:
                    pass
                else:
                    taxdict[r] = "NA"
            
            self.taxonomy.append(taxdict)
            #self.taxonomy is a list of dictionaries, each of which has an entry for each rank. eg 
            #self.taxonomy = [{cat}, {dog}, {wolf}]
            #self.taxonomy[0] = {cat} = {genus:felis, species:catus, phylum:chordata }
        return self.taxonomy
    #this should modify the seqIDS from Species_name|ACC##### to Euk|Meta|Chord|Species_name|ACC#####
    def AppendTaxonomy(self, ranklist = "NA"):
        for item in self.ids:
            index = self.ids.index(item)
            rankdict = self.taxonomy[index]
            if ranklist == "NA":
                newitem = rankdict["superkingdom"]+"|"+rankdict["kingdom"]+"|"+rankdict["phylum"]+"|"+rankdict["class"]+"|"+rankdict["order"]+"|"+rankdict["family"]+"|"+rankdict["genus"]+"|"+item
                self.ids[index] = newitem
            else:
                a = ""
                for rank in ranklist:
                    a = a + rankdict[rank]+"|"
                a = a + item
                self.ids[index] = a


    def gen_new_list(self, newname, ranklist = "NA"):
        with open(newname, "w") as new:
            for i in range(len(self.numbers)):
                rankdict = self.taxonomy[i]
                a = ""
                for rank in ranklist:
                    a = a + rankdict[rank]+"|"
                new.write(numbers[i]+"\t"+a)

    def gen_new_fasta(self, new_fasta_name):
        #this should print the changed seqids and changed AA sequences to file.
        newfasta = new_fasta_name
        # print(len(self.original_ids))
        # print(len(self.ids))
        # print(len(self.original_seqs))
        # print(len(self.seqs))
        with open (newfasta, "w") as new:
            for i in range(len(self.original_ids)):
                new.write(">"+self.ids[i].strip()+"\n")
                # print(i)      #
                #unclear if this needs a "\n" after it... check.#TODO
                                #print(self.seqs)
                                #print(type(self.seqs[i]))
                new.write(self.seqs[i]+"\n")
        print("Finished, your new fasta file is located at "+newfasta)
        #done
    def gen_species_lists(self):
        for item in self.ids:
            taxon = re.sub("([^_]*)([A-Z][a-z]*_[a-z]*)(.*)", "\\2", item)
            if "#" in taxon:
                print ("TAXON error in gen_species_lists():" + taxon)
            self.species_names.append(taxon)
    def PrintTaxonInfo(self):
        self.gen_species_lists()
        try:
            a=self.name.split(".")
            taxonfile = a[0]+"_taxinfo.txt"
        except:
            taxonfile = self.name+"_taxinfo.txt"
        with open(taxonfile, "w") as new:
            for item in self.ids:
                index = self.ids.index(item)
                seqid = item
                taxid = self.taxid[index]
                sciname = self.species_names[index]
                ranks = "superkingdom kingdom phylum class order family genus"
                ranklist = ranks.split()
                taxonomy = ""
                for r in ranklist:
                    try:
                        tax = self.taxonomy[index][r]
                    except:
                        pass
                    taxonomy = taxonomy+";"+tax
                new.write(item+"\t"+sciname+"\t"+taxid+"\t"+taxonomy)
                    
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = 'X'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    length = 20
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    #print("\r")
    print(' %s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="\r")
    # Print New Line on Complete
    if iteration == total: 
        print()

# 
# Sample Usage
# 

#parser

if __name__ == "__main__":

    print("Running in terminal")
    import subprocess
    import sys
    import argparse
    import os
    import re
    from xml.dom import minidom
    import urllib2
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("-r", "--ranks", action = "store", default = False, help="give specific ranks to append, if you dont want them all. eg -r phylum")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("-i", "--info", action = "store_true", default = False, help="toggle saves taxonomy information to new file")
    parser.add_argument("-n", "--numbers_to_names", action = "store_true", default = False, help="toggle prints a file of numbers___names")


    
    args = parser.parse_args()
    #change dir if given
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")
    #given a list of acc numbers
    if args.numbers_to_names is True:
        MyFasta = Fasta(args.fasta)
        MyFasta.NumbersList(args.fasta)
        MyFasta.SetTaxID()
        MyFasta.GetTaxonomy()
        try:
            base, ext = args.fasta.split(".")
            newfastaname = base+"_Taxo.fasta"
        except:
            newfastaname = args.fasta+"_Taxo.fasta"

        #make ranks work
        if args.ranks is False:
            print("error specify ranks eg species or genus")
        else:
            if " " in args.ranks:
                ranklist = args.ranks.split()
            else:
                ranklist = [args.ranks]        
        MyFasta.gen_new_list(newfastaname, ranklist)
        raise SystemExit 
    #run the thing
    MyFasta = Fasta(args.fasta)
    MyFasta.gen_original_lists()
    MyFasta.gen_numbers()
    MyFasta.SetTaxID()
    MyFasta.GetTaxonomy()
    if args.ranks is False:
        MyFasta.AppendTaxonomy()    
    else:
        if " " in args.ranks:
            ranklist = args.ranks.split()
        else:
            ranklist = [args.ranks]        
        MyFasta.AppendTaxonomy(ranklist)
    #get a name for the new fasta
    try:
        base, ext = args.fasta.split(".")
        newfastaname = base+"_Taxo.fasta"
    except:
        newfastaname = args.fasta+"_Taxo.fasta"
    MyFasta.gen_new_fasta(newfastaname)
    if args.info is True:
        MyFasta.PrintTaxonInfo()
    print("done, your new file it at: "+newfastaname)
