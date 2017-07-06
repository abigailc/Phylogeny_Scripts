#!/usr/bin/python



#this requires internet access, uses eutils to get each bit of data.

#last edit abigailc@ACTAEON March 20somthing.
#purpose: get accession numbers from a blast output (or gi number maybe?).
#make a new fasta file with the FULL SEQUENCE each is pulled from, instead of the ALIGNED SEQUENCE like blast+ outputs.

#RIGHT NOW i am optimizing this to work with the results of BlastOutFasta
#in the future it should be made compliant with other seqID formats.

#BLAST OUT FASTA headers are:


#>gi|949160082|gb|KRP31410.1|chromosome_segregation_protein_SMC_[Actinobacteria_bacterium_BACL2_MAG-120802-bin41]
#MTLKGFKSFAAPTTLKFEPGITCVVGPNGSGKSNVVDALSWVMGEQGAKSLRGGKMEDVIFAGTSGRAPLGRAEVSVTID
#SEQEUENCE

#1. get the accession number. split by "|", take the 4th.
#done gen_blastout_accnums
#2. build the queries. batch submit if possible, I am not sure what that would do though? lets not for now (oops)


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
        self.newsequences = []

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
    def gen_blastout_accnums(self, depth=3):
        depth = depth
        print("finding acc numbers")
        self.numbers = []
        for item in self.ids:
            il = item.split("|")
            num = il[depth]
            self.numbers.append(num)
    #gets number (gi num or acc. num from shortened format)
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
    #gets taxID

   
    

    def Remove_Duplicates(self):
        donenums = []
        to_remove = []
        for item in self.numbers:
            if item in donenums:
                to_remove.append(self.numbers.index(item))
            else:
                donenums.append(item)
        for item in sorted(to_remove, reverse = True):
            self.ids.pop(item)
            self.seqs.pop(item)
            self.numbers.pop(item)
        #returns ids, seqs, and numbers without any accession number duplicates.
            
            
    def GetSeqs(self):
        print("getting new full sequence information... this might take a while")
        print("total number of unique acc numbers: "+str(len(self.numbers)))
        total = len(self.numbers)
        current = 0
        self.newsequences = [] 
        error_urls = 0 
        for item in self.numbers:
            #printing progress
            done = float(current)/float(total)#will be, like, 0.2
            percent = done*100 #will be, like, 20
            if percent.is_integer() is True:                
                if int(percent)%10 == 0:
                    if int(percent) == 0:
                        pass
                    else:
                        print(str(int(percent))+" percent done")

            #doing it
            accnum = item
            
            #####new bit added nov 16: no longer requires XMLLINT, does require minidom, urllib2 (probably are standard)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+item+"&retmode=xml"
            # define XML location w. current accnum
            try:
                dom = minidom.parse(urllib2.urlopen(url)) # parse the data
            except:
                print("error with url:"+ url)
                print("trying once more, than will skip this sequence.... ")
                try:
                    dom = minidom.parse(urllib2.urlopen(url))
                except:
                    error_urls += 1
                    continue
                    #raise SystemExit
            staffs = dom.getElementsByTagName("GBSeq_sequence")
            #print (staffs)
            for staff in staffs:
                
                s=staff.firstChild.data
                s = s.upper()
                #print("should be a sequence: "+s)
               # sequence = staff.getElementsByTagName("GBSeq_sequence")[0]
               # print("seq"+sequence)
            
            self.newsequences.append(s)
            current += 1
            #self.taxonomy is a list of dictionaries, each of which has an entry for each rank.
        self.seqs = self.newsequences
        if error_urls != 0:
            print("there was an error in URL access for : "+str(error_urls)+" sequences")
        return self.newsequences
    def gen_new_fasta(self, new_fasta_name):
        #this should print the changed seqids and changed AA sequences to file.
        newfasta = new_fasta_name
        # print(len(self.original_ids))
        # print(len(self.ids))
        # print(len(self.original_seqs))
        # print(len(self.seqs))
        with open (newfasta, "w") as new:
            for i in range(len(self.ids)):
                new.write(">"+self.ids[i].strip()+"\n")
                # print(i)      #
                #unclear if this needs a "\n" after it... check.#TODO
                                #print(self.seqs)
                                #print(type(self.seqs[i]))
                new.write(self.seqs[i]+"\n")
        print("Finished, your new fasta file is located at "+newfasta)
        #done


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
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("-d", "--depth", action = "store", default = 3, help="give depth of accession number location if not 3. start counting at 0")



    
    args = parser.parse_args()
    #change dir if given
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")
    #run the thing
    MyFasta = Fasta(args.fasta)
    MyFasta.gen_original_lists()
    MyFasta.gen_blastout_accnums(int(args.depth))
    MyFasta.Remove_Duplicates()
    MyFasta.GetSeqs()
    #get a name for the new fasta
    try:
        base, ext = args.fasta.split(".")
        newfastaname = base+"_Full.fasta"
    except:
        newfastaname = args.fasta+"_Full.fasta"
    MyFasta.gen_new_fasta(newfastaname)
    print("done, your new file it at: "+newfastaname)
