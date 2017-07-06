# #!/usr/bin/python

# last edit abigailc@Actaeon on october 21 2016


#major speed increase... first run should take 2-4 hours, subsequent runs from same directory ~1 hour.
#usage
#assume you want to use species that are found in the single_gene_fasta file example.fasta, and genes as found in Euk_Ribo.fasta
#always use -r flag for now, it is faster.
#if you want to use species from a list file (one species per line) use -s SPECIESFILE instead of -f FASTAFILE
# $python MakesSpeciesTrees.py -f example.fasta -g Euk_Ribo.fasta -p MyTree1 -r

from __future__ import print_function

##SET
names_file = "/Users/abigailc/Documents/Taxonomy_Stuff/taxdump/names.dmp"
nodes_file = "/Users/abigailc/Documents/Taxonomy_Stuff/taxdump/nodes.dmp"



######### PERSONAL_SETTINGS #########
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe4.mit.edu"
clus_head = "abigailc@eofe4.mit.edu:/home/abigailc/"
Path_Blast = "/Users/abigailc/blast/"
New_Gene_List = "/Users/abigailc/Documents/TestForDani/Newest_Ribo.fasta"


#if running elsewhere there MIGHT be issues.




#from Parse_Taxonomy import *
from Classes_Standalone import *

a = Subtree("lol")

import sys
import argparse
import os
import re
import time
import subprocess
    
from xml.dom import minidom
try:
    import urllib2
except:
    pass

# #to use this script, you need to do several things.
# 1. set the ssh path so you will be able to run things on the cluster.
# 2. set the cluster_head path so files will be transferred to/from your partition (not mine!)
# 3. run from the command line
# a. cd to directory you are keeping this script.
# b. identify the genes you want your species tree to be built from, and create a .fasta file containing them. For ease of understanding some produced files, I am using the sequence id convention
# >GeneName|other_sequence_information_blah_blah
# You can use the files I have here (Archaeal_Ribo.fasta, Bacterial_Ribo.fasta, Eukaryal_Ribo.fasta) or create your own query file. Please provide the full path to the query file, or keep it in the same folder as the directory you indicate when calling MakeSpeciesTree.
# indicate the gene-queries to use with the -genes tag.
# c.  identify the species that should be included in your species tree. They should either be in a shorten-formatted .fasta file, or in a plaintext file with one species name on each line of the file eg.
# Cat
# Dog
# Rat
# indicate the species file with -species or -fasta if you are giving a shortened fasta instead.
# d. if you want your output to be named something reasonable (like SOD_cyanos_blah) give it a prefix using -p
# e. run something like $python MakesSpeciesTrees.py -s species_file.txt -g bacterial_ribo_gene_seqs.fasta -p MyProject_SpeciesTree
# 4. wait. it might take a long time to run all of the blast searches. then the program will send them to the cluster to align (per gene), wait until they are done checking every... 5 minutes, concatenate the genes for each species, and then send the concatenated file to the cluster to run raxml, and download when it is complete. This might take some time, do not close your terminal or turn your computer off while it is running. subsequent runs with the same gene/will not need to re-do the blast search, so will proceed more quickly.


##this is also usable as a module, when calling from within another program it is probably best to pass a list of trees to make like so [ [spec1, gene1, prefix1], [spec2, gene2, prefix2] ] to get_multiple_concat_alignments()
#& then also run the run_raxml_on_cluster thingy with the output muscle_aligned_list_of_files.

# idea to save time; implement before final run through of giant project
# make sure to save each hit in a current-project file, but ALSO append it to an OVERALL file, and
# when we read a new taxa/gene pair, consider searching the OVERALL file
# before running a new blast DONE

#for local blast use:
#i have nr database downloaded (via blast+ tool '$update_blastdb.pl nr' 'ls *.gz|xargs -n1 tar -xzvf' 'rm *.gz.*') and the download location added to my .profile (export BLASTDB=Users/abigailc/blast)

##########


def Check_if_output_already_exists(the_output_file):
    a = os.system.isfile(the_output_file)
    return a

def remove_slurm_files(ssh_inst, prefix, pattern):
    os.system(ssh_inst+" \' cd Species_Trees;cd "+prefix+"; rm "+pattern+"\'") 
    
class Fasta:
    def __init__(self, name="whatever"):
        # all ids should be stripped and have ">" removed for reasons.
        # for now, sequences do not have any stripping applied
        self.name = name
        self.ids = []
        self.original_ids = []
        self.original_seqs = []
        self.seqs = []
        self.species_list = []
        self.species_names = []
        self.numbers = []
        #the above two are the same thing. don't worry about it. fml.

    def gen_original_lists(self, fastaname):
        try:
            with open(fastaname) as fastafile:
                for line in fastafile:
                    if "\n" == line:
                        pass
                    if ">" in line:
                        # write the previous AA seq
                        try:
                            AAseq = AAseq.strip()
                            self.seqs.append(AAseq)
                            self.original_seqs.append(AAseq)
                        except:
                            pass
                            # initialize a new AAseq
                        AAseq = ""
                        # format the seqID
                        newline = line.strip()
                        newline = newline.strip(">")
                        # write the seqID
                        self.ids.append(newline)
                        self.original_ids.append(newline)
                    else:
                        AAseq = AAseq + line
                AAseq=AAseq.strip()
                # catch the last AAseq pass
                self.seqs.append(AAseq)
                self.original_seqs.append(AAseq)
            #print("Initial sequence and ID lists created. Contains " + str(len(self.ids)) + " sequences")
        except UnboundLocalError:
            print("probably this file :" + fastaname +
                  " has nothing in it. skipping.")
            pass
        except IOError:
            print(os.getcwd())
            print("no file named: " + fastaname +
                  " exists... creating a blank file")
            with open(fastaname, "w") as new:
                pass
            print("hopefully you intended that!")
    def number_seqs(self):
        a = len(self.ids)
        return a
    def gen_species_blast(self):
        self.species_names = []
        #add each thing - the bit between start and | in each id.
        for item in self.ids:
            nam = re.sub("([^\|]*)(\|)(.*)", "\\1", item)
            nam.strip()
            self.species_names.append(nam)
    def gen_raw_blast_lists(self, fastaname):
        bf = open(fastaname, 'r')
        self.species_names = []
        for line in bf:
            gis = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\1", line)
            names = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\3", line)
            seq = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\5", line)
            if "\t" in gis:
                print("ERROR in blast parsing: " + line)
                continue
            else:
                gilist = gis.split(";")
                namelist = names.split("<>")
            for name in namelist:
                index = namelist.index(name)
                gi = gilist[index]
                gi = re.sub("(gi\|)([0-9]*)(.*)", "\\2", gi)
                gi = gi.strip()
                seqid = re.sub("[ ]", "_", name)
                seqid = seqid.strip()
                #test that the identified species is the same as you want it to be... if not, skip and try again.
                #because some dummy likes to name "gene like in [arabidopsis thaliana] [mus musculus]"
                species_as_found_in_gsl = re.sub("([^\[]*)(.*)", "\\2", seqid)
                species_as_found_in_gsl = re.sub("[\[\]]", "",  species_as_found_in_gsl)
                species_name = species_as_found_in_gsl
                if species_name == "":
                    #then we are just dropping this sequence.
                    #not ok!
                    #instead, we will use the gi number to get the proper species name.

                    ###################
                   
                    species_name = gi_to_species_name(gi)
                    if species_name == "error":
                        print("dropping one :"+str(gi))
                        continue
                    else:
                        #print("fixed "+species_name)
                        species_as_found_in_gsl = species_name

                #this will return a species name of form "Danio rerio"
                #lists
                try:
                    species_name = species_name.strip()
                except:
                    print("strip error? dropping: "+str(gi))
                    print(species_name)
                    continue

                species_name = re.sub("[\[\]:;=,/\+'\.\-\(\)", "_", species_name)
                species_name = re.sub(" ", "_", species_name)
                species_name = re.sub("__", "_", species_name)
                species_name = re.sub("__", "_", species_name)
                if species_name[0] == "_":
                    species_name = species_name[1:]
                #keep only one sequence per species names 
                if species_name in self.species_names:
                    pass
               
                else:
                    newID = species_name+"|gi#|"+gi
                    self.ids.append(newID)
                    self.original_ids.append(newID)
                    self.seqs.append(seq)
                    self.original_seqs.append(seq)
                    self.species_names.append(species_name)
    def write_one_seq_per_file(self):
        geneflist = []
        genenames = []
        for i in range(len(self.ids)):
            with open("Seq" + str(i), "w") as new:
                new.write(">" + self.ids[i]+"\n")
                new.write(self.seqs[i]+"\n")
                name = re.sub("([^\|]*)(\|)(.*)", "\\1", self.ids[i])
                geneflist.append("Seq" + str(i))
                genenames.append(name)
        return geneflist, genenames
        print("one per file generated")
    def number_of_sites(self, num = 0):
        try:
            testseq = self.original_seqs[num]
        except:
            print("wtf no sequences in this fasta object?")
            print(self.original_seqs)
            print(num)
        testseq = re.sub("\n", "", testseq)
        #print(testseq)
        #print (len(self.original_seqs[num]))
        #print (len(testseq))
        return len(testseq)
    def shorten(self):
        unk = "no"
        normal = 0
        ucount = 0
        for line in self.ids:
            index = self.ids.index(line)

            # this removes words in brackets that aren't Species_name
            # and then changes NCBI's default naming scheme to be
            #>Species_name|gi#|#########
            # and makes a list of all gi nums and all
            # duplicates
            number = re.sub(
                "(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\3", line)
            num = number.strip()
            edit1 = re.sub(
                "(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\8\\2\\1#|\\3", line)
            if "[" in edit1:
                unk = "no"
                normal += 1
            edit2 = re.sub("[\[\]]", "", edit1)
            edit3 = re.sub("[:;=,/\+'\.\(\)]", "_", edit2)
            edit4 = re.sub(" ", "_", edit3)
            edit4 = re.sub("__", "_", edit4)
            if unk == "no":
                self.ids[index] = edit4
            else:
                print("Unknown Species in ID:" + line)
    def single_shorten(self):
        unk = "no"
        normal = 0
        ucount = 0
        for line in self.ids:
            index = self.ids.index(line)

            # this removes words in brackets that aren't Species_name
            # and then changes NCBI's default naming scheme to be
            #>Species_name|gi#|#########
            # and makes a list of all gi nums and all
            # duplicates
            seplist = line.split("|")
            newid = "gi#|"+seplist[1]
            self.ids[index] = newid
    def blast2fasta(self, blastlist, ENTREZ=False, num=False):
        #returns "bad" if bad. else returns True
        # entrez is used to ensure that sequence saved uses correct TAXON, esp. if sequence is a MULTISPECIES entry.
        # entrex should be somethin like "Mycobacterium triplex"
        # num is how many sequences to write. for species trees, we almost certainly only want one.
        # for converting full downloaded .fastas, we will want all of them (default = False means to do all of them)
        # Converts blast outfmt "6 sseqid stitle sseq" to original lists if
        # entrez = false

        #... now converting outfmt "6 sallseqid salltitles sseq" to sh fasta with selection of proper gi/acc/taxon
        # this should take format " " blast names and replace them with the proper
        # fasta shit
    
        # we open each file in a unique call to blast2fasta. files should be
        # deleted afterwards.
        #print(blastlist)
        #print(os.getcwd())
        bf = open(blastlist, 'r')
        bf2 = open(blastlist, 'r')
        #print(bf)
        error = 0
        end = "no"
        bad = "no"
        length_bf = 0
        for line in bf2:
            length_bf+=1
        errnum = length_bf
        #print("Maxerror: "+str(length_bf))
        run = 0
        #this should be 5 if called from current setup, but might be less if less hits are found.
        for line in bf:
            run +=1
            if end == "yes":
                break
            # gi|738518257|ref|WP_036466735.1|;gi|620038207|emb|CDO87046.1|   50S
            # ribosomal protein L15 [Mycobacterium triplex]<>50S ribosomal protein L15
            # [Mycobacterium triplex]
            #i just removed a ] from group3/... used to be (.*]) but occasionally species name isn't given at end of thing causes errors.
            gis = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\1", line)
            names = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\3", line)
            seq = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\5", line)
            # this removes sequences with no Species_name given, so as to avoid errors
            # downstream
            if "\t" in gis:
                error += 1
                #print(error)
                print("ERROR in blast parsing: " + line)
                continue
            #check if entrez is in the thing
            else:
                #THIS IS GOING TO ERROR A GOOD AMOUNT I THINK. MAYBE AVOID USING THIS
                gilist = gis.split(";")
                namelist = names.split("<>")
                if ENTREZ is False:
                    index = 0
                else:
                    ENTREZ = ENTREZ.strip("\"")
                    found = "no"
                    for item in namelist:
                        if ENTREZ+"]" in item:
                            index = namelist.index(item)
                            found = "yes"
                    if found == "no":
                        #entrez not found in this thing
                        error += 1
                        #print("Name error... might fix")
                        print(namelist)
                        if error == errnum:
                            print("Serious ENTREZ error:")
                            print(ENTREZ)
                            print("This gene wasn't found in this taxon, skipping")
                            bad = "yes"
                            break
                        continue
                        
                        
                try:
                    seqi = gilist[index].strip() + namelist[index].strip()
                except UnboundLocalError:
                    #print(error)
                    #print("Name error... might fix")
                    if error == errnum:
                        print("Serious ENTREZ error:")
                        print(ENTREZ)
                        print(namelist)
                        print("This gene wasn't found in this taxon, skipping")
                        bad = "yes"
                        break
                    continue
                except:
                    #print(index)
                    error += 1
                    print("unknown index error")
                    if error == errnum:
                        break
                    continue
                #print("passed unbound error")
                    # goes to next line, abandoning this one
                
                # strips for .fasta format
                seqid = re.sub("[ ]", "_", seqi)
                seqid = seqid.strip()
                seqid = seqid.strip(">")
                #test that the identified species is the same as you want it to be... if not, skip and try again.
                #because some dummy likes to name "gene like in [arabidopsis thaliana] [mus musculus]"
                species_as_found_in_gsl = re.sub("([^\[]*)(.*)", "\\2", seqid)
                species_as_found_in_gsl = re.sub("[\[\]]", "",  species_as_found_in_gsl)
                species_as_found_in_gsl = re.sub("[_]", " ",  species_as_found_in_gsl)
                if ENTREZ is False:
                    pass
                else:
                    if ENTREZ == species_as_found_in_gsl:
                        pass
                    else:
                        #print(ENTREZ)
                        #print(species_as_found_in_gsl)
                        #print(seqid)
                        print("Odd, species identified did not match Entrez. Trying next seq, might fix.")
                        error += 1
                        #print(error)
                        end = "no"
                        if error == errnum:
                            print("Serious ENTREZ error:")
                            print(ENTREZ)
                            print(namelist)
                            print("This gene not found in this taxon, skipping")
                            print("Trying to label as non-existant to prevent re-searching later... ")
                            bad = "yes"
                            break
                        continue
               
                #print("passed entrez error")
                # add the new sequence id to the list.
                self.ids.append(seqid)
                self.original_ids.append(seqid)
                # the new sequence
                slist = []
                count = 0
                newseq = ""
                for letter in seq:
                    if count > 79:
                        count = 0
                        newseq = newseq + ("\n")
                    newseq = newseq + letter
                    count += 1
                end = "yes"
                #print(newseq)
                #print(seqid)
                self.seqs.append(newseq.strip())
                self.original_seqs.append(newseq.strip())
        #print(error)
        #print(run)
        #how do you manage to get to here, without either being set "Bad" or set good...
        if bad == "yes":
            #print("Has been defined as bad")
            entrez = ENTREZ.replace(" ", "_")
            self.ids.append(entrez+"|gi#|000000000")
            self.original_ids.append(entrez+"|gi#|000000000")
            self.original_seqs.append("----")
            self.seqs.append("----")
            return "bad"
        else:
            #print("returning now")
            return True
    def load_info_swap(self, info_file_in):
        # reads a file of form
        #   originalID
        #   changedID
        # and generates self.ids from that file.
        kid = "no"
        vid = "no"
        CTdict = {}
        with open(infofile) as old:
            for line in old:
                # first pass: gets key (original ID)
                # second pass: gets value (new ID)
                # if we have no info, get key
                if kid == "no":
                    key = line.strip()
                    kid = "yes"
                    continue
                elif kid == "yes":
                    # if we have key and value, record.
                    if vid == "yes":
                        CTdict[key] = value
                        vid = "no"
                        kid = "no"
                        continue
                    # if we have key but no value, get value.
                    if vid == "no":
                        value = line.strip()
                        vid = "yes"
            # catch the final pass
            CTdict[key] = value
        for item in self.original_ids:
            index = self.original_ids.index(item)
            newid = CTdict[item]
            self.ids[index] = newestid
        # done
        # troubleshooting: do not preform this operation after any that change
        # self.ids. this op must be done first, or in a seperate command.
    def gen_new_fasta(self, new_fasta_name):
        # this should print the changed seqids and changed AA sequences to
        # file.
        newfasta = new_fasta_name
        # print(len(self.original_ids))
        # print(len(self.ids))
        # print(len(self.original_seqs))
        # print(len(self.seqs))
        
        with open(newfasta, "w") as new:
            for i in range(len(self.original_ids)):
                new.write(">" + self.ids[i].strip() + "\n")
                # print(i)  #
                # unclear if this needs a "\n" after it... check.#TODO
                new.write(self.seqs[i].strip()+"\n")
        #print("Finished, your new fasta file is located at :" + newfasta)
        # done
        return newfasta
    def swap_in_newick(self, old_newick_name, new_file_name):
        # this replaces the tip names in a newick file. sometimes works on nexus
        # files too, but I havent extensively tested it.
        newick = old_newick_name
        newnewick = new_file_name
        with open(newick) as old:
            with open(newnewick, "w") as new:
                for line in old:
                    for item in self.original_ids:
                        index = self.original_ids.index(item)
                        line = line.replace(item, self.ids[index])
                        new.write(line)
        print("finished, tip-replaced-newick file at: " + newnewick)
        # done
    def swap_in_nexus(self):
        print (
            "You didn't implement this yet. try using newick replace, it might work")
        pass
        # something
        # to-do, try nexus replace in the meantime, it should work
    def gen_info(self, info_file_name):
        # writes a file of form
        #   originalID
        #   changedID
        with open(info_file_name, "w") as inf:
            listlength = len(self.original_ids)
            if listlength != len(self.ids):
                print ("List lengths do not match! FATAL ERROR")
                print (self.original_ids)
                print (self.ids)
                raiseSystemExit
            for i in range(listlength):
                inf.write(self.original_ids[i])
                inf.write(self.ids[i] + "\n")
        print("Info file was generated. Named " + info_file_name)
        # done
# this ought to be deleted - but needs to be tested... hopefully it is identicle to the version in Classes_Standalone
    def gen_species_lists(self):
        speclist = []
        for item in self.ids:
            taxon = re.sub("([^_]*)([A-Z][a-z]*_[A-Z]?[a-z]*[^\|]*)(.*)", "\\2", item)
            if "|" in taxon:
                tlist = item.split("|")
                taxon = tlist[-1]
                if "|" in taxon:
                    print ("TAXON error in gen_species_lists():" + taxon)
            speclist.append(taxon)
        self.species_list = speclist
        self.species_names = speclist

        return speclist
    def SetTaxID(self):
  
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
        if len(self.taxid) == len(self.ids):
            pass
        else:
            print("taxid list is not the same length as overall id list")
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
            first = True
            for staff in staffs:
                    tid = staff.getElementsByTagName("TaxId")[0]
                   
                    tname = staff.getElementsByTagName("ScientificName")[0]
                    trank = staff.getElementsByTagName("Rank")[0]
                    taxid = tid.firstChild.data
                    taxname = tname.firstChild.data
                    taxrank = trank.firstChild.data
                    if first is True:
                        first = False
                        taxname = taxname.replace(".", "_")
                        taxname = taxname.replace("-", "_")
                        taxname = taxname.replace("=", "_")
                        taxname = taxname.replace(",", "_")
                        taxname = taxname.replace(" ", "_")
                        taxname = taxname.replace("__", "_")
                        taxname = taxname.replace("__", "_")
                        taxdict["ScientificName"] = taxname
                        
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



# this will correlated IDs based on Genus_species across datasets.
def correlate_ids(list_of_id_lists):
    # so now we have a list of species from each Fasta, and we need to correlate them. should be it's own matching function?
    # go though first list, save name (index list1 ,index list2, indexlist3), and add name to "used" list
    # then go through second list, save same thing but skip any that are already in "used" list.
    # continue etc.
    numlists = len(list_of_id_lists)
    output_list = []
    used = []
    for item in list_of_id_lists:
        # list1
        for ids in item:

            # id1 of list1
            name = ids
            if name in used:
                pass
            else:
                used.append(name)
                id_index_list = []
                id_index_list.append(ids)
                # check the index of that id in each list and append to "id_index_list"
                # if the id is not in that list, should append "NA"
                for eachlist in list_of_id_lists:
                    try:
                        index = eachlist.index(ids)
                    except ValueError:
                        index = "NA"
                    id_index_list.append(index)
                # add the result of scanning that id to overall output list.
                output_list.append(id_index_list)
    # output list looks like:
    # outputlist = [ ["Cat",1,2,3,4] , ["Dog", 2,1,13,14] ]
    print("Concat")
    print(output_list)
    return output_list


# this should work but hasn't been tested. requires stuff above.
def Concatenate(listoffastafiles, new_cc_fasta_name, prefix, strain_limit=False):
    # for each thing in listoffastas, create a fasta object in a list.
    fasta_class_list = []
    list_species_lists = []
    # this will create an instance of Fasta from each fasta file
    for i in listoffastafiles:
        fasta_class_list.append(Fasta(i))
    # this will create the original lists (seqID and sequence) for each Fasta object
    # this will create a list of the "Genus_species" from each seq in each
    # Fasta object (pass true/false of strain control if need be)
    for f in fasta_class_list:
        f.gen_original_lists(f.name)
        species_list = f.gen_species_lists()
        list_species_lists.append(species_list)
    # so now we have a list of species from each Fasta, and we need to correlate them. should be it's own matching function?
    # go though first list, save name (index list1 ,index list2, indexlist3), and add name to "used" list
    # then go through second list, save same thing but skip any that are already in "used" list.
    # continue etc.
    # check that things i already did for greg... swapping placement in file
    # was similar.. where did i save that.
    indexed_ids_list = correlate_ids(list_species_lists)
    # indexed_ids_list in form [ ["Cat", 1,2,3],["dog",2,1,2] ]
    # this part will create a new, concatenated .fasta file
    # requires indexed_ids_list, fasta_class_list, function that returns
    # number of sites. self.number_of_sites
    #debugging:
    
    lenlist_final = []
     #ensure that for each fas in fasta class list, all sequences are of the same length.
    for fas in fasta_class_list:
        ndash1 = fas.number_of_sites()
    with open(new_cc_fasta_name, "w") as new:
        for item in indexed_ids_list:
            lenlist = []
            # do i want to print species_name to file, or full
            # taxonomy? for now, just species_name is gunna happen,
            # but easy to swap i think
            new.write(">" + item[0].strip()+"\n")
            fas_num = 0
            allseq = ""
            #ensure that for each fas in fasta class list, all sequences are of the same length.
            for fas in fasta_class_list:
                fas_num += 1
                # fas_num keeps track of what number fasta we are on, which
                # correlates to the index of the index in indexed_ids_list
                search_index = item[fas_num]
                # search_index will be something like "22"

                # if search_index is NA, generate a str of "-" that is n
                # characters long, where n is the return from
                # fas.number_of_sites
                if search_index == "NA":
                    ndash = fas.number_of_sites()
                    retreived_seq = ""
                    for i in range(int(ndash)):
                        retreived_seq = retreived_seq + ("-")
                else:
                    retreived_seq = fas.seqs[search_index]
                    # retreived_seq wil be something like "the 22nd sequence in
                    # object Fas's sequence list... " or "BLAHSEQUENCEDATA"
                retreived_seq = re.sub("\n", "", retreived_seq)
                lenlist.append(len(retreived_seq))
                count = 0
                allseq = allseq + retreived_seq
            #print(lenlist)
            #print(len(allseq))
            lenlist_final.append(len(allseq))
            newseq = ""
            for letter in allseq:
                if count > 79:
                    count = 0
                    newseq = newseq + ("\n")
                newseq = newseq + letter
                count += 1
                
            new.write(newseq.strip()+"\n")
                # there might be too many line-breaks in this fashion. if so,
                # concatenated all retreived_seqs and then only call write()
                # once.
        for length in lenlist_final:
            if length == lenlist_final[0]:
                pass
            else:
                print("ERROR your concat sequences are not all of the same length something has gone horribly wrong! aborting.")
                raise SystemExit
        print("Should be finished generating new concatenated fasta at: " + new_cc_fasta_name)
    print("done w cc gen!")
    return new_cc_fasta_name

def muscle_align_on_cluster(end_file_list, prefix):
    #this creates dir you will use on the cluster.
    aligned_list = []
    remove_list = []
    for item in end_file_list:
        #print("testing if enough data to continue...")
        testf = Fasta(item)
        testf.gen_original_lists(item)
        a = testf.number_seqs()
        #print("fasta has "+str(a)+"sequences....")
        if a < 5:
            #print(item+" has too few sequences to be aligned and run through raxml. removing from analysis")
            #if len(end_file_list) is 1:
            #    print("aborting job... no gene has enough sequences!")
            #    raise SystemExit
            remove_list.append(item)
        else:
            aligned_list.append(item+"_Muscle.fasta")
    for athing in remove_list:
        end_file_list.remove(athing)
    if end_file_list == []:
        print("Too few files for to make tree")
        return "NA"
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Species_Trees"
    a = gen_muscle_script(prefix+"_Sc.sh", "~"+clus_path+"/"+prefix+"/"+prefix+"_Corr.txt", str(len(end_file_list)), prefix+"job")
    b = gen_correlate_file(end_file_list, prefix+"_Corr.txt")
    end_file_list.append(a)
    end_file_list.append(b)
    direct = os.getcwd()
    move_to_cluster(end_file_list, clus_path+"/"+prefix)
    print("everything should be generated and on the cluster... beginning species_tree_muscle alignment")
    os.system(ssh_inst+" 'cd ~/Species_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    movehome = []
    for i in aligned_list:
        movehome.append(i)
    anum = len(movehome)
    if anum < 15:
        timen = 60
        times = "1"
    elif anum <25:
        timen = 120
        times = "2"
    else:
        times = "3"
        timen = 180
    time.sleep(30)
    
    while finished is not True:
      
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename[2:]+" "+direct+"/"+prefix)
        finished = "yes"
        #and will remain so if everything is home this time.
        for item in aligned_list:
            if item in movehome:
            #see if it got moved home.
                exists = os.path.isfile(item)
                #see what size it is (sometimes, we move home in the split second between file creation and writing).
                if exists is True:
                    size = os.path.getsize(item)
                    if size > 0:
                        movehome.remove(item)
                    else:
                        finished = False
                        print("waiting some time "+times+" mins and trying again")
                        print("didn't find: "+item)
                else:
                    finished = False
            else:
                pass
        if len(movehome) > 0:
            print("still to move home")
            print(movehome)
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait some minutes and then try again.
            time.sleep(60)
            print("test")
            print("Rax not done yet. could not locate : "+item+"checking again in 5 minutes")
            finished = "yes"
    c = Concatenate(aligned_list, prefix+"_CC.fasta", prefix)
    #print(c)
    print("muscle alignments finished. _CC files should exist locally")
    print("removing .fasta, muscle, and slurm.out files from cluster...")
    remove_slurm_files(ssh_inst, prefix, "slurm*")
    remove_slurm_files(ssh_inst, prefix, "*.fasta*")
    #if list_of_roots == False:
    return c
  #  else:
   #     return c, list_of_roots


def write_species_list(species_list, prefix):
    with open("./"+prefix+"/"+prefix+"_Species_List.txt", "w") as new:
        for item in species_list:
            item = item.strip()
            item = item.strip("\"")
            new.write(item+"\n")








 



def Set_Fasta_Species_List_And_Category(list_of_big_clades, Euk_Gene_List, Bac_Gene_List, Arc_Gene_List, All_Gene_List):
    sp_tree_gen_list = []
    for subtree_ob in list_of_big_clades:
        #set the species list
        b = subtree_ob.set_fasta_object()
        a =  subtree_ob.ret_fasta_object()
        species_li = a.gen_species_lists()
        species_list = []
        for item in species_li:
            item = item.strip()
            item = "\"" + item + "\""
            species_list.append(item)
        subtree_ob.species_list_original = species_li 
        subtree_ob.species_list_gene_to_species = species_list
        #now set the category -> genes file.
        #get tip list
        #see if string in first tip
        #if so, see if Arc/Euk/Bac also in.
        strcheck = subtree_ob.ret_string()
        tiplist = subtree_ob.ret_tips()
        for item in tiplist:
            if strcheck in item:
                if "Arc|" in item:
                    cat = "Archaea"
                    break
                if "Archaea|" in item:
                    cat = "Archaea"
                    break
                if "Bac|" in item:
                    cat = "Bacteria"
                    break
                if "Bacteria|" in item:
                    cat = "Bacteria"
                    break
                if "Euk|" in item:
                    cat = "Eukaryota"
                    break
                if "Eukaryota|" in item:
                    cat = "Eukaryota"
                    break
              
                else:
                    break
        try:
            if cat == "":
                print("Error with clade: "+str(subtree_ob.ret_string()))
        except:
            cat = "All"
        else:
            subtree_ob.set_category(cat)
        #pick the file based on category to send on to species tree generator
        if cat == "Eukaryota":
            cat_file = Euk_Gene_List
        if cat == "Archaea":
            cat_file = Arc_Gene_List
        if cat == "Bacteria":
            cat_file = Bac_Gene_List
        if cat == "All":
            cat_file = All_Gene_List
        print("set :"+subtree_ob.ret_string()+" to category file: "+cat_file)
        subtree_ob.cat_file = cat_file

def Fix_Your_Mistakes(list_of_big_clades):
    #TO FIX A CURRENT AND TEMPORARY NAMING ERROR
    for item in list_of_big_clades:
        cladestr = item.string_name
        if "myc" in cladestr[-3:]:          
            cladestr = re.sub("myc", "mycetes", cladestr)
        
        if "mycl" in cladestr[-4:]:          
            cladestr = re.sub("mycl", "mycetales", cladestr)
        
        if "bacl" in cladestr[-4:]:
            cladestr = re.sub("bacl", "bacteriales", cladestr)
        
        elif "bac" in cladestr[-3:]:
            cladestr = re.sub("bac", "bacteria", cladestr)
       
        item.string_name = cladestr
        #done




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

#stuff for -b
def gi_to_species_name(gi):
    #this takes a gi number and gets from it a species name, by first getting a taxid and then taxid -> name + formatting.
    gi = gi.strip()
    taxid = gi_to_taxid(gi)
    if taxid == "error":
        return "error"
    #else:
    #    print (taxid)
    scientific_name = taxid_to_scinam(taxid)

    return scientific_name

def taxid_to_scinam(taxid):
 

    #####new bit added nov 16: no longer requires XMLLINT, does require minidom, urllib2 (probably are standard)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id='+taxid # define XML location w. current taxid
    try:
        dom = minidom.parse(urllib2.urlopen(url)) # parse the data
    except:
        print("error with url:"+ url)
        try:
            time.sleep(10)
            dom = minidom.parse(urllib2.urlopen(url)) 
        except:
            print("double couldnt parse url: "+url)
            return "error"

    staffs = dom.getElementsByTagName("Taxon")
            #print (staffs)
    first = True
    for staff in staffs:
            tid = staff.getElementsByTagName("TaxId")[0]     
            tname = staff.getElementsByTagName("ScientificName")[0]
            #trank = staff.getElementsByTagName("Rank")[0]
            taxid = tid.firstChild.data
            taxname = tname.firstChild.data
            #taxrank = trank.firstChild.data
            if first is True:
                first = False
                taxname = taxname.replace(".", "_")
                taxname = taxname.replace("-", "_")
                taxname = taxname.replace("=", "_")
                taxname = taxname.replace(",", "_")
                taxname = taxname.replace(" ", "_")
                taxname = taxname.replace("__", "_")
                taxname = taxname.replace("__", "_")
                if taxname == "":
                    return "error"
                return taxname
    #should return scientific name to you
           
def gi_to_taxid(gi):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+gi+"&retmode=xml" # define XML location w. current taxid
    try:
        dom = minidom.parse(urllib2.urlopen(url)) # parse the data
    except:
        try:
            time.sleep(10)
            dom = minidom.parse(urllib2.urlopen(url)) 
        except:
            print("couldnt parse url: "+url)
            return "error"
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
        print("Error with "+url)
        return "error"
    return taxonid

#-b
def Single_Blast_Search(QUERY, ENTREZ, STORED_BLAST_FILE, MYFASTA):
    #this does a blast for each missing bit. might take some time.

    #WHEN YOU HAVE THE SEQ DO THE ABOVE ^

    print("Searching")
    OUTPUT = "blastp_intermediate.fasta"

        # run the blast
        #do the thing remotely
    blast_query = "blastp -remote -query " + QUERY + " -db nr -out " + OUTPUT +" -max_target_seqs 10 -entrez_query " + ENTREZ + " -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\""
                #print(blast_query)
    #probably the loss in information is do to malformed titles.
    #so we need to "fix" it by getting the species_name from the given gi or accession numbers.
    #so lets see how the shortening goes? and the call Get Taxonomy.
    os.system(blast_query)
    #a = raw_input("does blastp_intermediate exist?")
    #temporary to deal ith DNS issues / pausing for canceling.
    tempfasta = Fasta(OUTPUT)


    a = tempfasta.blast2fasta(OUTPUT)

    tempfasta.gen_original_lists(OUTPUT)
    try:
        a = tempfasta.ids[0]
    except:
        print("couldn't make an original list from result")
        return "None"
    tempfasta.gen_new_fasta("postb2ftest.txt")
    #print(tempfasta.ids[0])
    a = tempfasta.single_shorten()
    #print(tempfasta.ids[0])
    gis = tempfasta.gen_numbers()
    #print(tempfasta.numbers[0])
    taxids = tempfasta.SetTaxID()
    #print(tempfasta.taxid[0])
    rankdict = tempfasta.GetTaxonomy()
    #print(tempfasta.taxonomy[0])
    withtaxonomy = tempfasta.AppendTaxonomy(["ScientificName"])
    print(tempfasta.ids[0])
    a = tempfasta.gen_species_lists()
    #print(tempfasta.species_list[0])
    found = False
    for i in range(len(tempfasta.original_ids)):
        print(ENTREZ)
        print(tempfasta.ids[i])
        if ENTREZ == tempfasta.species_names[i]:
            found = True
            #add this sequence and id and then exit
            MYFASTA.original_seqs.append(tempfasta.seqs[i])
            MYFASTA.original_ids.append(tempfasta.ids[i])
            MYFASTA.seqs.append(tempfasta.seqs[i])
            MYFASTA.ids.append(tempfasta.ids[i])
            STORED_BLAST_FILE.original_seqs.append(tempfasta.seqs[i])
            STORED_BLAST_FILE.original_ids.append(tempfasta.ids[i])
            STORED_BLAST_FILE.seqs.append(tempfasta.seqs[i])
            STORED_BLAST_FILE.ids.append(tempfasta.ids[i])

            return "Done"
    if found is False:
        return "None"

    #what is this outptu

def Gen_New_Species_list():
    for item in list_of_big_clades:
        self.species_fasta_object
    
def gen_muscle_script(scriptfile, indexname, n, Jobname):
    #currently assuming you are running in the dir that files are in and should be returned to.
    # direct = os.getcwd()
    # host = socket.gethostname()
    # user = getpass.getuser()
    # addr = user+"@"+host+":"+direct
    #figure out how many need to be run. n = len(listoffilesinindex). n
    #figure out a name - scriptfile
    #figure out path to index file. indexname
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    addr = "PLACEHOLDER"
    

    ##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 0-10:00:00                                                                                   
#SBATCH -J Mus"""+Jobname+"""                                                                                         
##SBATCH -o Mus"""+Jobname+""".out
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/muscle/3.8.31
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
ENDING=_Muscle.fasta
echo $THE_INPUT_FILE$ENDING

muscle -in $THE_INPUT_FILE -out $THE_INPUT_FILE$ENDING

##scp $THE_INPUT_FILE$ENDING """+addr+"""

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

def gen_raxml_script(scriptfile, indexname, n, Jobname, list_of_root_tips = False, number_nonstandard_bootstraps = 100):
    #is jobname == prefix?
    #add second /prefix/ to indexname.
    #currently assuming you are running in the dir that files are in and should be returned to.
    # direct = os.getcwd()
    # host = socket.gethostname()
    # user = getpass.getuser()
    # addr = user+"@"+host+":"+direct
    #figure out how many need to be run. n = len(listoffilesinindex). n
    #figure out a name - scriptfile
    #figure out path to index file. indexname
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    if list_of_root_tips != False:
            
        a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 7-00:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
# #SBATCH --exclusive                                                                                   
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
THE_OUTGROUP_TIP=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $3}' )
echo $THE_INPUT_FILE
ENDING=_Muscle.fasta
echo $THE_INPUT_FILE$ENDING

NEW=${THE_INPUT_FILE%%.*}
echo $NEW
  
raxmlHPC-PTHREADS-AVX -T 20 -o $THE_OUTGROUP_TIP -f a -m PROTGAMMALGF -p 12345 -x 12345 -#"""+str(number_nonstandard_bootstraps)+""" -n $NEW -s $THE_INPUT_FILE         

exit"""
    else:

##example script
        a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 7-00:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
# #SBATCH --exclusive                                                                                   
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
ENDING=_Muscle.fasta
echo $THE_INPUT_FILE$ENDING

NEW=${THE_INPUT_FILE%%.*}
echo $NEW
  
raxmlHPC-PTHREADS-AVX -T 20 -f a -m PROTGAMMALGF -p 12345 -x 12345 -#"""+str(number_nonstandard_bootstraps)+""" -n $NEW -s $THE_INPUT_FILE         

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

def gen_correlate_file(list_of_input_files, corr_file, list_of_root_tips=False):
    #this should be in form
    #1 name
    #2 name
    #3 name
    #requires 1: list of files 2. name for corr_file.
    i = 0
    with open(corr_file, "w") as corr:
        for item in list_of_input_files:
            if "/" in item:
                itemlist = item.split("/")
                item = itemlist[-1]
            i += 1
            #make sure the \n spacing works correctly.
            if list_of_root_tips != False:
                corr.write(str(i)+" "+item+" "+list_of_root_tips[(i-1)]+"\n")
            else:
                corr.write(str(i)+" "+item+"\n")
    return corr_file

def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)
    
def write_species_list(species_list, prefix):
    with open("./"+prefix+"/"+prefix+"_Species_List.txt", "w") as new:
        for item in species_list:
            item = item.strip()
            item = item.strip("\"")
            new.write(item+"\n")

#13
def raxml_run_on_cluster(cc_file_list, prefix, list_of_roots = False, number_nonstandard_bootstraps = 100):
    #this creates dir you will use on the cluster.
    tree_list = []
    besttree_list = []
    to_remove = []
    remove_indices = []
    print("potentials to align:")
    print(cc_file_list)
    ind = 0
    for finam in cc_file_list:
        if finam is "NA":
            print(finam+" had too few sequences to be aligned. removing from analysis")
            remove_indices.append(ind)
            print(remove_indices)
            ind+=1
            continue
        #remove any that will error due to too few sequences
        ind +=1
        Fas = Fasta(finam)
        Fas.gen_original_lists(finam)
        numseq = Fas.number_seqs()
        print(finam+" has sequences = "+str(numseq))
        if numseq >3:
            alist = finam.split(".")
            athing = alist[0]
            tree_list.append("RAxML_bipartitions."+athing)
            besttree_list.append("RAxML_bestTree."+athing)
        else:
            
            print(finam+" has too few sequences to be aligned and run through raxml. removing from analysis")
            #            to_remove.append(finam)

            remove_indices.append(cc_file_list.index(finam))
    
    remove_indices.sort()
    
    remove_indices.reverse()
    for itemr in remove_indices:
        #why does this think that list_of_roots is a string?
        #how to avoid error if we do not pass in a list of roots? FALSE cannot be removed either.
        cc_file_list.pop(itemr)
        if list_of_roots != False:
            list_of_roots.pop(itemr)
        #if len(cc_file_list) is 1:
         #   print("aborting job...")
          #  raise SystemExit
        #this is going to be something different... like raxml_bipartitions.blah
        #print(tree_list)
    if cc_file_list == []:
        return []
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Species_Trees"
    if list_of_roots == False:
        a = gen_raxml_script(prefix+"_Rax_Sc.sh", "~"+clus_path+"/"+prefix+"/"+prefix+"_Rax_Corr.txt", str(len(cc_file_list)), prefix+"job", False, number_nonstandard_bootstraps)
        b = gen_correlate_file(cc_file_list, prefix+"_Rax_Corr.txt")
    else:
        list_roots_2 = []
        for root_tip in list_of_roots:
            root_tip = root_tip.strip("\"")
            root_tip = root_tip.replace(" ", "_")
            if " " in root_tip:
                print("still space!")
                re.sub(" ", "_", root_tip)
                print(root_tip)
            list_roots_2.append(root_tip)
        a = gen_raxml_script(prefix+"_Rax_Sc.sh", "~"+clus_path+"/"+prefix+"/"+prefix+"_Rax_Corr.txt", str(len(cc_file_list)), prefix+"job", list_roots_2, number_nonstandard_bootstraps)
        b = gen_correlate_file(cc_file_list, prefix+"_Rax_Corr.txt", list_roots_2)
    cc_file_list.append(a)
    cc_file_list.append(b)
    direct = os.getcwd()
    move_to_cluster(cc_file_list, clus_path+"/"+prefix)
    print("everything should be generated and on the cluster. starting raxml.")
    os.system(ssh_inst+" 'cd ~/Species_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    movehome = []
    ret_list = []
    bothlist = []
    for i in tree_list:
        movehome.append(i)
        bothlist.append(i)
    for i in besttree_list:
        bothlist.append(i)
        movehome.append(i)
    time.sleep(20)
    while finished is not True:
        finished = "yes"
        #and will remain so unless otherwise is true
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+prefix+"/"+filename+" "+direct)
        for item in bothlist:
            #see if it got moved home.
            exists = os.path.isfile("./"+item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)                
                rax, ext = item.split(".")
                os.system("cp "+item+" "+ext+"_Rax_Bipart.newick")
                ret_list.append(item)
            else:
                finished = False
                print("Rax not done yet. could not locate :./"+item+"checking again in 5 minutes")
        if len(movehome) > 0:
            print("all items still not found at home:")
            print(movehome)
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait ten minutes and then try again.
            time.sleep(300)
            finished = "yes"
    c = tree_list
    #print(c)
    print("RAXML finished, tree(s) should exist locally")
    print("removing .fasta, muscle, and slurm.out files from cluster...")
    remove_slurm_files(ssh_inst, prefix, "slurm*")
    remove_slurm_files(ssh_inst, prefix, "*.fasta*")
    #if list_of_roots == False:
    return c
  #  else:
   #     return c, list_of_roots

#12.1
def convert_to_genus_species(sciname):
    ilist = sciname.split("_")
    try:
        item = ilist[0]+"_"+ilist[1]
        return item
    except:
        return sciname

#12
def Get_Sequences(clade_list, projectname, fixmissing = False, strain = False):
    #also call the missing thing within here
    cc_files_out = []
    for clade in clade_list:
        #in this version will be only one -- your given fasta.
        end_file_list = []
        prefix = clade.prefix
        gene_name_list = clade.gene_names_list
        #make the new files. these will contain all the species from the gene fasta in all the new genes specified, aligned, and then cc'd
        i = 0
        list_overall_gene_fastas = clade.blast_stored
        list_raw_fasta = clade.blast_raw
        #gene 1, 2, 3, 4 etc. the ribos
        for QUERY in gene_name_list:
            query_seq_name = "Seq"+str(i)
            print(QUERY)
            #for each gene you want in the overall list....
            current_overall_fasta = list_overall_gene_fastas[i]
            if current_overall_fasta.species_names == []:
                current_overall_fasta.gen_species_lists()   
            #print("Current overall: "+current_overall_fasta)
            current_raw_fasta = list_raw_fasta[i]
            if current_raw_fasta.species_names == []:
                current_raw_fasta.gen_species_lists()
            #print("current raw: "+current_raw_fasta)
            NewFasta = Fasta()
            NewFasta.projectname = projectname
            #set up the new fasta to write to (NewFasta)
            if prefix is not False:
                NewFastaName = "./"+prefix+"/"+prefix + gene_name_list[i].strip() + ".fasta"
            else:
                NewFastaName = "./"+prefix+"/"+gene_name_list[i].strip() + ".fasta"
            i += 1
            donum = str(i)
            totalnum = str(len(gene_name_list))
            print("Beginning on: " + NewFastaName + "which is number: " +donum + " of a total :" + totalnum)
            #this is where were getting into the species level... has gone too far.
            #for each species....
            #print(clade.species_list_plus_og_loss)

            list_of_taxa_to_get = clade.species_list_plus_og_loss
            if strain is True:
                templist = []
                #convert this list to only have genus_species.
                for item in list_of_taxa_to_get:
                    a = convert_to_genus_species(item)
                    templist.append(a)
                list_of_taxa_to_get = templist

            print("first taxa to get: "+list_of_taxa_to_get[0])
            for ENTREZ in list_of_taxa_to_get:  
                in_overall = "no"
                no_under = re.sub(" ", "_", ENTREZ.strip("\""))     
                #see if its already listed in current_overall_fasta
                for iteration in range(len(current_overall_fasta.species_names)):
                    item = current_overall_fasta.species_names[iteration]
             
                    if strain is True:
                        item = convert_to_genus_species(item)
                    #print(item)
                    if no_under == item:
                        in_overall = "yes"
                        
                        nseq = current_overall_fasta.seqs[iteration]
                        item = item.strip()
                        nseq = nseq.lstrip()
                        if nseq == "----":
                        #was searched before and got no valid hits
                            break
                        else:
                            NewFasta.original_seqs.append(nseq)
                            NewFasta.original_ids.append(item)
                            NewFasta.seqs.append(nseq)
                            NewFasta.ids.append(item)
                            break
                if in_overall == "yes":
                    continue
                #if not in current overall, look in raw_blast
                else:
                    match = False
                    #see if its in current_raw_fasta
                    for iteration in range(len(current_raw_fasta.species_names)):
                        spec = current_raw_fasta.species_names[iteration]
                   
                        if strain is True:
                            spec = convert_to_genus_species(spec)
                        if no_under == spec:
                            #if no_under != spec:
                                #print("subbing "+spec+" for query: "+no_under+" ... sanity check?")
                            match = True
                            
                            nseq = current_raw_fasta.seqs[iteration]
                            nid = current_raw_fasta.ids[iteration]
                            nseq = nseq.strip()
                            nid = nid.strip()
                            NewFasta.original_seqs.append(nseq)
                            NewFasta.original_ids.append(no_under)
                            NewFasta.seqs.append(nseq)
                            NewFasta.ids.append(no_under)

                            current_overall_fasta.original_seqs.append(nseq)
                            current_overall_fasta.original_ids.append(nid)
                            current_overall_fasta.seqs.append(nseq)
                            current_overall_fasta.ids.append(nid)
                            break
                    #temp
                    #match = False
                    if match == False:
                        if fixmissing == False:
                            pass
                        else:
                            print("Could not find: "+no_under+" in existing files. going to try an independent query")
                            ret = Single_Blast_Search(query_seq_name, no_under, current_overall_fasta, NewFasta)
                            if ret == "None":
                                print("No sequence for "+no_under)
                            elif ret == "Done":
                                print("....fixed!")

            current_overall_fasta.shorten()
            #add any newly-found sequences to the list of stored current overall's
            current_overall_fasta.gen_new_fasta(current_overall_fasta.name)
            NewFasta.shorten()
            end_fas_num = NewFasta.number_seqs()
            #print(end_fas_num)
            #print(NewFasta.ids)
            if end_fas_num > 2:
                end_file_list.append(NewFasta.gen_new_fasta(NewFastaName))
            else:
                print("Too few sequences found in fasta associated with "+QUERY+", this gene will be excluded")
        print("Should be finished creating one fasta file for each gene in your input! They are called: ")
        #print(end_file_list)
        #for removethis in gene_f_list:
         #   os.system("rm " + removethis)
        # after everything else...
        for d in list_overall_gene_fastas:
            nam = d.name
            d.gen_new_fasta(nam)
        a = muscle_align_on_cluster(end_file_list, prefix)
        #print (a)
        #a will be a concatenated file, including all species specified w/ all 30 genes aligned independently and than concatenated.
        if a == "NA":
            clade.cc_file = "NA"
            pass
        else:
            cc_files_out.append(a)
            clade.cc_file = a
            clade.cc_file_obj = Fasta(a)
            clade.cc_file_obj.gen_original_lists(a)
            cc_sp_list = clade.cc_file_obj.gen_species_lists()
            #send the species_names from the new concatenated file (cc_sp_list) and the list of taxa we were searching for (list_of_taxa_to_get)
            #compares and prints discrepancy (and saves) 
            compare_the_species_tree_taxa_with_the_original_taxa(cc_sp_list, list_of_taxa_to_get, projectname)
    return cc_files_out

#11
def Assign_Blasts_To_Clades(list_of_clade, tuple_blast_objects):
    for item in list_of_clade:
        item.blast_raw = "NONE"
        genefile = item.cat_file.strip()
        for tup in tuple_blast_objects:
            if tup[0] == genefile:
                item.blast_raw = tup[2]
                item.blast_stored = tup[1]
                item.gene_names_list = tup[3]
        if item.blast_raw == "NONE":
            print("ERROR in Assign_Blasts_To_Clades.")
            print(genefile)
            print(tuple_blast_objects)
    print("all clades were assigned a blast object")


#10
def Create_Blast_Objects(list_of_clades, genefile_stored_raw_tuples):
    #generate the blasts in Fasta form with names?
    tuple_blast_objects = []
    for item in genefile_stored_raw_tuples:
        name_of_genefile = item[0]
        list_of_stored_filenames = item[1]
        list_of_raw_filenames = item[2]
        list_of_stored_fasobjects = []
        for name in list_of_stored_filenames:
            obj = Fasta(name)
            obj.gen_original_lists(name)
            list_of_stored_fasobjects.append(obj)
        list_of_raw_fasobjects = []
        for name in list_of_raw_filenames:
            #print("looking for raw fasta of this name")
            #print(name)
            obj = Fasta(name)
            obj.gen_original_lists(name)
            obj.gen_species_blast()
            list_of_raw_fasobjects.append(obj)
        list_gene_names = item[3]
        tuple_blast_objects.append( (name_of_genefile, list_of_stored_fasobjects, list_of_raw_fasobjects, list_gene_names) )
    return tuple_blast_objects

#9 Do The Blast    
def Check_Blasts_Exist(clade_align):
    #probably i need to output something like --- which blast file lists need to be looked at per category.
    #make directory
    list_of_genelists = []
    os.system("mkdir Stored_Blasts")
    #Put species list in dir for future reference
    #write_species_list(species_list, prefix)
    for item in clade_align:
        prefix = item.prefix
        os.system("mkdir "+prefix)
        write_species_list(item.species_list_original, prefix)
        genelist = item.cat_file
        if genelist in list_of_genelists:
            pass
        else:
            list_of_genelists.append(genelist)
    print(len(list_of_genelists))
    #genelist right now will be three - arc, bac, euk.
    #second pass will be only bac
    #actual verification
    track_stored = []
    #^is a list of lists. the inner lists will be the names of each stored_blast eg [[euk1, euk2, euk3], [bac1,bac2.bac3]]
    track_raw = []
    #^is a list of lists. the inner lists will be the names of each stored_blast eg [[euk1, euk2, euk3], [bac1,bac2.bac3]]
    list_of_gene_name_lists = []
    for gene_sequences_file in list_of_genelists:
        #make sure gene_seq_file exists
        if os.path.isfile(gene_sequences_file) is False:
            print("Error, cannot find specified gene_sequences_file:"+gene_sequences_file)
            raise SystemExit
        
        #make list of stored blasts for this gene.
        GeneFasta = Fasta(gene_sequences_file)
        GeneFasta.gen_original_lists(gene_sequences_file)
        gene_fasta_list, gene_name_list = GeneFasta.write_one_seq_per_file()
        list_overall_gene_fastas = []
        list_raw_gene_fastas = []
        gene_overall_fasta_names = []
        for b in range(len(gene_name_list)):
            gene_name = gene_name_list[b].strip()
            gene_overall_fasta = "./Stored_Blasts/"+gene_name+"_All.fasta"
            gene_overall_fasta_names.append(gene_overall_fasta)
        track_stored.append(gene_overall_fasta_names)
        n = 0
        i = 0
        output_list = []
        # first, initialize (or check existance of) files containing all genes we care about.
        #similar to listoverallgenefastas but instead of containing things we've used, they contain raw blast data.
        #1. check existance of blast-result-files
        #2. if no, do the blast
        #3. now, check if each query/species pair are in the overall files. yes - grab them. no - trawl through the raw blast data.
        #this is up-to-date...
        
        #make dir
        os.system("mkdir Raw_Blasts")
        #set up lists
        raw_fasta_list = []
        gene_raw_fasta_names = []
        #get names of raw_blasts (all)
        gene_raw_blast_conv_names = []
        for d in range(len(gene_name_list)):
            gene_name = gene_name_list[d].strip()
            gene_raw_fasta = "./Raw_Blasts/"+gene_name + "_raw.fasta"
            gene_raw_fasta_names.append(gene_raw_fasta)
            gene_raw_blast_conv_names.append(gene_raw_fasta+"_raw_blast")
        track_raw.append(gene_raw_blast_conv_names)
            
        #determine which raw_blasts still need to be run:
        run_blast = []
        query_runblast = []
        convert_blast=[]
        timer = 0
        #        print(gene_raw_fasta_names)
        for f in gene_raw_fasta_names:
            #check if the blast has already been run
            g = os.path.isfile(f+"_blast.txt")
            #if not, add it to the list of blasts to run
            if g == False:
               run_blast.append(f)
               query_runblast.append(gene_fasta_list[timer])
            #if blast exists, check that conversion has happened
            else:
                j = os.path.isfile(f+"_raw_blast")
                #if not, add to list of needs to be converted
                if j == False:
                    #add needs to be converted, not run_blast
                    convert_blast.append(f)
            timer+=1
                #we need to send which query to run blast as well as which filename? else blast_)query will fail
                #timer will set this correctly (probably)
        h = []
        timer2 = 0

        #run blast is a list of 30 29's right now. why.
        if len(run_blast) > 0:
            print("Need to run blast on "+str(len(run_blast))+" items")
            print(run_blast)
        time1 = time.clock()
        everything = ""
        #we will run 4 at a time for efficiency's sake, tracked by variable iteration
        iteration = 0
        for item in run_blast:
            QUERY = query_runblast[timer2]
            timer2+=1
            h.append(Fasta(item+"_raw_blast"))
            this_blast_query = "blastp -query " + QUERY + " -remote -db nr -out " + item+"_blast.txt -max_target_seqs 30000 -evalue 1e-4 -outfmt \"6 sallseqid salltitles sseq\""
            this_blast_query = this_blast_query+" & sleep 5 ; "
            iteration += 1
            #if we have 4 loaded: go
            if iteration > 3:
            #reset iteration
                iteration = 0
                everything = everything+this_blast_query
                everything = everything+"wait"
                print("Beginning 4 blast searches")
                print(everything)
                c = os.system(everything)
                everything = ""
            else:
                everything = everything+this_blast_query
        #catch the final blast(s)
        everything = everything+" wait"
        if len(run_blast) > 0 :
            print("Beginning final blast searches")
            print(everything)
            c = os.system(everything)
            print("Done with all blast searches")
            time2 = time.clock()
            print("That took time = ")
            print(time2-time1)
        for newthing in convert_blast:
            h.append(Fasta(newthing+"_raw_blast"))
        if len(convert_blast)>0:
            print("Need to do conversion on "+str(len(h))+" items")
        #print(h)
        for i in h:            
            #gen the lists including species list from the raw blast data
            #print("this")
            i.gen_raw_blast_lists(i.name[:-10]+"_blast.txt")
            #write the parsed blast to fasta. this fasta should have each species_name listed once. 
            i.gen_new_fasta(i.name)

            #finally, new in APRIL, replace that version will a version that has the FULL SEQUENCE for each protein ; no partials
            #get a full version of the current shortened blast.
            #TEST THIS UNIT
        #    print("EXPERIMENTAL FULL PROTEIN CONVERSION ASPECT RUNNING")
            
        #    os.system("full -d 2 -f "+i.name)
            #should have produced a file called PREFIX_raw_Full.fasta
        #    os.system("mv "+i.name+" "+i.name+".aligned_version")
            #turns PREFIX_raw.fasta_raw_blast into PREFIX_raw.fasta_raw_blast.aligned_version
        #    os.system("mv "+i.name[:-16]+"_Full.fasta "+i.name)
            #turns PREFIX_raw_Full.fasta into PREFIX_raw.fasta_raw_blast so it matches with previously set values.
            #i know it's a stupid name blame past abby
            
            
        # for e in list_raw_gene_fastas:  
        #     e.gen_original_lists(e.name+"_raw_blast")
        #     e.gen_species_blast()
        #DONE
        #correlate the generated blast_raw and blast_used files to the names??
        list_of_gene_name_lists.append(gene_name_list)
        print("Blast files for: "+gene_sequences_file+" exist")
    genefile_stored_raw_tuples = []
    for i in range(len(list_of_genelists)):
        genefile_stored_raw_tuples.append( (list_of_genelists[i],track_stored[i],track_raw[i], list_of_gene_name_lists[i]) )
    tuple_blast_objects = Create_Blast_Objects(clade_align, genefile_stored_raw_tuples)
    return tuple_blast_objects

#8
def Get_Attach_Boots(list_of_big_clades, clade_to_align, projectname):
    print("Attaching boots!")
    return 
    boots_list = []
    for item in list_of_big_clades:
        if item in clade_to_align:
            boot_name = "RAxML_bootstrap."+item.prefix+"_CC"
            item.boot_species_file = boot_name
            os.system("scp "+clus_head+"Species_Trees/"+projectname+"/"+boot_name+" ~/Documents/MakeSpeciesTrees/"+projectname+"/Species_Trees/boots/"+boot_name)
            #scp the boot home
            boots_list.append(boot_name)
            #add the boot name to a boots list
        else:
            #look for a boot_file already in place
            boot_name = "RAxML_bootstrap."+item.prefix+"_CC"
            ex_b = os.path.isfile("./"+projectname+"/Species_Trees/boots/"+boot_name)
            if ex_b is True:
                item.boot_species_file = boot_name
                boots_list.append(boot_name)
            else:
                #check online:
                 boot_name = "RAxML_bootstrap."+item.prefix+"_CC"
                 item.boot_species_file = boot_name
                 os.system("scp "+clus_head+"Species_Trees/"+projectname+"/"+boot_name+" ~/Documents/MakeSpeciesTrees/"+projectname+"/Species_Trees/boots/"+boot_name)
                 #scp the boot home
                 boots_list.append(boot_name)
    print("should have moved home:")
    print(boots_list)
    return boots_list


#7
def Correlate_Tree_Subtree(rax_sp_out, list_of_big_clades, cc_files):
#correlate the species trees with our subtree objects
    raxset = 0
    total = len(list_of_big_clades)
    totrax = len(rax_sp_out)
    cc_set = 0
    for tree in cc_files:
        #print(rax_sp_out)
        try:
            pref,junk = tree.split("_CC")
        except:
            pref= tree[:-3]
        if "." in pref:
            rax, pref = pref.split(".")
        if "/" in pref:
            p_l = pref.split("/")
            pref = p_l[-1]
        #v = len(projectname)
        #pref = pref[v:]
        #print(pref)
        for st_o in list_of_big_clades:
            #print(pref)
            #print(st_o.ret_prefix())
          #  print(st_o.ret_prefix()+"  :   "+pref)
            if st_o.ret_prefix() == pref:
                st_o.species_fasta = tree
                sp_f_ob = Fasta(tree)
                sp_f_ob.gen_original_lists(tree)
                sp_f_ob.gen_species_lists()
                st_o.species_fasta_object = sp_f_ob
                #print(st_o.set_species_fasta_object)
                cc_set += 1
    print("there were "+str(cc_set)+" cc_fastas matched with clade objects")
    for tree in rax_sp_out:
        #print (tree)
        #print(rax_sp_out)
        pref,junk = tree.split("_CC")
        if "." in pref:
            rax, pref = pref.split(".")
        #v = len(projectname)
        #pref = pref[v:]
        if "/" in pref:
            p_l = pref.split("/")
            pref = p_l[-1]
        #print(pref)
        for st_o in list_of_big_clades:
            #print(pref)
            #print(st_o.ret_prefix())
           # print st_o.prefix
            if st_o.prefix == pref:
                st_o.set_species_tree_name(tree)
                st_o.species_besttree_name = "RAxML_bestTree."+st_o.prefix+"_CC"
                raxset +=1
    if raxset == total:
        print("All species trees are complete!")
        incomplete = []
    else:
        print("of "+str(total)+" clade objects, there were: "+str(totrax)+" raxml trees generated and :"+str(raxset)+" were successfully matched with the correct subtree object!")
        incomplete = []
        for item in list_of_big_clades:
            if item.species_tree_name == "":
                incomplete.append(item)
    return incomplete

#6
def compare_the_species_tree_taxa_with_the_original_taxa(cc_sp_list, old_sp_list, projectname):
    print("looking at the species tree vs the original request")
    print("#################################################")
    new_sp_list = cc_sp_list
    list_in_both = []
    list_in_old = []

    for item in old_sp_list:
        if item in new_sp_list:
            list_in_both.append(item)
        else:
            list_in_old.append(item)

    print("taxa not found, thus not in species tree:")
    for taxa in list_in_old:
        print(taxa)
    with open("species_tree_generation_info_"+projectname+".txt", "w") as info:
        info.write("species NOT FOUND:")
        for item in list_in_old:
            info.write(item)
        info.write("\n\nspecies that were found and included:")
        for item in list_in_both:
            info.write(item)


#5 main
def Run_Everything_On_SS_Tree(clade_to_align, to_add_to_cc_files, to_add_to_rax_lists, projectname, fixmissing, strain):
    #need to edit this to do only 10 bootstraps
    #
    #currently a wrapper for all the running things functions.
    #merge w/ get multiple alignments?
    if clade_to_align == []:
        print("nothing to align, moving on")
        rax_sp_out = []
        cc_files = []
    else:
        #make sure all blasts are done, and present, for all gene_files
        #species_lists_original need to  exist by here.
        #make them when we initialize the fasta object?
        print("Checking if blasts for all query genes have been run")
        tuple_blast_objects = Check_Blasts_Exist(clade_to_align)
        print("Assigning gene-specific blast files to clades")
        Assign_Blasts_To_Clades(clade_to_align, tuple_blast_objects)
        print("Initializing fastas containing given genes from relevant species")
        #here we need to get clade.species_list_plus_og_loss to exist
        #do we even have species_list existing yet?????
        cc_files = Get_Sequences(clade_to_align, projectname, fixmissing, strain)     
        print("Muscle alignment has finished... running raxml on the resulting files)")
        #raxspone is bipartitions TREES.
        #problem: this is running it in SOD7 not SOD7_SS.
        #this needs to be sent "10" instead of "100"
        rax_sp_one  = raxml_run_on_cluster(cc_files, projectname, False, 10)
        rax_sp_out = rax_sp_one
    for rem in to_add_to_cc_files:
        cc_files.append(rem)
    for rem in to_add_to_rax_lists:
        rax_sp_out.append(rem)        
        #rax_sp_out = ReRoot_Trees(rax_sp_one, list_of_roots)
    try:
        for item in results_list:
            pass
    except:
        results_list = []
    print(to_add_to_rax_lists)
    return cc_files, rax_sp_out, results_list


#4
def Check_For_Existance(list_of_big_clades, projectname):
    #this prevents wasting resources re-making a tree that already exists.
    clade_to_align = []
    existing_rax_clades = []
    to_remove_lists = []
    to_remove_lists_err = []
    to_add_to_rax_lists = []
    to_add_to_cc_files = []
    for item in list_of_big_clades:
        prefix = item.ret_prefix()
        outputfile_exists = os.path.isfile(prefix+"_CC_Rax_Bipart.newick")
        #maybe is the renaming not happening anymore? check this in terminal.
        if outputfile_exists is True:
            print("found: "+prefix+"bipart file")
            bestfile = os.path.isfile("RAxML_bestTree."+prefix+"_CC")
            if bestfile is True:
                existing_rax_clades.append(item)
                to_add_to_rax_lists.append(prefix+"_CC_Rax_Bipart.newick")
                to_add_to_cc_files.append(prefix+"_CC.fasta")
            else:
                best_stored_exists = os.path.isfile("./"+projectname+"/Species_Trees/trees/RAxML_bestTree."+prefix+"_CC")
                if best_stored_exists is True:
                     existing_rax_clades.append(item)
                     to_add_to_rax_lists.append(prefix+"_CC_Rax_Bipart.newick")
                     to_add_to_cc_files.append(prefix+"_CC.fasta")
                else:
                    clade_to_align.append(item)

        else:
            #check the proper sub_folder for it.
            output_stored_exists = os.path.isfile("./"+projectname+"/Species_Trees/trees/"+prefix+"_CC_Rax_Bipart.newick")

            if output_stored_exists is True:
                print("got bipart "+prefix)
                best_stored_exists = os.path.isfile("./"+projectname+"/Species_Trees/trees/RAxML_bestTree."+prefix+"_CC")
                if best_stored_exists is True:
                    existing_rax_clades.append(item)
                    to_add_to_rax_lists.append("./"+projectname+"/Species_Trees/trees/"+prefix+"_CC_Rax_Bipart.newick")
                    to_add_to_cc_files.append("./"+projectname+"/Species_Trees/concat/"+prefix+"_CC.fasta")
                else:
                    print("missing besttree? final check:")
                    here_bestfile = os.path.isfile("RAxML_bestTree." + prefix + "_CC")
                    if here_bestfile is True:
                        print("gotcha!")
                        existing_rax_clades.append(item)
                        os.system("cp RAxML_bestTree." + prefix + "_CC ./"+projectname+"/Species_Trees/trees/")
                        to_add_to_rax_lists.append("./" + projectname + "/Species_Trees/trees/" + prefix + "_CC_Rax_Bipart.newick")
                        to_add_to_cc_files.append("./" + projectname + "/Species_Trees/concat/" + prefix + "_CC.fasta")
                    else:
                        clade_to_align.append(item)

            else:
                print("didn't find: ./"+projectname+"/Species_Trees/trees/"+prefix+"_CC_Rax_Bipart.newick")
                print(os.getcwd)
                clade_to_align.append(item)
                #check if we already make the concat file, but didn't finish the job.
                # outputfile_CC_only_exists = os.path.isfile(item[3]+"_CC.fasta")
                # if outputfile_CC_only_exists is True:
                #     print("found CC_only of: "+item[3]+" assuming removed for a reason.")
                #     to_remove_lists_err.append(item)
                #     to_add_to_cc_files.append(prefix+"_CC.fasta")
                # else:
    #this is only checking if species trees are done, not the gene trees.
    print("Still_to_align:")
    print(clade_to_align)
    return to_add_to_rax_lists, to_add_to_cc_files, clade_to_align, existing_rax_clades
    #done

#3
def check_directory_existance(prefix, ssh_inst):
    import os
    #print(type(ssh_inst))
    #print(type(prefix))
    #print(prefix)
    os.system(ssh_inst+" \' mkdir Species_Trees;cd Species_Trees;mkdir "+prefix+"\'")

#2 main
def Species_Tree_For_The_Subsampled_Fasta(ss_subtree_object, projectname, cat_file, individual_party = False, fixmissing = False, strain = False):
    #attempting projectname fix here
    #ss_subtree_object.projectname = projectname+"_SS"
    #should already be set (?)
   # print("should be making")
    check_directory_existance(projectname, ssh_inst)
   # raise  SystemExit
    ss_subtree_object.fasta_object.gen_original_lists(ss_subtree_object.fasta)
    sp_list = ss_subtree_object.fasta_object.gen_species_lists()
    ss_subtree_object.species_list_original = sp_list
    #this is the same, since we are not adding a rooting species or any loss candidates this time around.
    ss_subtree_object.species_list_plus_og_loss = sp_list
    ss_subtree_object.cat_file = cat_file

    #for SOD old version only
    #Fix_Your_Mistakes([ss_subtree_object])
    
    #currently a wrapper for all the running things functions.
    print("Checking if files already exist")
    to_add_to_rax_lists, to_add_to_cc_files, subtree_object_in_list, existing_rax_clades = Check_For_Existance([ss_subtree_object], projectname)

    #potential error if species trees completed but gene trees didnt!!!
    print("Beginning tree-making")


    #change this to not do any taxonomy / rooting / loss addition BS
    cc_files, rax_out, results_list = Run_Everything_On_SS_Tree(subtree_object_in_list, to_add_to_cc_files, to_add_to_rax_lists, projectname, fixmissing, strain)
    
    print("Matching trees to corresponding clades")
    incomplete = Correlate_Tree_Subtree(rax_out, [ss_subtree_object], cc_files)

    for item in clade_to_align:
        if item in incomplete:
            clade_to_align.remove(item)

    boots_list = Get_Attach_Boots([ss_subtree_object], clade_to_align, projectname)
    print("Done with everything")
    return rax_out, cc_files, incomplete, boots_list



#1 main
def set_up_an_independent_run(fasta, projectname, directory_to_run_in, fixmissing, strain):
    os.chdir(directory_to_run_in)

    #set up my subtree and fasta objects
    MySubtree = Subtree(projectname)
    MySubtree.fasta = fasta
    MySubtree.fasta_object = Fasta(fasta)
    MySubtree.fasta_object.gen_original_lists(fasta)

    MySubtree.fasta_object.prefix = projectname
    MySubtree.prefix = projectname
    MySubtree.projectname = projectname
    
    #make sure all directories exist 
    os.system("mkdir "+projectname+"; cd "+projectname+"; mkdir Species_Trees; mkdir Gene_Trees; cd Species_Trees; mkdir trees; mkdir construction; mkdir boots; mkdir concat; cd -; cd Gene_Trees; mkdir fasta; mkdir muscle; mkdir boots;  mkdir trees; cd ..; cd ..")

    # for now, use all 30 hopefully homolopgous queries to make the SubSampled Species Tree
    cat_file = New_Gene_List
    #should be defined above
    MySubtree.cat_file = cat_file

    #everything runs here.
    rax_out, cc_files, incomplete, boots_species  = Species_Tree_For_The_Subsampled_Fasta(MySubtree, projectname, cat_file, True, fixmissing, strain)
    print("species trees were probably made and should exist")
    print(rax_out)
    print("DONE HAHAHAHAHAH MAYBE THIS WORKED")

    #directions


#to use:

#1. make a new folder somewhere, and call it something easy, like "SpeciesTreeGen"
#2. put into that folder a. your fasta file with gene sequences
#                        b. the query file (New_Ribo.fasta)
#                        c. optionally, pre-downloaded blasts in a folder called "Raw_Blasts"
#3. replace the strings at the top of this document with your personal credentials eg path to this folder, cluster sign on materials, etc
#4. try and do the thing!



if __name__ == "__main__":
    print("Running in terminal")
    parser = argparse.ArgumentParser(description="All")

    #optional directory
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in eg MakeSpeciesTree")
    parser.add_argument("-p", "--projectname", action = "store", help="type projectname eg SOD8")
    parser.add_argument("-f", "--fasta", action = "store", help="type fasta eg SOD.fasta")
    parser.add_argument("-b", "--blastmissing", action = "store_true", default = False, help="type -b to preform a blast to search for any sequence data that is missing, in case it is available on NCBI but was missing or badly formatted in the original search.")
    parser.add_argument("-s", "--strain", action = "store_true", default = False, help="type -s to ignore strains -- that is, only consider genus_species when building the tree. no subspecies or strains")
   
    args = parser.parse_args()


    print("single tree runthrough only option")
    if args.projectname is False:
        print("...specificity?")
        raise SystemExit
    set_up_an_independent_run(args.fasta, args.projectname, args.directory, args.blastmissing, args.strain)


#strain... set as argument, passed to 1, 2, 5, 12(GET SEQUENCES)


