#!/usr/bin/python
#abigailc@Actaeon Jan 4 2017
#classes used in OVerall_DTL_Detector.py
import sys
import argparse
import os
import re
import time
class Fasta:
        def __init__(self, name):
                #all ids should be stripped and have ">" removed for reasons.
                #for now, sequences do not have any stripping applied
                self.name = name
                self.ids = []
                self.original_ids = []
                self.original_seqs = []
                self.seqs = []
                self.species_names = []
        def gen_original_lists(self, fastaname):
                if self.original_ids != []:
                        self.ids = []
                        self.original_ids = []
                        self.original_seqs = []
                        self.seqs = []
                        self.species_names = []
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
                                # print(i)              #
                                #unclear if this needs a "\n" after it... check.#TODO
                                #print(self.seqs)
                                #print(type(self.seqs[i]))
                                new.write(self.seqs[i]+"\n")
                print("Finished, your new fasta file is located at "+newfasta)
                #done
        def extract(self, list_of_keeps):
                keep_ids = []
                keep_seq = []
                success = 0
                suc_num = len(list_of_keeps)
                or_num = len(self.original_ids)
                for item in list_of_keeps:
                        item = item.strip()
                        found = "n"
                        for thing in self.original_ids:
                                if thing.strip() == item:
                                        keep_ids.append(thing)
                                        index = self.original_ids.index(item)
                                        seq = self.original_seqs[index]
                                        keep_seq.append(seq)
                                        success += 1
                                        #print("matched:"+item+":with:"+thing.strip())
                                        found ="y"
                        if found == "n":
                            print ("could not find in .fasta the tip:"+item)
                if suc_num == success:
                        print("100% complete extract")
                else:
                        print(str(success)+"out of "+str(suc_num)+" sequences extracted")
                        #print("looked for")
                        #print(list_of_keeps)
                        #print("in")
                        #print(self.original_ids)
                self.ids = keep_ids
                self.seqs = keep_seq

        def gen_gis_list(self):
            gilist = []
            for item in self.ids:
                #print(item)
                taxon = re.sub("(.*)(gi#\|?)([0-9]*)(.*)", "\\3", item)
                #print(taxon)
                if "|" in taxon:
                    print("TAXON error in gen_gis_lists():" + taxon)
                    bleh, taxon = taxon.split("#")
                #print(taxon)
                gilist.append(taxon)
            self.gis_list = gilist
            return gilist

        def gen_species_lists(self):
            speclist = []
            for item in self.ids:
                taxon = re.sub("([^_]*)([A-Z][a-z]*_?[A-Z]?[a-z]*[^\|]*)(.*)", "\\2", item)
                if "|" in taxon:
                    tlist = item.split("|")
                    taxon = tlist[-2]
                    if "|" in taxon:
                        print ("TAXON error in gen_species_lists():" + taxon)
                speclist.append(taxon)
            self.species_names = speclist
            return speclist

        def ret_speclist(self):
            return self.species_names


class Subtree:
    def __init__(self, name):
    #all ids should be stripped and have ">" removed for reasons.
    #for now, sequences do not have any stripping applied
        # 123
        self.number_name = name
        #Cyanobacteria
        self.string_name = ""
        #Bacterua
        self.category = ""
        #[A,B,C]
        self.tips = []
        # cyano.fasta
        self.fasta = ""
        # object<Fasta(cyano)>
        self.fasta_object = ""
        # object<Fasta(species_cyano)>
        self.species_fasta_object = ""
        # object<Fasta(gene_cyano)>
        self.gene_fasta_object = ""
        #the subtree as generated in figtree
        self.fasttree_st = ""
        #i expect this to be the RAxML_BestTree but it looks like its a fasta instead?
        self.gene_tree = ""
        #RAxML_besttree
        self.gene_tree_name = ""
        self.species_tree_name = ""
        self.gene_tree_species_tips = ""
        self.species_tree = ""
        self.species_list_original = []
        self.species_list_gene_to_species = []
        self.species_list_after_removal = []
        self.species_list_plus_og_loss = []
        self.besttree_str = ""
        self.prefix = ""
        self.cladetype = ""
        self.projectname = ""
    def set_gene_tree_species_tips(self):
        a = self.gene_tree
        if self.fasta_object == "":
            if self.fasta == "":
                print("this isn't going to work")
                raise SystemExit
            self.fasta_object = Fasta(self.fasta)
        b = self.fasta_object
        c = b.gen_species_lists()
        d = b.gen_gis_list()
        print(d[0])
        old = b.ids
        newt = a
        new = []
        length = len(b.species_names)
        print(length)
        for i in range(length):
            new.append(b.species_names[i]+"_"+d[i])
        for item in old:
            index = old.index(item)
            newt = newt.replace(item, new[index])
        self.gene_tree_species_tips = newt
    def set_alignment_with_species_tips(self):
        #make a new .fasta object and then change its tips.
        a = Fasta(self.fasta)
        try:
            a.gen_original_lists(self.fasta)
        except:
            a.gen_original_lists("./"+self.projectname+"/Gene_Trees/muscle/"+self.fasta+"_Muscle.fasta")
        self.fasta_object_with_species_names = a
        
        b = a.gen_species_lists()
        c = a.gen_gis_list()
        old = a.original_ids
        new = []
        for i in range(len(a.species_names)):
            new.append(a.species_names[i]+"_"+c[i])

        for item in old:
            index = old.index(item)
            a.ids[index] = new[index]
        fasta_with_species_names = a.gen_new_fasta("./"+self.projectname+"/Gene_Trees/fasta/"+self.prefix+".gene_sp_names.fasta")
        return fasta_with_species_names
        #edit april ^
    def ret_gene_tree_species_tips(self):
        return self.gene_tree_species_tips
    def set_fasta_object(self):
        MyFasta = Fasta(self.fasta)
        MyFasta.gen_original_lists(self.fasta)
        MyFasta.gen_species_lists()
        self.fasta_object = MyFasta
        
    def set_species_fasta_object(self, input_fasta):
        MyFasta = Fasta(input_fasta)
        MyFasta.gen_original_lists(self.fasta)
        MyFasta.gen_species_lists()
        self.species_fasta_object = MyFasta
        self.species_list_after_removal = MyFasta.gen_species_lists()
    def set_species_list(self, splist):
        self.species_list = splist
    def set_type(self, typ):
        self.cladetype = typ
    def set_prefix(self, pref):
        self.prefix = pref
    def set_species_tree_name(self,speciest):
        self.species_tree_name = speciest
        with open(speciest) as old:
            self.species_tree = old.read().strip()
    def set_gene_tree_name(self, gene_tree):
        self.gene_tree_name = gene_tree
        with open(gene_tree) as old:
            self.gene_tree = old.read().strip()
    def set_fasttree(self, ft):
        self.fasttree_st=ft
    def set_string(self, stringn):
        self.string_name = stringn
    def set_category(self, cat):
        self.category = cat
    def set_fasta(self, fasta):
        self.fasta = fasta
    def set_tips(self, tips):
        if type(tips) == str:
            self.tips.append(str)
        else:
            for item in tips:
                self.tips.append(item)
    def ret_fasta_object(self):
        return self.fasta_object
    def ret_number(self):
        return self.number_name
    def ret_type(self):
        return self.cladetype
    def ret_species_list(self):
        return self.species_list
    def ret_string(self):
        return self.string_name
    def ret_name(self):
        return self.number_name
    def ret_prefix(self):
        return self.prefix
    def ret_cat(self):
        return self.category
    def ret_tips(self):
        return self.tips
    def ret_fasta(self):
        return self.fasta
    def ret_fasttree(self):
        return self.fasttree_st
    def ret_gene_tree(self):
        return self.gene_tree
    def ret_species_tree(self):
        return self.species_tree
