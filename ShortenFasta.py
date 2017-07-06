#!/usr/bin/python
#last edit abigailc@Actaeon Oct 18 2016
#this is a very short script - a piece of FISH that just shortens the output from blast on NCBI's online interface
#it will also work on local blasts run through blast+ if your outfmt is "6 sallseqid salltitles sseq" or "6 sseqid stitle sseq"


#here is a class called Fasta. it includes:
#  a function to add sequence data and sequenceID data (open_seqs) as self.original_(ids/seqs)
#  a function to modify the self.original_ids, and save the output at the same index in self.new_ids
#  a function to write that modified data to .fasta

class Fasta:
    def __init__(self, name):
        #name can be anything; I usually store the name of the file here.
        #to get the name returned, type FastaName.name (if out of class def) or self.name (if within definition)
        self.name = name
        #initialize lists
        # ALWAYS keep IDS and SEQS the same length. id[1] should ALWAYS correspond to seq[1].
        # this is implemented as lists due to ordering properties and the ease of adding more list (eg, one that keeps track of just gi nums or species names can be easily added later. dictionaries make that more difficult, though harder to misorder error.
        #these will be the modified versions of seq and seqID
        self.new_ids = []
        self.new_seqs = []
        # these are the original SEQids and Sequences. They should never be modified after generation in open_seqs or blast_to_fasta
        self.original_ids = []
        self.original_seqs = []
    def open_seqs(self, fastaname):
        #open the fasta name. if you saved it as self.name, can call FastaName.open_seqs(FastaName.name) or FastaName.open_seqs("example.fasta")
        with open(fastaname) as fastafile:
            #parse per line
            for line in fastafile:
                #avoid badly written start or end lines
                if "\n" == line:
                    pass
                #this marks the beginning of a seqid
                if ">" in line:
                    #write the previous AA seq (if one exists)
                    try:
                        #this will fail on the first iteration, since AAseq has not yet been defined.
                        AAseq=AAseq.strip()
                        self.new_seqs.append(AAseq)
                        self.original_seqs.append(AAseq)
                    except:
                        pass
                    #initialize a new AAseq
                    AAseq = ""
                    #format the seqID to have no whitespace or linebreaks in front/back, and remove the ">" during storage
                    newline = line.strip()
                    newline = line.strip(">")
                    #save the seqID
                    self.new_ids.append(newline.strip())
                    #populate new_ids even though these are identicle to the originals right now, so that it will be the same length, and if any are NOT modified, they will still print.
                    self.original_ids.append(newline.strip())
                else:
                    #this line is not a new seqID, so it is part of an AA sequence, store it as a continutation of AAseq (and remove whitespace/linebreaks)
                    AAseq = AAseq+line
                    AAseq= AAseq.strip()
            #usually AAseq will be written on the next instance of ">", but need to also catch the last pass
            self.new_seqs.append(AAseq)
            self.original_seqs.append(AAseq)
        #tell me if it worked, and ret the number of seqs we opened just fyi
        print("Initial sequence and ID lists created. Contains "+str(len(self.original_ids))+" sequences")
    #EDITED to work with new NCBI download format
    def shorten(self):
        print("shortening ids...")
        #set defaults
        unk = "no"
        normal = 0
        ucount = 0
        #go line by line through ids
        for line in self.original_ids:
            #get index so we know which new_id to modify
            index = self.original_ids.index(line)
            # changes NCBI's default naming scheme to be
            #>Species_name|#########
            #where ####is either gi number or accession number, depending on new or old input
            # AAH91460.1 Ribosomal protein L3 [Danio rerio] (new format)
            #this works w old format
            if "gi|" in line:
                #get gi num
                number = re.sub("(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\3", line)
                num = number.strip()
                #get species name
                edit1 = re.sub("(gi)(\|)([0-9]*)(\|)([A-Za-z]*)(\|)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\8|", line)
            #get acc number
            else:
                #get accession num
                number = re.sub("([^ ]*)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\1", line)
                num = number.strip()
                #get species name
                edit1 = re.sub("([^ ]*)(.*)(\[\'?[A-Z]?[a-z]* ?.*\])(.*)", "\\3|", line)
            #check if species name recovered, if not toggle unknown, if so add to normal counter
            if "[" in edit1:
                unk = "no"
                normal += 1
            else:
                unk = "yes"
            #substitution to remove any weird shit in species name.
            edit2 = re.sub("[\[\]]", "", edit1)
            #for now, leave periods in number (if accession num) but not name
            edit3 = re.sub("[:;\.=,/\+'\(\)]", "_", edit2)
            edit4 = re.sub(" ", "_", edit3)
            edit4 = re.sub("__", "_", edit4)
            edit4 = edit4+num
            #if its good, change the entry in new_ids at given index
            if unk == "no":
                self.new_ids[index] = edit4
            else:
                print("Unknown Species in ID:" + line)
        print("shortened: "+str(normal)+" sequence")
    def gen_new_fasta(self, new_fasta_name):
        #this should print the changed seqids and changed AA sequences to file.
        newfasta = new_fasta_name
        #open the new file
        with open (newfasta, "w") as new:
            #for each sequence...
            for i in range(len(self.original_ids)):
                #write the id + linebreak
                new.write(">"+self.new_ids[i].strip()+"\n")
                #write the sequence + linebreak
                new.write(self.new_seqs[i]+"\n")
        print("Finished, your new fasta file is located at "+newfasta)
        #done
## adding in support for ten characters conversion and unconversion.
    def ten_char(self):
        #something
        #this one should be done in a seperate loop
        CTdict = {}
        iteration = 0
        for line in self.ids:
            iteration +=1
            line = line.strip()

##                  #i have something like
##                  >Methanococcoides_burtonii|gi|909890
##                  #i want
##                  >MethBurt00
            GenusSpecies = re.sub("([A-Z][a-z]*)(_)([A-Z]*[a-z]*)(.*)", "\\1~\\3", line)
            try:
                Genus, Species = GenusSpecies.split("~")
                g4 = Genus[:4]
                try:
                    s4 = Species[:4]
                    s3 = Species[:3]
                except:
                    s4 = Species[:2]
                    s3 = Species[:2]
                if iteration < 10:
                    newid = g4+s4.capitalize()+"0"+str(iteration)
                elif iteration > 99:
                    newid = g4+s3.capitalize()+str(iteration)
                else:
                    newid = g4+s4.capitalize()+str(iteration)
            except:
##                      print(GenusSpecies)
                gs8 = GenusSpecies[1:9]
                if iteration < 10:
                    newid = gs8+"0"+str(iteration)
                elif iteration > 99:
                    newid = gs8[:-1]+str(iteration)
                else:
                    newid = gs8+str(iteration)
##                      print(newid)

            CTdict[line] = newid
        for line in self.ids:
            index = self.ids(line)
            newestid = CTdict[line]
            self.ids[index] = newestid
        print("ten char done")

    # def load_info_swap(self, info_file_in):
    #     #reads a file of form
    #     #   originalID
    #     #   changedID
    #     #and generates self.ids from that file.
    #     kid = "no"
    #     vid = "no"
    #     CTdict = {}
    #     with open (info_file_in) as old:
    #         for line in old:
    #             #first pass: gets key (original ID)
    #             #second pass: gets value (new ID)
    #             #if we have no info, get key
    #             if kid == "no":
    #                 key = line.strip()
    #                 kid = "yes"
    #                 continue
    #             elif kid == "yes":
    #                 #if we have key and value, record.
    #                 if vid == "yes":
    #                     CTdict[key]=value
    #                     vid = "no"
    #                     kid = "no"
    #                     continue
    #                 #if we have key but no value, get value.
    #                 if vid == "no":
    #                     value = line.strip()
    #                     vid = "yes"
    #         #catch the final pass
    #         CTdict[key]=value

    #     if self.original_ids == []:
    #         for thing in CTdict:
    #             self.ids.append(thing)
    #             self.original_ids.append(CTdict[thing])
    #     else:
    #         for item in self.original_ids:
    #             index = self.original_ids.index(item)
    #             newid = CTdict[item]
    #             self.ids[index] = newid
    #     print("original ids:")
    #     print(self.original_ids)
    #     print("new ids:")
    #     print(self.ids)

    # def gen_info(self, info_file_name):
    #     #writes a file of form
    #     #   originalID
    #     #   changedID
    #     with open(info_file_name, "w") as inf:
    #         listlength = len(self.original_ids)
    #         if listlength != len(self.ids):
    #             print ("List lengths do not match! FATAL ERROR")
    #             print (self.original_ids)
    #             print (self.ids)
    #             raiseSystemExit
    #         for i in range(listlength):
    #             inf.write(self.original_ids[i])
    #             inf.write(self.ids[i]+"\n")
    #     print("Info file was generated. Named "+info_file_name)
    #     #done

    # def swap_in_newick(self, old_newick_name, new_file_name):
    #     #this replaces the tip names in a newick file. sometimes works on nexus files too, but I havent extensively tested it.
    #     newick = old_newick_name
    #     newnewick = new_file_name
    #     with open (newick) as old:
    #         with open (newnewick, "w") as new:
    #             for line in old:
    #                 for item in self.original_ids:
    #                     index = self.original_ids.index(item)
    #                     line = line.replace(item, self.ids[index])
    #                 new.write(line)
    #     print("finished, tip-replaced-newick file at: "+newnewick)
    #     #done



        
#here is the parser... this bit will only run if called from command line
if __name__ == "__main__":

    print("Running in terminal")
    #imports (not sure we use them all)
    import sys
    import argparse
    import os
    import re
    parser = argparse.ArgumentParser(description="All")

    #optional directory
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in where fasta resides, if not pwd")
    parser.add_argument("Fasta", type=str, help="type name your .fasta")
    parser.add_argument("-t", "--tenchars", action = "store_true", help="flag -t means makes a new fasta file where names are 10chars long")

    args = parser.parse_args()

    
        
        
    ##begin commands
    
    #change dir if desired
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")


    #create the object
    ExampleFasta = Fasta(args.Fasta)
    
    #populate the lists
    ExampleFasta.open_seqs(args.Fasta)
    
    #shorten the ids
    ExampleFasta.shorten()

    #gen an output name
    try:
        start, ext = args.Fasta.split(".")
        outname = start+"Sh.fasta"
    except:
        outname = args.Fasta+"Sh.fasta"
    #write the new fasta
    ExampleFasta.gen_new_fasta(outname)


    if args.tenchars is True:
        ExampleFasta = Fasta(outname)
    #populate the lists
        ExampleFasta.open_seqs(outname)
    
    #shorten the ids
        ExampleFasta.shorten(tenchars)
        try:
            start, ext = outname.split(".")
            outname2 = start+"Tenchars.fasta"
        except:
            outname2 = outname+"Tenchars.fasta"
    #write the new fasta
        ExampleFasta.gen_new_fasta(outname2)
    print("Operation finished, closing!")
