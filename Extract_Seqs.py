#!/usr/bin/python

#abigailc@Actaeon June 22 2017

#this script assumes you have a fasta file containing a LOT of sequences, and you want to extract a subset of them
#also, you have either :
#    a NEXUS file (eg copy pasted clade from figtree)
            ##NEXUS
            #begin trees;
            #tree tree_1 = [&R] ((Methanosarcina_mazei:0.030681,(Methanosarcina_acetivorans_C2A:0.010784,Methanosarcina_siciliae_C2J:0.017418):0.022073):0.017515,Methanosarcina_lacustris_Z_7289:0.049195);
            #end;
#    a TEXT file (line by line including seqids or keywords) 
            #Methanosarcina_mazei
            #Methanosarcina_acetivorans_C2A
            #Methanosarcina_siciliae_C2J
            #Methanosarcina_lacustris_Z_7289
#    or will MANUALLY enter a keyword
            #Methanosarcina
#    ....for the purposes of constraining the extraction

#finally, specify if you want EXACT MATCHING or SUBSTRING matching.
            # given "Cats" and "Cats2", substring query "Cats" matches both. exact matching "Cats" only matches the first.


import os
class Fasta:
    def __init__(self, name):
        #all ids and seqs should be stripped of leading and trailing whitespace and have ">" removed for reasons.
        #this is the name of the fasta, it can be anything, i'm not really using it right now.
        self.name = name
        #this is the to-be-modified version of sequence IDs and sequence-Data
        # ALWAYS keep IDS and SEQS the same length. id[1] should ALWAYS correspond to seq[1].
        self.ids = []
        self.seqs = []
        # these are the original SEQids and Sequences. They should never be modified after generation in gen_original_lists or blast_to_fasta
        self.original_ids = []
        self.original_seqs = []
        # initialize yourself
        self.gen_original_lists(name)

    #basically, opens the file. loads IDS and SEQS
    def gen_original_lists(self, fastaname):
        print("original lists for : "+fastaname)
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

    #writes a new fasta file
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

    #extract_exact requires perfect matching between list given and fasta_files_sequence_IDS
    #eg "bacteria" will extract >bacteria but not >bacteria1 or >Cyanobacteria
    def extract_exact(self, list_of_keeps):
        keep_ids = []
        keep_seq = []
        success = 0
        suc_num = len(list_of_keeps)
        for item in list_of_keeps:
            item = item.strip()
            for iteration in range(len(self.original_ids)):
                thing = self.original_ids[iteration]
                if thing == item:
                    keep_ids.append(thing)
                    seq = self.original_seqs[iteration]
                    keep_seq.append(seq)
                    success += 1
        if suc_num == success:
            print("100% complete extract")
        else:
            print(str(success)+"out of "+str(suc_num)+" sequences extracted")
        self.ids = keep_ids
        self.seqs = keep_seq

    #extract within requires something from listgiven be a substring of fasta_files_sequence_ID
    #eg: "bacteria" will extract both >Cyanobacteria and >Actinobacteria
    #also useful given gi numbers or other unique identifiers that may or may not be the full seqid
    def extract_within(self, list_of_keeps):
        keep_ids = []
        keep_seq = []
        success = 0
        suc_num = len(list_of_keeps)
        for item in list_of_keeps:
            item = item.strip()
            for iteration in range(len(self.original_ids)):
                seqid = self.original_ids[iteration]
                #eg if "cats" in "this_sequence_from_cats"
                if item in seqid:
                    keep_ids.append(seqid)
                    seq = self.original_seqs[iteration]
                    keep_seq.append(seq)
                    success += 1
        else:
            print(str(success)+" sequences extracted")
        self.ids = keep_ids
        self.seqs = keep_seq
    



def parse_nexus(nexus):
    tips_list = []
    tree = ""
    with open (nexus) as opennexus:
        for line in opennexus:
            if "\ttree " in line:
                tree = line
    if tree == "":
        print("error parsing the nexus file: "+nexus)
        raise SystemExit
    begin = "no"
    for char in tree:
        #start reading seqid if char is ( or ,(and not currently reading)
        if char == "(":
            new = ""
            begin = "yes"
        elif char == ",":
            if begin == "no":
                new = ""
                begin = "yes"
            elif begin == "yes":
                begin = "no"
                tips_list.append(new)
        #seqid finished if : or ,(and currently reading) -- so add it to the list

        elif char == ":":
            if begin == "no":
                pass
            else:
                begin = "no"
            tips_list.append(new)
        elif begin == "yes":
            new = new+char
    if tips_list == []:
        print("error building list of sequence id tips parsing nexus: "+nexus)
        raise SystemExit
    if tips_list[0][0] == "'":
        print("probably all the tips are coated in quotes from processing... removing them")
        fixed = []
        for item in tips_list:
            fitem = item.replace("'", "")
            fixed.append(fitem)
        tips_list = fixed
    #will be a list of seqids (minus the >)
    print(tips_list)
    return tips_list

def parse_text(textfile):
    tips_list = []
    with open (textfile) as text:
        for line in text:
            line = line.strip()
            tips_list.append(line)
    return tips_list

   
if __name__ == "__main__":

    print("Running in terminal")
    import argparse
    parser = argparse.ArgumentParser(description="All")
    #necessary bits
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in where fasta resides, if not pwd")
    parser.add_argument("-fas", "--fasta", action = "store", default = False, help="type the name of your .fasta file (")
    #what to extract
    parser.add_argument("-nex", "--nexus", action = "store", default = False, help="type name of nexus constraining file")
    parser.add_argument("-txt", "--txtfile", action = "store", default = False, help="type name of text file containing constraining seqids or strings on seperate lines")
    parser.add_argument("-man", "--manual_input", action="store", default= False, help="type string for constraining")
    #flags
    parser.add_argument("-ex", "--exact", action = "store_true", default = False, help="toggle for exact seqid matching")
    parser.add_argument("-sub", "--substring", action="store_true", default= False, help="toggle for substring matching")
    #output
    parser.add_argument("-wf", "--write_fasta", action="store", default= False, help="give a name for the output file")

    args = parser.parse_args()


    #change dir if desired
    try:
        os.chdir(args.directory)
        if args.verbose == True:
            print("moved to dir: "+args.directory)
    except:
        print ("didn't change dir")

    #load the base fasta (that contains all the sequences)
    myfasta = Fasta(args.fasta)

    #generate the list of constraining things
    if args.nexus != False:
        constraint = parse_nexus(args.nexus)
    if args.txtfile != False:
        constraint = parse_text(args.txtfile)
    if args.manual_input != False:
        constraint = args.manual_input.split(" ")

    #do the extraction, either exact or substring version
    if args.substring is True:
        myfasta.extract_within(constraint)
    if args.exact is True:
        myfasta.extract_exact(constraint)



    #write the result to file
    if args.write_fasta is False:
        print("silly you, I can't write an output without a name for it! try using -wf OUTPUT")
        raise SystemExit
    else:
        myfasta.gen_new_fasta(args.write_fasta)
    print("All done! hopefully! your output file is at: "+args.write_fasta)












