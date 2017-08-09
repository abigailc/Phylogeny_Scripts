#local blast helper


#local blast of gene = gene
#with conversion to readableformat
#and seperation of multispecies hits

#local blast for SPECIFIC species_name:


def Blast_Query_Get_All(query, evalue, max_hits):
    #get output name
    op = query+"_"+str(max_hits)+"_blast.txt"

    #the blast
    os.system("blastp -query "+query+" -remote -db nr -out "+op+" -max_target_seqs "+max_hits+" -evalue "+evalue+" -outfmt \"6 sseqid stitle sseq\"")

    return op

def Blast_Query_Get_Specific(query, sp, evalue, MYFASTA):
    #getoutput_name
    OUTPUT = "tempout"
    #the blast, returns a file -- going to need to merge them --- is merge a script? 
    blast_query = "blastp -remote -query " + query + " -db nr -out " + OUTPUT +" -max_target_seqs 10 -entrez_query " + sp + " -evalue "+evalue+" -outfmt \"6 sallseqid salltitles sseq\""
    os.system(blast_query)
    tempfasta = Fasta(OUTPUT)
    a = tempfasta.gen_raw_blast_lists(OUTPUT)
    #a = tempfasta.blast2fasta(OUTPUT+"b")
    #tempfasta.gen_original_lists(OUTPUT)
    try:
        a = tempfasta.ids[0]
    except:
        print("couldn't make an original list from result")
        return "None"
    #tempfasta.gen_new_fasta("postb2f_b.txt")
    #a = tempfasta.single_shorten()
    gis = tempfasta.gen_numbers()
    a = tempfasta.gen_species_lists()
    found = False
    if " " in sp:
        sp = sp.replace(" ", "_")
    if "\"" in sp:
        sp = sp.replace("\"", "")

    print("looking for: "+sp)
    for i in range(len(tempfasta.original_ids)):
        if sp == tempfasta.species_names[i]:
            found = True
            #add this sequence and id and then exit
            MYFASTA.original_seqs.append(tempfasta.seqs[i])
            MYFASTA.original_ids.append(tempfasta.ids[i])
            MYFASTA.seqs.append(tempfasta.seqs[i])
            MYFASTA.ids.append(tempfasta.ids[i])
    if found == False:
        gis = tempfasta.gen_numbers()
        taxids = tempfasta.SetTaxID()
        rankdict = tempfasta.GetTaxonomy()
        a = tempfasta.single_shorten()
        withtaxonomy = tempfasta.AppendTaxonomy(["ScientificName"])
        a = tempfasta.gen_species_lists()
        if sp == tempfasta.species_names[i]:
            found = True
            #add this sequence and id and then exit
            MYFASTA.original_seqs.append(tempfasta.seqs[i])
            MYFASTA.original_ids.append(tempfasta.ids[i])
            MYFASTA.seqs.append(tempfasta.seqs[i])
            MYFASTA.ids.append(tempfasta.ids[i])
    if found == False:
        print("couldnt match:")
        print(sp)

def Run_Multiple_Specifics(query, species_names, evalue):
    list_all = []
    for item in query:
        output = item+"_specifics.fasta"
        MYFASTA = Fasta(output)
        for sp in species_names:
            a = Blast_Query_Get_Specific(item, sp, evalue, MYFASTA)
        MYFASTA.gen_new_fasta(output)
        list_all.append(output)
    return(list_all)

def Blast_Query_Get_Entrez(query, ez, evalue, max_hits):
    #getoutput name
    op = query+"_"+ez+"_"+str(max_hits)+"_blast.txt"
    #theblast
    os.system("blastp -query "+query+" -remote -db nr -out "+op+" -max_target_seqs "+max_hits+" -evalue "+evalue+" -entrez_query "+ez+" -outfmt \"6 sseqid stitle sseq\"")
    return op

def Run_Multiple_Queries(query_list, species_names, entrez, evalue, max_hits):
    list_of_blast_outputs = []
    print(type(query_list))
    if species_names != False:
        list_of_blast_outputs = Run_Multiple_Specifics(query_list, species_names, evalue)
    elif entrez_list != False:
        for item in query_list:
            a = Run_Multiple_Entrez(item, entrez_list, evalue, max_hits)
            for output in a:
                list_of_blast_outputs.append(output)
                print("A")
                print(list_of_blast_outputs)
    else:
        for item in query_list:
            a = Blast_Query_Get_All(item, evalue, max_hits)
            list_of_blast_outputs.append(a)
    return list_of_blast_outputs

def Run_Multiple_Entrez(query, entrez_names, evalue, max_hits):
    list_specifics = []
    for ez in entrez_names:
        a = Blast_Query_Get_Entrez(query, ez, evalue, max_hits)
        list_specifics.append(a)
    return list_specifics

def Convert_Blast_Result_To_Fasta(list_of_files_to_convert):
    #blast results will be in
    MYFASTA = Fasta(item)
    MYFASTA.gen_raw_blast_lists(item)
    MYFASTA.gen_new_fasta(item+"_converted.fasta")
    return item+"_converted.fasta"


#this is only used in Run_Specific_Blast
class Fasta:
    #
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
        #???
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

        #???
    def gen_species_blast(self):
        self.species_names = []
        #add each thing - the bit between start and | in each id.
        for item in self.ids:
            nam = re.sub("([^\|]*)(\|)(.*)", "\\1", item)
            nam.strip()
            self.species_names.append(nam)
            #


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
                if "[" in species_as_found_in_gsl:
                    species_as_found_in_gsl = ""
                species_as_found_in_gsl = re.sub("[\[\]]", "",  species_as_found_in_gsl)
                species_name = species_as_found_in_gsl
                if species_name == "":

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

                species_name = re.sub("[\[\]:;=,/\+'\.\-\(\)]", "_", species_name)
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

                #
    def single_shorten(self):
        unk = "no"
        normal = 0
        ucount = 0
        for i in range(len(self.ids)):
            line = self.ids[i]
            index = i
            print(line)
            seplist = line.split("|")
            newid = "gi#|"+seplist[1]
            self.ids[index] = newid
    # def blast2fasta(self, blastlist, ENTREZ=False, num=False):
    #     #returns "bad" if bad. else returns True
    #     # entrez is used to ensure that sequence saved uses correct TAXON, esp. if sequence is a MULTISPECIES entry.
    #     # entrex should be somethin like "Mycobacterium triplex"
    #     # num is how many sequences to write. for species trees, we almost certainly only want one.
    #     # for converting full downloaded .fastas, we will want all of them (default = False means to do all of them)
    #     # Converts blast outfmt "6 sseqid stitle sseq" to original lists if
    #     # entrez = false

    #     #... now converting outfmt "6 sallseqid salltitles sseq" to sh fasta with selection of proper gi/acc/taxon
    #     # this should take format " " blast names and replace them with the proper
    #     # fasta shit
    
    #     # we open each file in a unique call to blast2fasta. files should be
    #     # deleted afterwards.
    #     #print(blastlist)
    #     #print(os.getcwd())
    #     bf = open(blastlist, 'r')
    #     bf2 = open(blastlist, 'r')
    #     #print(bf)
    #     error = 0
    #     end = "no"
    #     bad = "no"
    #     length_bf = 0
    #     for line in bf2:
    #         length_bf+=1
    #     errnum = length_bf
    #     #print("Maxerror: "+str(length_bf))
    #     run = 0
    #     #this should be 5 if called from current setup, but might be less if less hits are found.
    #     for line in bf:
    #         run +=1
    #         if end == "yes":
    #             break
    #         # gi|738518257|ref|WP_036466735.1|;gi|620038207|emb|CDO87046.1|   50S
    #         # ribosomal protein L15 [Mycobacterium triplex]<>50S ribosomal protein L15
    #         # [Mycobacterium triplex]
    #         #i just removed a ] from group3/... used to be (.*]) but occasionally species name isn't given at end of thing causes errors.
    #         gis = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\1", line)
    #         names = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\3", line)
    #         seq = re.sub("(.*)(\t)(.*)(\t)([A-Z-]*)", "\\5", line)
    #         # this removes sequences with no Species_name given, so as to avoid errors
    #         # downstream
    #         if "\t" in gis:
    #             error += 1
    #             #print(error)
    #             print("ERROR in blast parsing: " + line)
    #             continue
    #         #check if entrez is in the thing
    #         else:
    #             #THIS IS GOING TO ERROR A GOOD AMOUNT I THINK. MAYBE AVOID USING THIS
    #             gilist = gis.split(";")
    #             namelist = names.split("<>")
    #             if ENTREZ is False:
    #                 index = 0
    #             else:
    #                 ENTREZ = ENTREZ.strip("\"")
    #                 found = "no"
    #                 for item in namelist:
    #                     if ENTREZ+"]" in item:
    #                         index = namelist.index(item)
    #                         found = "yes"
    #                 if found == "no":
    #                     #entrez not found in this thing
    #                     error += 1
    #                     #print("Name error... might fix")
    #                     print(namelist)
    #                     if error == errnum:
    #                         print("Serious ENTREZ error:")
    #                         print(ENTREZ)
    #                         print("This gene wasn't found in this taxon, skipping")
    #                         bad = "yes"
    #                         break
    #                     continue
                        
                        
    #             try:
    #                 seqi = gilist[index].strip() + namelist[index].strip()
    #             except UnboundLocalError:
    #                 #print(error)
    #                 #print("Name error... might fix")
    #                 if error == errnum:
    #                     print("Serious ENTREZ error:")
    #                     print(ENTREZ)
    #                     print(namelist)
    #                     print("This gene wasn't found in this taxon, skipping")
    #                     bad = "yes"
    #                     break
    #                 continue
    #             except:
    #                 #print(index)
    #                 error += 1
    #                 print("unknown index error")
    #                 if error == errnum:
    #                     break
    #                 continue
    #             #print("passed unbound error")
    #                 # goes to next line, abandoning this one
                
    #             # strips for .fasta format
    #             seqid = re.sub("[ ]", "_", seqi)
    #             seqid = seqid.strip()
    #             seqid = seqid.strip(">")
    #             #test that the identified species is the same as you want it to be... if not, skip and try again.
    #             #because some dummy likes to name "gene like in [arabidopsis thaliana] [mus musculus]"
    #             species_as_found_in_gsl = re.sub("([^\[]*)(.*)", "\\2", seqid)
    #             species_as_found_in_gsl = re.sub("[\[\]]", "",  species_as_found_in_gsl)
    #             species_as_found_in_gsl = re.sub("[_]", " ",  species_as_found_in_gsl)
    #             if ENTREZ is False:
    #                 pass
    #             else:
    #                 if ENTREZ == species_as_found_in_gsl:
    #                     pass
    #                 else:
    #                     #print(ENTREZ)
    #                     #print(species_as_found_in_gsl)
    #                     #print(seqid)
    #                     print("Odd, species identified did not match Entrez. Trying next seq, might fix.")
    #                     error += 1
    #                     #print(error)
    #                     end = "no"
    #                     if error == errnum:
    #                         print("Serious ENTREZ error:")
    #                         print(ENTREZ)
    #                         print(namelist)
    #                         print("This gene not found in this taxon, skipping")
    #                         print("Trying to label as non-existant to prevent re-searching later... ")
    #                         bad = "yes"
    #                         break
    #                     continue
               
    #             #print("passed entrez error")
    #             # add the new sequence id to the list.
    #             self.ids.append(seqid)
    #             self.original_ids.append(seqid)
    #             # the new sequence
    #             slist = []
    #             count = 0
    #             newseq = ""
    #             for letter in seq:
    #                 if count > 79:
    #                     count = 0
    #                     newseq = newseq + ("\n")
    #                 newseq = newseq + letter
    #                 count += 1
    #             end = "yes"
    #             #print(newseq)
    #             #print(seqid)
    #             self.seqs.append(newseq.strip())
    #             self.original_seqs.append(newseq.strip())
    #     #print(error)
    #     #print(run)
    #     #how do you manage to get to here, without either being set "Bad" or set good...
    #     if bad == "yes":
    #         #print("Has been defined as bad")
    #         entrez = ENTREZ.replace(" ", "_")
    #         self.ids.append(entrez+"|gi#|000000000")
    #         self.original_ids.append(entrez+"|gi#|000000000")
    #         self.original_seqs.append("----")
    #         self.seqs.append("----")
    #         return "bad"
    #     else:
    #         #print("returning now")
    #         return True
    
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
        #
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
        #
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
            #####new bit added nov 16: no longer requires XMLLINT, does require minidom, urllib2 (probably are standard)
            #print(item)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+item+"&retmode=xml" # define XML location w. current taxid
            #print(url)
        
            try:
                dom = minidom.parse(urllib2.urlopen(url)) # parse the data
            except:
                try:
                    time.sleep(10)
                    dom = minidom.parse(urllib2.urlopen(url)) 
                except:
                    print("triple couldnt parse url: "+url)
                    self.taxid.append("NA")
                    print("error with"+item) 
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
        print(total)
        print(len(self.taxid))
            #
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
           # printProgressBar(current, total, prefix = 'Progress:', suffix = 'Complete', length = 100)
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
        #
    def gen_numbers(self):
        print("Finding identification numbers...")
        for item in self.ids:
        
            if "gi#" in item:
              

                number = re.sub("(.*)(\|)(.*)(\|)(.*)","\\5", item)
         
         
            else:
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
            #
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

def gi_to_species_name(gi):
    #this takes a gi number and gets from it a species name, by first getting a taxid and then taxid -> name + formatting.
    gi = gi.strip()
    taxid = gi_to_taxid(gi)
    if taxid == "error":
        return "error"
    scientific_name = taxid_to_scinam(taxid)
    return scientific_name

def gi_to_taxid(gi):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+gi+"&retmode=xml" # define XML location w. current taxid
    #print(url)
    #b = urllib2.urlopen(url)
    #print(b)
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
           

if __name__ == "__main__":

    print("Running in terminal")
    print("you need to have an internet connection and blastp installed locally. blasts against remote nr database.")
    import sys
    import os
    import argparse
    import re
    from xml.dom import minidom
    import urllib2

    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("-s", "--species_names", action="store", default=False, help ="Type species_file name")
    parser.add_argument("-e", "--entrez", action="store", default=False, help ="Type entrez, or list of entrez in quotes")
    parser.add_argument("-q", "--query", action="store", help="Type query_file name or list of names in quotes")
    parser.add_argument("-d", "--directory", action="store", default=os.getcwd(), help="type path to your query&species files. results will be go here.")
    parser.add_argument("-c", "--cutoff_evalue", action="store", default="1e-5", help="e-value cutoff. defaut 1e-5")
    parser.add_argument("-m", "--max_hits", action="store", default="30000", help="max_hits to download. default 30,000 may take a while. doesn't apply to species_names mode")
    args = parser.parse_args()

    #change dir
    os.chdir(args.directory)

    #set up specifics
    #must be one species_name per line of species_names_file.
    #USE EITHER QUOTES OR UNDERSCORES!!!!!!
    species_names = []
    if args.species_names != False:
        with open(args.species_names) as sp:
            for line in sp:
                line = line.strip()
                species_names.append(line)
    else:
        species_names = False

    #set up entrez
    if args.entrez != False:
        entrez_list = args.entrez.split()
    else:
        entrez_list = False

    #set up queries
    query_list = args.query.split()

    #run the blast
    list_of_output_blasts = Run_Multiple_Queries(query_list, species_names, entrez_list, args.cutoff_evalue, args.max_hits)
    #if this was a specific run, you're done. 
    if args.species_names != False:
        print("finished! outputs are at "+args.directory+" : ")
        print(list_of_output_blasts)
    #else convert the output files to shortened form.
    else:
        list_outputs = []
        #turn to fasta
        for item in list_of_output_blasts:
            newname = Convert_Blast_Result_To_Fasta(item)
            list_outputs.append(newname)
        print("finished! outputs are at "+args.directory+" : ")
        print(list_outputs)

