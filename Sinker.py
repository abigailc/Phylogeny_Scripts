#!/usr/bin/python

# created abigailc@Artemis on november 8 2016

# this script will take a folder, and run all it's contents though muscle and raxml on the engaging cluster.

#set these yourself
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe5.mit.edu"
clus_head = "abigailc@eofe5.mit.edu:/home/abigailc/"

#ssh_inst = "ssh dgruen eofe5.mit.edu"
#clus_head = "dgruen@eofe5.mit.edu:/home/dgruen/"

#imports
import sys
import argparse
import os
import re
import time

#makes dir if need be
def check_directory_existance(prefix, ssh_inst):
    import os
    print("checking dirs")
    os.system(ssh_inst+" \'mkdir "+prefix+"\'")


def gen_correlate_file(list_of_input_files, corr_file):
    #this should be in form
    #1 name
    #2 name
    #3 name
    #requires 1: list of files 2. name for corr_file.
    i = 0
    with open(corr_file, "w") as corr:
        for item in list_of_input_files:
            i += 1
            #make sure the \n spacing works correctly.
            corr.write(str(i)+" "+item+"\n")
    return corr_file



def gen_musrax_array_script(scriptfile, indexname, n, Jobname, AA_MODEL):
    #currently assuming you are running in the dir that files are in and should be returned to.
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 5-00:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                                                                                 
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""


##add the modules
. /etc/profile.d/modules.sh
module add engaging/muscle/3.8.31
module add engaging/RAxML/8.2.9

##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
TAIL=_Muscle.fasta
HEAD=${THE_INPUT_FILE%%.*}
MUSOUT=$HEAD$TAIL
RAXOUT=${MUSOUT%%.*}

muscle -in $THE_INPUT_FILE -out $MUSOUT

raxmlHPC-PTHREADS-AVX -T 20 -f a -m """+AA_MODEL+""" -p 12345 -x 12345 -#100 -n $RAXOUT -s $MUSOUT

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

    

#does everything
def musrax_array_on_cluster(end_file_list, prefix, AA_MODEL):
    #this creates dir you will use on the cluster.
    aligned_list = []
    remove_list = []
    print(end_file_list)
    for item in end_file_list:
        exts = item.split(".")
        tail = exts[0]
        head = "RAxML_bipartitions."
        new = head+tail+"_Muscle"
        aligned_list.append(new)
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/"+prefix
    a = gen_musrax_array_script(prefix+"_Sc.sh", "~"+clus_path+"/"+prefix+"_Corr.txt", str(len(end_file_list)), prefix+"job", AA_MODEL)
    b = gen_correlate_file(end_file_list, prefix+"_Corr.txt")
    end_file_list.append(a)
    end_file_list.append(b)
    direct = os.getcwd()
    move_to_cluster(end_file_list, clus_path)
    print("everything should be generated and on the cluster")
    n = str(len(end_file_list))
    print("there are "+n+" files that should be aligning in muscle and then RAXML-ing right now")
    os.system(ssh_inst+" 'cd ~/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    #initialize list of things to move home. initially equal to aligned_list.
    movehome = []
    time.sleep(180)
   
    for i in aligned_list:
        movehome.append(i)
    while finished is not True:
        #try and move each home.
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        for item in aligned_list:
            #see if it got moved home.
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
                finished = "yes"
            else:
                finished = False
                print(item+" not found")
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait five minutes and then try again.
            print("checking.... some alignment outputs do not exist yet. sleeping for 10 minutes.")
            time.sleep(600)
    #files should        
    print("Your files have been aligned! They are located at "+direct)
    return aligned_list



#moves your .fasta and script/corr to cluster
def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)

# parser

if __name__ == "__main__":
    print("Running in terminal")  
    import sys
    import argparse
    import os
    import re
    import time
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("-p", "--projectname", action = "store", default = "Unnamed_Job", help="give a name for your script/job")
    parser.add_argument("-f", "--folder", action = "store", default = False, help="give a folder name. should contain ONLY .fasta files.")
    parser.add_argument("-ss", "--shorten_species", action = "store_true", default = False, help = "toggles FISH-2 shortening and subsampling one per species")
    parser.add_argument("-sg", "--shorten_genus", action = "store_true", default = False, help = "toggles FISH-2 shortening and subsampling one per genus")
    parser.add_argument("-so", "--shorten_only", action = "store_true", default = False, help = "only shortens, doesn't run on cluster ")
    parser.add_argument("-ns", "--no_shorten", action = "store_true", default = False, help = "only runs on cluster, no shortening ")
    parser.add_argument("-aa", "--aa_model", action = "store_true", default = False, help = "specify model if you hate PROTGAMMALGF ")
  
    
    args = parser.parse_args()
    #create input file list

    #HAHAHAH I DELETED THIS SHIT ACCIDENTALLY SMOL FIX BELOW
    args.keep_gene_info = False
    args.aa_model = False


    Input_List = os.listdir(args.folder)
    if args.keep_gene_info is True:
        short = "-ski"
    else:
        short = "-sh"
    for item in Input_List:
        print("found file :"+item)
        if item[0] == ".":
            print(". remove")
            Input_List.remove(item)
        elif "_Sc.sh" in item:
            print("script remove")
            Input_List.remove(item)
        elif "_Corr.txt" in item:
            print("corr remove")
            Input_List.remove(item)
    print(Input_List)
    #will be [x.fasta, y.fasta]
    
    #for each input, shorten and subsample one per species using FISH_2 ( fish_2 needs to be in same folder as this script)
    if args.shorten_only is True:
        fish_out_list = []
        for item in Input_List:
            #item = x.fasta
            #vis_item = ./folder/x.fasta
            vis_item = "./"+args.folder+"/"+item
            vis_outfile = "./"+args.folder+"/Sh_1pS_"+item
            outfile = "Sh_1pS_"+item
            os.system("python FISH_2.py "+short+" -fas "+vis_item+" -wf "+vis_outfile)
            #flag for shorten to keep text information!!!
            fish_out_list.append(outfile)
        #redefine the input list to be the newly generated bait files.
        Input_List = fish_out_list
        print(Input_List)
        print("Your files have been shortened but not subsampled or sent to the cluster")
        raise SystemExit
    if args.no_shorten is True:
        pass
    elif args.shorten_species is True:
        fish_out_list = []
        for item in Input_List:
            #item = x.fasta
            #vis_item = ./folder/x.fasta
            vis_item = "./"+args.folder+"/"+item
            vis_outfile = "./"+args.folder+"/Sh_1pS_"+item
            outfile = "Sh_1pS_"+item
            os.system("python FISH_2.py -os "+short+" -fas "+vis_item+" -wf "+vis_outfile)
            fish_out_list.append(outfile)
        #redefine the input list to be the newly generated bait files.
        Input_List = fish_out_list
        print(Input_List)
    #there should now be a bunch of sh_1ps_x.fasta in the folder args.folder, and also in the list "fish_out_list"
    elif args.shorten_genus is True:
        fish_out_list = []
        for item in Input_List:
            #item = x.fasta
            #vis_item = ./folder/x.fasta
            vis_item = "./"+args.folder+"/"+item
            vis_outfile = "./"+args.folder+"/Sh_1pG_"+item
            outfile = "Sh_1pG_"+item
            os.system("python FISH_2.py -og -sh -fas "+vis_item+" -wf "+vis_outfile)
            fish_out_list.append(outfile)
        #redefine the input list to be the newly generated bait files.
        Input_List = fish_out_list
        print(Input_List)
    else:
        fish_out_list = []
        for item in Input_List:
            #item = x.fasta
            #vis_item = ./folder/x.fasta
            vis_item = "./"+args.folder+"/"+item
            vis_outfile = "./"+args.folder+"/Sh_"+item
            outfile = "Sh_"+item
            #change this to python FISH_2.py if running elseqherer
            os.system("python FISH_2.py "+short+" -fas "+vis_item+" -wf "+vis_outfile)
            fish_out_list.append(outfile)
        #redefine the input list to be the newly generated bait files.
        Input_List = fish_out_list
        print(Input_List)
        
    #there should now be a bunch of sh_1ps_x.fasta in the folder args.folder, and also in the list "fish_out_list"
    
    #change the directory so we are within args.folder
    try:
        os.chdir(args.folder)
    except:
        print ("didn't change dir")
        raise SystemExit
    if args.aa_model != False:
        AA_MODEL = args.aa_model
    else:
        AA_MODEL = "PROTGAMMALGF"
    #run the thing
    musrax_array_on_cluster(Input_List, args.projectname, AA_MODEL)
    print("done")



