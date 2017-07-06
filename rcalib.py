import re
def replace(filen, pattern, subst):
    # Read contents from file as a single string
    file_handle = open(filen, 'r')
    file_string = file_handle.read()
    file_handle.close()

    # Use RE package to allow for replacement (also allowing for (multiline) REGEX)
    file_string = (re.sub(pattern, subst, file_string))
##    print(file_string)
    # Write contents to file.
    # Using mode 'w' truncates the file.
    file_handle = open(filen, 'w')
    file_handle.write(file_string)
    file_handle.close()

def MakeRCal(fasta):
    nlist = []
    maxlist = []
    minlist = []
    with open(fasta) as old:
        for line in old:
            mini = re.sub("([0-9.]*)(.*)(\n)", "\\1", line)
            maxi = re.sub("([0-9.]*)(\s)([0-9.]*)(.*)(\n)", "\\3", line)
            
            node = re.sub("([0-9.]*)(\s)([0-9.]*)(\s)([0-9]*)(.*)(\n)", "\\5", line)
            nlist.append(node+",")
            maxlist.append(maxi+",")
            minlist.append(mini+",")
    with open(fasta, "a") as app:
        app.write("node=c(")
        for no in nlist:
            app.write(no)
        app.write(", , age.min=c(")
        for no in minlist:
            app.write(no)
        app.write(", , age.max=c(")
        for no in maxlist:
            app.write(no)
        app.write(", ")
        app.flush()
        app.close()
    replace(fasta, ",, ", ")")
    
##MakeRCal("C:/Users/Abby/Dropbox/Phylo/Clock_Analysis/FinalDataset/Current/CulledDataset/SQMO/SecondRun/coelcalib.txt")
##        


if __name__ == "__main__":
    print("CalFile should be a tab separated list:\nminage\tmaxage\tnode#\n")
    import sys
    import argparse
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("calfile", action="store", help ="list file to be calibrated.")
    parser.add_argument("-p", "--printf", action="store_true", help = "prints file after making calib block")
    args = parser.parse_args()

MakeRCal(args.calfile)
if args.printf:
    with open (args.calfile) as this:
        for line in this:
            print (line)


    
