# -*- coding: utf-8 -*-
import os, sys, subprocess, argparse, transposer, random, numpy, csv, scipy, gzip

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(usage ='''
VCF_prune.py [<args>]

The VCF_prune.py args are:
  -v        STR     path to vcfs (required)
  -w        INT     size of scaffold window (required)
  -d        INT     distance between any 2 windows (required)

  -m        FLOAT   amount of missing data to allow per site (between 0-1; required)
  -minf     FLOAT   minimum allele frequency, ALT or REF (between 0-1) [0]
  -maxf     FLOAT   maximum allele frequency, ALT or REF (between 0-1) [1]
  -minSnps  INT     Minimal amount of SNPs in window. If less, window will be skipped [1]

  -o        STR     output prefix (required)
  -r        INT     number of replicate data sets [1]
  
  -p        INT     length of population name (required)
  -n        INT     ploidy of output Structure file (required)

  -s                subsample polyploid data to create psuedo-diploid data
  -gz               use if vcfs are gzipped
  -vcf              use if you want to print also pruned VCF files 
''')
parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to vcfs')
parser.add_argument('-w', type = int, metavar = 'Window_Size', required = True, help = 'size of scaffold window')
parser.add_argument('-r', type = int, metavar = 'Number_Replications', required = False, default = '1', help = 'Number of replicate data sets')
parser.add_argument('-d', type = int, metavar = 'Window_Distance', required = True, help = 'distance between any 2 windows')
parser.add_argument('-m', type = float, metavar = 'Missing_Data', required = True, help = 'amount of missing data to allow per site (between 0-1)')
parser.add_argument('-minf', type = float, metavar='Minimum_ALT_REF_Frequency', required = False, default = '0', help='minimum allele frequency, ALT or REF (between 0-1)')
parser.add_argument('-maxf', type = float, metavar = 'Maximum_ALT_REF_Frequency', required = False, default = '1', help = 'maximum allele frequency, ALT or REF (between 0-1)')
parser.add_argument('-minSnps', type = int, metavar = 'Minimal amount of SNPs in window. If less, window will be skipped.', required = False, default = '1', help = 'minimal amount of SNPs in window. If less, window will be skipped.')
parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original scaffold name but the concatenated Structure input file will be a text file with specified by output and within the VCF_Pruned directory')
parser.add_argument('-p', type = int, metavar = 'PopFlag_Length', required = True, help = 'length of population name')
parser.add_argument('-n', type=int, metavar='Ploidy', required = True, help='ploidy of Structure file')
parser.add_argument('-s',  action="store_true", help = 'if true, this will subsample polyploid data to create psuedo-diploid data')
parser.add_argument('-gz', action="store_true", help='are vcfs gzipped (true) or not (false)')
parser.add_argument('-vcf', action="store_true", help='if used, pruned VCF files will be printed')




#output population as 2nd column (must be integer for STRUCTURE*)
#

args = parser.parse_args()

if args.v.endswith("/") is False:
    args.v += "/"
if os.path.exists(args.v + 'VCF_Pruned/') == False: #Create folder for output if it doesn't already exist
    os.mkdir(args.v + 'VCF_Pruned/')
vcf_list = []

for file in os.listdir(args.v): #get names of vcf files in args.v directory
    if args.gz:
        if file[-3:] == '.gz':
            vcf_list.append(file)
    else:   
        if file[-3:] == 'vcf':
            vcf_list.append(file)


count = 0   


def TestSnpQuality():
    info = cols[7].split(";")
    AN = float(info[2].split("=")[1])
    AC = float(info[0].split("=")[1])
    currentAlleles = cols[9:]
    missingData = sum(map(lambda x : x[0] == '.', currentAlleles))
                                                                        # Missing_Data            # min. ALT allele frequency     # MAX. ALT allele frequency     # min. ALT allele frequency           # MAX. ALT allele frequency
    result = cols[6] == 'PASS' and (float(missingData) / len(currentAlleles)) <= float(args.m) and AC / AN >= float(args.minf) and AC / AN <= float(args.maxf) and (AN-AC) / AN >= float(args.minf) and (AN-AC) / AN <= float(args.maxf)    
    return result


def ProcessCurrentWindow(current_window, vcf_sites, current_lines):
    print '  SNPs in window: ',len(current_window)  
    if len(current_window) < args.minSnps: # SKIP
        print '  skipping window'
        current_window, current_lines = [], []
        return current_window, current_lines, vcf_sites

    for rep in range(int(args.r)):
        rx = random.randrange(0,len(current_window))
        site = current_window[rx]
        genos, dgenos = [], []
        
        #Convert alleles
        ref_base = ConvertAllele(site[3])
        alt_base = ConvertAllele(site[4])
                                
        for geno in site[9:]:
            geno=geno.split(":")[0]
            geno=geno.split("/")
            if (len(geno)==1):
                geno=geno[0].split("|")

            for allele in geno:
                if allele == '.':
                    genos.append(-9)
                elif allele == '0':
                    genos.append(ref_base)
                elif allele == '1':
                    genos.append(alt_base)
                else:
                    print("allele not matched")

            genos.extend(-9 for x in range(args.n - len(geno)))

            if args.s:
                for allele in numpy.random.choice(geno,2,replace = False):
                    if allele == '.':
                        dgenos.append(-9)
                    elif allele == '0':
                        dgenos.append(ref_base)
                    elif allele == '1':
                        dgenos.append(alt_base)
                    else:
                        print("allele not matched")
        
        exec('markernames' + str(rep+1) + '.append(str(site[0])+"_"+str(site[1]))')
        # markernames.append(str(site[0])+"_"+str(site[1]))

        exec('structtempfile' + str(rep+1) + '.write("\t".join(str(item) for item in genos))')
        exec('structtempfile' + str(rep+1) + '.write("""\n""")')
        # structtempfile.write('\t'.join(str(item) for item in genos))
        # structtempfile.write("\n")

        if args.s:
            exec('subtempfile' + str(rep+1) + '.write("\t".join(str(item) for item in dgenos))')
            exec('subtempfile' + str(rep+1) + '.write("""\n""")')
            # subtempfile.write('\t'.join(str(item) for item in dgenos))
            # subtempfile.write("\n")

        if args.vcf:
            exec('newVCF' + str(rep+1) + '.write(current_lines[rx])')
            # newVCF.write(current_lines[rx])

    vcf_sites+=1
    current_window, current_lines = [], []
    return current_window, current_lines, vcf_sites

def ConvertAllele(base):
    switcher={
                'A': 1,
                'T': 2,
                'G': 3,
                'C': 4,
             }
    return switcher.get(base,"-99")


###############x
## BEGIN

#for rep in range(int(args.r)): 

for rep in range(int(args.r)):
    exec('markernames' + str(rep+1) + ' = []')
tot_sites = 0

for rep in range(int(args.r)): 
    exec('structtempfile' + str(rep+1) + '= open("' + args.v + 'VCF_Pruned/'+args.o+'.rep'+str(rep+1)+'.VCF_Pruned.TransposedStruct.txt",' + ' \'w\', 0)')   # VYHOD ", 0"
    exec('subtempfile' + str(rep+1) + '= open("' + args.v + 'VCF_Pruned/'+args.o+'.rep'+str(rep+1)+'.VCF_Pruned.TransposedStructSubSample.txt",' + ' \'w\', 0)')  # VYHOD ", 0"
    exec('structfile' + str(rep+1) + '= open("' + args.v + 'VCF_Pruned/'+args.o+'.StructureInput.rep'+str(rep+1)+'.VCF_Pruned.txt",' + ' \'w\', 0)')   # VYHOD ", 0"
    exec('structfile' + str(rep+1) + '.write("\t")')
    # structtempfile= open(args.v+ 'VCF_Pruned/'+args.o+".rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt",'w')
    # subtempfile= open(args.v+ 'VCF_Pruned/'+args.o+"               .rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt",'w')
    # structfile= open(args.v+ 'VCF_Pruned/'+args.o+".StructureInput.rep"+str(rep+1)+".VCF_Pruned.txt",'w')

# structfile.write("Ind\tPopFlag\t")

if args.s: #Create files if subset is true.
    for rep in range(int(args.r)): 
        exec('subfile' + str(rep+1) + '= open("' + args.v + 'VCF_Pruned/'+args.o+'.StructureInput.rep'+str(rep+1)+'.VCF_Pruned.Diploidized.txt",' + ' \'w\')')
        exec('subfile' + str(rep+1) + '.write("\t")')
        # subfile= open(args.v+ 'VCF_Pruned/'+args.o+".StructureInput.rep"+str(rep+1)+".VCF_Pruned.Diploidized.txt",'w')
        # subfile.write("\t")

for iii,vcf in enumerate(vcf_list):
    print vcf
    if args.gz:
        if args.vcf:
            for rep in range(int(args.r)): 
                exec('newVCF' + str(rep+1) + '= open("' + args.v + 'VCF_Pruned/' +vcf[:-6] + "rep" + str(rep+1) + '.VCF_Pruned.vcf",' + ' \'w\')')
            # newVCF = open(args.v + 'VCF_Pruned/' +vcf[:-6] + "rep" + str(rep+1) + ".VCF_Pruned.vcf", "w") # new vcf if gzipped previously 
        src = gzip.open(args.v + vcf)
    else:
        if args.vcf:
            for rep in range(int(args.r)): 
                exec('newVCF' + str(rep+1) + '= open("' + args.v + 'VCF_Pruned/' +vcf[:-3] + "rep" + str(rep+1) + '.VCF_Pruned.vcf",' + ' \'w\')')
            # newVCF = open(args.v+ 'VCF_Pruned/'+vcf[:-3]+"rep"+str(rep+1)+".VCF_Pruned.vcf",'w') # new vcf for storing info from the random draws
        src = open(args.v + vcf)

    current_window, current_lines, chosen_lines, names = [], [], [], []
    site_num, start, end  = 0, 0, 0 - int(args.d) # end is seh as negative value, due to check: position >= end + int(args.d)   -  to catch first snp of vcf file
    first_site=True
    vcf_sites = 0
    
    # evaluate contents of each line of input file
    for line_idx, line in enumerate(src): #Cycle over lines in the VCF file
        cols = line.replace('\n', '').split('\t')  #Split each line of vcf
        if len(cols) < 2:               ## This should be info just before header
            if args.vcf:
                for rep in range(int(args.r)): 
                    exec('newVCF' + str(rep+1) + '.write(line)')
                # newVCF.write(line)
        elif cols[0] == "#CHROM": #This should be header
            if args.vcf:
                for rep in range(int(args.r)): 
                    exec('newVCF' + str(rep+1) + '.write(line)')
                # newVCF.write(line)
            names.extend(cols[9:])
            # for j in cols[9:]: #get names of individuals in vcf
            #     names.append()

        else: 
            position = int(cols[1])
            
            #Convert alleles to Structure input
            alt_base = ConvertAllele(cols[4])
            ref_base = ConvertAllele(cols[3])
            
            if ((alt_base == '-99') or (ref_base == '-99')):
                pass

            elif first_site == True and iii == 0: #Initial setup of info for output files.

                #Write individual name information for the temporary file that is to be transposed.  Names need to be repeated for each observed allele.
                for rep in range(int(args.r)):
                    exec('structtempfile' + str(rep+1) + '.write("\t".join(item for item in names for i in range(args.n)))')
                    exec('structtempfile' + str(rep+1) + '.write("""\n""")')
                    # structtempfile.write("\t".join(item for item in names for i in range(args.n)))
                    # structtempfile.write("\n")

                # Write PopFlag for each individual which is a unique integer for each population.
                if len(names)>0:
                        oldname=names[1][:args.p]
                popcount=1
                for j,item in enumerate(names):  
                    if oldname!=item[:args.p]:
                        oldname=item[:args.p]
                        popcount+=1
                    for rep in range(int(args.r)):
                        exec('structtempfile' + str(rep+1) + '.write(str("%s\t" % str(popcount)) * args.n)')
                    #structtempfile.write(str("%s\t" % str(popcount)) * args.n)
                for rep in range(int(args.r)):
                    exec('structtempfile' + str(rep+1) + '.write("""\n""")')
                #structtempfile.write("\n") 

                if args.s:  #Same task as steps immediately above, but adjusted to accomodate subsetting.  
                    for rep in range(int(args.r)):
                        exec('subtempfile' + str(rep+1) + '.write("\t".join(item for item in names for i in range(2)))')
                        exec('subtempfile' + str(rep+1) + '.write("""\n""")')
                    # subtempfile.write("\t".join(item for item in names for i in range(2)))
                    # subtempfile.write("\n")
                    popcount = 1
                    
                    if len(names)>0:
                            oldname=names[1][:args.p]
                    for j,item in enumerate(names):
                        if oldname!=item[:args.p]:
                            oldname=item[:args.p]
                            popcount+=1
                        for rep in range(int(args.r)):
                            exec('subtempfile' + str(rep+1) + '.write(str("%s\t" % str(popcount)) * 2)')
                            # subtempfile.write(str("%s\t" % str(popcount)) * 2)
                    for rep in range(int(args.r)):
                        exec('subtempfile' + str(rep+1) + '.write("""\n""")')
                        # subtempfile.write("\n")

                first_site=False

                #if first line of vcf passes filters, we go ahead and use it
                if TestSnpQuality():
                    current_window.append(cols)
                    current_lines.append(line)
                    start=position
                    end = position + args.w # MK: END je uz zaciatok dalsieho okna
                    # line_num += 1  # MK NACO???
                    site_num = 1

            #All lines caught by this statement are within the current window
            elif position >= start and position < end and site_num!=0:
                if TestSnpQuality():
                    current_window.append(cols)
                    current_lines.append(line)
                    site_num += 1

            #Here, we have moved past the current window and now need to select a site from all sites within the current window before resetting window bounds and moving on to next window.
            #if first line of vcf does not pass filter, then this statement will catch the first line that does
            elif position >= end:  # MK neries kraviny, ak prelezie end, tak vyber jeden ztych OK  # MK: END je uz zaciatok dalsieho okna

                if len(current_window) !=0: # MK ak tam nieco je, tak z toho vyber

                    current_window, current_lines, vcf_sites = ProcessCurrentWindow(current_window, vcf_sites, current_lines)
                    site_num += 1

                if TestSnpQuality() and position >= end + int(args.d): # ak SNP ktory mame je odst daleko a vyhovuje, tak ho spracuj
                    start = position
                    end = position + args.w

                    current_window.append(cols)
                    current_lines.append(line)
                    site_num += 1
    
    if len(current_window) !=0: # MK ak tam nieco ostalo, tak z toho vyber
        current_window, current_lines, vcf_sites = ProcessCurrentWindow(current_window, vcf_sites, current_lines)
    
    if args.vcf:
        for rep in range(int(args.r)): 
                exec('newVCF' + str(rep+1) + '.close()')
        # newVCF.close()

    print '  sites for vcf: ', vcf_sites, '\n'
    tot_sites = tot_sites + vcf_sites



#This section transposes the temporary files created during this replicate and adds headers, so that they are formatted according to structure
for rep in range(int(args.r)):
    exec('structtempfile' + str(rep+1) + '.close()')
    exec('subtempfile' + str(rep+1) + '.close()')
    # structtempfile.close()
    # subtempfile.close()

    #Transposes file
    jj=transposer.transpose(i=args.v+ 'VCF_Pruned/'+args.o+".rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt",d="\t",)

    #Write header names for each marker
    exec('structfile' + str(rep+1) + '.write("\t".join(str(marker) for marker in markernames' + str(rep+1) + '))')
    # structfile.write('\t'.join(str(marker) for marker in markernames))
    
    exec('structfile' + str(rep+1) + '.write("""\n""")')
    # structfile.write("\n")
    exec('structfile' + str(rep+1) + '.write("""\n""".join(j for j in jj))')
    # structfile.write('\n'.join(j for j in jj))
    
    exec('structfile' + str(rep+1) + '.close()')
    # structfile.close()
    if args.s:
        kk=transposer.transpose(i=args.v+ 'VCF_Pruned/'+args.o+".rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt",d="\t",)
        exec('subfile' + str(rep+1) + '.write("\t".join(str(marker) for marker in markernames' + str(rep+1) + '))')
        # subfile.write('\t'.join(str(marker) for marker in markernames))
        exec('subfile' + str(rep+1) + '.write("""\n""")')
        # subfile.write("\n")
        exec('subfile' + str(rep+1) + '.write("""\n""".join(j for j in jj))')
        # subfile.write('\n'.join(k for k in kk))
        exec('subfile' + str(rep+1) + '.close()')
        # subfile.close()

    #remove the temporary files that contained the info that needed to be transposed
    os.remove(args.v+ 'VCF_Pruned/'+args.o+".rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt")
    os.remove(args.v+ 'VCF_Pruned/'+args.o+".rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt")

        
print 'Finished'



print '\nStructure params:'
print '    NUMINDS  ', len(names)
print '    NUMLOCI  ', tot_sites
print '    PLOIDY   ', args.n
print '    LABEL     1'
print '    POPDATA   1'

