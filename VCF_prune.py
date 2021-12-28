# -*- coding: utf-8 -*-
import os, sys, subprocess, argparse, transposer, random, numpy, csv, scipy, gzip

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(usage ='''
VCF_prune.py [<args>]

The VCF_prune.py args are:
  --inVcfs        STR     path to vcfs (required)
  --winSize       INT     size of scaffold window (required, if --regions is not specified)
  --winDist       INT     distance between any 2 windows (required, if --regions is not specified)
  --regions       STR     path to file defining regions from which to take SNPs (1 SNP per region) format: "CHR:FROM-TO" (required, if --winSize and --winDist is not specified)

  --missing        FLOAT   amount of missing data to allow per site (between 0-1; required)
  --minf     FLOAT   minimum allele frequency, ALT or REF (between 0-1) [0]
  --maxf     FLOAT   maximum allele frequency, ALT or REF (between 0-1) [1]
  --minSnps  INT     Minimal amount of SNPs in window. If less, window will be skipped [1]

  --prefix        STR     output prefix (required)
  --reps        INT     number of replicate data sets [1]
  
  --popFlagLength        INT     length of population name (required)
  --ploidy        INT     ploidy of output Structure file (required)

  --subsample                subsample polyploid data to create psuedo-diploid data
  --gz               use if vcfs are gzipped
  --vcf              use if you want to print also pruned VCF files 
''')
parser.add_argument('--inVcfs', help = 'path to vcfs',
                    type = str, metavar = 'vcf_path',
                    required = True)
parser.add_argument('--winSize', help = 'size of scaffold window',
                    type = int, metavar = 'Window_Size', 
                    required = True)
parser.add_argument('--winDist', help = 'distance between any 2 windows',
                    type = int, metavar = 'Window_Distance',
                    required = True)
parser.add_argument('--reps', help = 'Number of replicate data sets', 
                    type = int, metavar = 'Number_Replications', 
                    required = False, default = '1')
parser.add_argument('--missing', help = 'amount of missing data to allow per site (between 0-1)', 
                    type = float, metavar = 'Missing_Data',
                    required = True)
parser.add_argument('--minf', help='minimum allele frequency, ALT or REF (between 0-1)', 
                    type = float, metavar='Minimum_ALT_REF_Frequency', 
                    required = False, default = '0')
parser.add_argument('--maxf', help = 'maximum allele frequency, ALT or REF (between 0-1)',
                    type = float, metavar = 'Maximum_ALT_REF_Frequency',
                    required = False, default = '1')
parser.add_argument('--minSnps', help = 'minimal amount of SNPs in window. If less, window will be skipped.',
                    type = int, metavar = 'Minimal amount of SNPs in window. If less, window will be skipped.',
                    required = False, default = '1' )
parser.add_argument('--prefix', help = 'Vcfs retain original scaffold name but the concatenated Structure input file will be a text file with specified by output and within the VCF_Pruned directory',
                    type = str, metavar = 'Output_Prefix', 
                    required = False, default = '')
parser.add_argument('--popFlagLength', help = 'length of population name', 
                    type = int, metavar = 'PopFlag_Length',
                    required = False, default = '0' )
parser.add_argument('--ploidy', help='ploidy of Structure file',
                    type=int, metavar='Ploidy',
                    required = True)
parser.add_argument('--subsample', help = 'if true, this will subsample polyploid data to create psuedo-diploid data',
                    action="store_true")
parser.add_argument('--gz', help='are vcfs gzipped (true) or not (false)',
                    action="store_true")
parser.add_argument('--vcf', help='if used, pruned VCF files will be printed',
                    action="store_true")




#output population as 2nd column (must be integer for STRUCTURE*)
#

args = parser.parse_args()

if args.inVcfs.endswith("/") is False:
    args.inVcfs += "/"
if os.path.exists(args.inVcfs + 'VCF_Pruned/') == False: #Create folder for output if it doesn't already exist
    os.mkdir(args.inVcfs + 'VCF_Pruned/')
vcf_list = []

for file in os.listdir(args.inVcfs): #get names of vcf files in args.inVcfs directory
    if args.gz:
        if file[-3:] == '.gz':
            vcf_list.append(file)
    else:   
        if file[-3:] == 'vcf':
            vcf_list.append(file)


count = 0   


def TestSnpQuality(colsToProcess):
    info = colsToProcess[7].split(";")
    AN = float(info[2].split("=")[1])
    AC = float(info[0].split("=")[1])
    currentAlleles = colsToProcess[9:]
    missingData = sum(map(lambda x : x[0] == '.', currentAlleles))
                                                                        # Missing_Data            # min. ALT allele frequency     # MAX. ALT allele frequency     # min. ALT allele frequency           # MAX. ALT allele frequency
    result = colsToProcess[6] == 'PASS' and (float(missingData) / len(currentAlleles)) <= float(args.missing) and AC / AN >= float(args.minf) and AC / AN <= float(args.maxf) and (AN-AC) / AN >= float(args.minf) and (AN-AC) / AN <= float(args.maxf)    
    return result


def ProcessCurrentWindow(current_window, vcf_sites, current_lines):
    print '\n  SNPs in window: ', len(current_window),
    if len(current_window) < args.minSnps: # SKIP
        print '  too few SNPs to select any',
        current_window, current_lines = [], []
        return current_window, vcf_sites, current_lines
    
    for rep in range(int(args.reps)):
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

            genos.extend(-9 for x in range(args.ploidy - len(geno)))

            if args.subsample:
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

        if args.subsample:
            exec('subtempfile' + str(rep+1) + '.write("\t".join(str(item) for item in dgenos))')
            exec('subtempfile' + str(rep+1) + '.write("""\n""")')
            # subtempfile.write('\t'.join(str(item) for item in dgenos))
            # subtempfile.write("\n")

        if args.vcf:
            exec('newVCF' + str(rep+1) + '.write(current_lines[rx])')
            # newVCF.write(current_lines[rx])

    vcf_sites+=1
    current_window, current_lines = [], []
    return current_window, vcf_sites, current_lines

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

#for rep in range(int(args.reps)): 

for rep in range(int(args.reps)):
    exec('markernames' + str(rep+1) + ' = []')
tot_sites = 0

for rep in range(int(args.reps)): 
    exec('structtempfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'.rep'+str(rep+1)+'.VCF_Pruned.TransposedStruct.txt",' + ' \'w\', 0)')   # VYHOD ", 0"
    exec('subtempfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'.rep'+str(rep+1)+'.VCF_Pruned.TransposedStructSubSample.txt",' + ' \'w\', 0)')  # VYHOD ", 0"
    exec('structfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'.StructureInput.rep'+str(rep+1)+'.VCF_Pruned.txt",' + ' \'w\', 0)')   # VYHOD ", 0"
    exec('structfile' + str(rep+1) + '.write("\t")')
    # structtempfile= open(args.inVcfs+ 'VCF_Pruned/'+args.prefix+".rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt",'w')
    # subtempfile= open(args.inVcfs+ 'VCF_Pruned/'+args.prefix+"               .rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt",'w')
    # structfile= open(args.inVcfs+ 'VCF_Pruned/'+args.prefix+".StructureInput.rep"+str(rep+1)+".VCF_Pruned.txt",'w')

# structfile.write("Ind\tPopFlag\t")

if args.subsample: #Create files if subset is true.
    for rep in range(int(args.reps)): 
        exec('subfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'.StructureInput.rep'+str(rep+1)+'.VCF_Pruned.Diploidized.txt",' + ' \'w\')')
        exec('subfile' + str(rep+1) + '.write("\t")')
        # subfile= open(args.inVcfs+ 'VCF_Pruned/'+args.prefix+".StructureInput.rep"+str(rep+1)+".VCF_Pruned.Diploidized.txt",'w')
        # subfile.write("\t")

for iii,vcf in enumerate(vcf_list):
    print vcf,
    if args.gz:
        if args.vcf:
            for rep in range(int(args.reps)): 
                exec('newVCF' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/' +vcf[:-6] + "rep" + str(rep+1) + '.VCF_Pruned.vcf",' + ' \'w\')')
            # newVCF = open(args.inVcfs + 'VCF_Pruned/' +vcf[:-6] + "rep" + str(rep+1) + ".VCF_Pruned.vcf", "w") # new vcf if gzipped previously 
        src = gzip.open(args.inVcfs + vcf)
    else:
        if args.vcf:
            for rep in range(int(args.reps)): 
                exec('newVCF' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/' +vcf[:-3] + "rep" + str(rep+1) + '.VCF_Pruned.vcf",' + ' \'w\')')
            # newVCF = open(args.inVcfs+ 'VCF_Pruned/'+vcf[:-3]+"rep"+str(rep+1)+".VCF_Pruned.vcf",'w') # new vcf for storing info from the random draws
        src = open(args.inVcfs + vcf)

    current_window, current_lines, chosen_lines, names = [], [], [], []
    site_num, start, end  = 0, 0, 0 - int(args.winDist) # end is seh as negative value, due to check: position >= end + int(args.winDist)   -  to catch first snp of vcf file
    first_site=True
    vcf_sites = 0
    
    # evaluate contents of each line of input file
    for line_idx, line in enumerate(src): #Cycle over lines in the VCF file
        cols = line.replace('\n', '').split('\t')  #Split each line of vcf
        if len(cols) < 2:               ## This should be info just before header
            if args.vcf:
                for rep in range(int(args.reps)): 
                    exec('newVCF' + str(rep+1) + '.write(line)')
                # newVCF.write(line)
        elif cols[0] == "#CHROM": #This should be header
            if args.vcf:
                for rep in range(int(args.reps)): 
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
                for rep in range(int(args.reps)):
                    exec('structtempfile' + str(rep+1) + '.write("\t".join(item for item in names for i in range(args.ploidy)))')
                    exec('structtempfile' + str(rep+1) + '.write("""\n""")')
                    # structtempfile.write("\t".join(item for item in names for i in range(args.ploidy)))
                    # structtempfile.write("\n")

                # if popFlagLength > 0; Write PopFlag for each individual which is a unique integer for each population.
                if (args.popFlagLength > 0):
                    if len(names)>0:
                            oldname=names[1][:args.popFlagLength]
                    popcount=1
                    for j,item in enumerate(names):  
                        if oldname!=item[:args.popFlagLength]:
                            oldname=item[:args.popFlagLength]
                            popcount+=1
                        for rep in range(int(args.reps)):
                            exec('structtempfile' + str(rep+1) + '.write(str("%s\t" % str(popcount)) * args.ploidy)')
                        #structtempfile.write(str("%s\t" % str(popcount)) * args.ploidy)
                    for rep in range(int(args.reps)):
                        exec('structtempfile' + str(rep+1) + '.write("""\n""")')
                    #structtempfile.write("\n") 

                    if args.subsample:  #Same task as steps immediately above, but adjusted to accomodate subsetting.  
                        for rep in range(int(args.reps)):
                            exec('subtempfile' + str(rep+1) + '.write("\t".join(item for item in names for i in range(2)))')
                            exec('subtempfile' + str(rep+1) + '.write("""\n""")')
                        # subtempfile.write("\t".join(item for item in names for i in range(2)))
                        # subtempfile.write("\n")
                        popcount = 1
                        
                        if len(names)>0:
                                oldname=names[1][:args.popFlagLength]
                        for j,item in enumerate(names):
                            if oldname!=item[:args.popFlagLength]:
                                oldname=item[:args.popFlagLength]
                                popcount+=1
                            for rep in range(int(args.reps)):
                                exec('subtempfile' + str(rep+1) + '.write(str("%s\t" % str(popcount)) * 2)')
                                # subtempfile.write(str("%s\t" % str(popcount)) * 2)
                        for rep in range(int(args.reps)):
                            exec('subtempfile' + str(rep+1) + '.write("""\n""")')
                            # subtempfile.write("\n")

                first_site=False

                #if first line of vcf passes filters, we go ahead and use it
                if TestSnpQuality(cols):
                    current_window.append(cols)
                    current_lines.append(line)
                    start=position
                    end = position + args.winSize # MK: END je uz zaciatok dalsieho okna
                    # line_num += 1  # MK NACO???
                    site_num = 1

            #All lines caught by this statement are within the current window
            elif position >= start and position < end and site_num!=0:
                if TestSnpQuality(cols):
                    current_window.append(cols)
                    current_lines.append(line)
                    site_num += 1

            #Here, we have moved past the current window and now need to select a site from all sites within the current window before resetting window bounds and moving on to next window.
            #if first line of vcf does not pass filter, then this statement will catch the first line that does
            elif position >= end:  # MK neries kraviny, ak prelezie end, tak vyber jeden ztych OK  # MK: END je uz zaciatok dalsieho okna

                if len(current_window) !=0: # MK ak tam nieco je, tak z toho vyber

                    current_window, vcf_sites, current_lines = ProcessCurrentWindow(current_window, vcf_sites, current_lines)
                    site_num += 1

                if TestSnpQuality(cols) and position >= end + int(args.winDist): # ak SNP ktory mame je odst daleko a vyhovuje, tak ho spracuj
                    start = position
                    end = position + args.winSize

                    current_window.append(cols)
                    current_lines.append(line)
                    site_num += 1
    
    if len(current_window) !=0: # MK ak tam nieco ostalo, tak z toho vyber
        current_window, vcf_sites, current_lines = ProcessCurrentWindow(current_window, vcf_sites, current_lines)
    
    if args.vcf:
        for rep in range(int(args.reps)): 
                exec('newVCF' + str(rep+1) + '.close()')
        # newVCF.close()

    print '\n  sites for vcf: ', vcf_sites, '\n'
    tot_sites = tot_sites + vcf_sites



#This section transposes the temporary files created during this replicate and adds headers, so that they are formatted according to structure
for rep in range(int(args.reps)):
    exec('structtempfile' + str(rep+1) + '.close()')
    exec('subtempfile' + str(rep+1) + '.close()')
    # structtempfile.close()
    # subtempfile.close()

    #Transposes file
    jj=transposer.transpose(i=args.inVcfs+ 'VCF_Pruned/'+args.prefix+".rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt",d="\t",)

    #Write header names for each marker
    exec('structfile' + str(rep+1) + '.write("\t".join(str(marker) for marker in markernames' + str(rep+1) + '))')
    # structfile.write('\t'.join(str(marker) for marker in markernames))
    
    exec('structfile' + str(rep+1) + '.write("""\n""")')
    # structfile.write("\n")
    exec('structfile' + str(rep+1) + '.write("""\n""".join(j for j in jj))')
    # structfile.write('\n'.join(j for j in jj))
    
    exec('structfile' + str(rep+1) + '.close()')
    # structfile.close()
    if args.subsample:
        kk=transposer.transpose(i=args.inVcfs+ 'VCF_Pruned/'+args.prefix+".rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt",d="\t",)
        exec('subfile' + str(rep+1) + '.write("\t".join(str(marker) for marker in markernames' + str(rep+1) + '))')
        # subfile.write('\t'.join(str(marker) for marker in markernames))
        exec('subfile' + str(rep+1) + '.write("""\n""")')
        # subfile.write("\n")
        exec('subfile' + str(rep+1) + '.write("""\n""".join(j for j in jj))')
        # subfile.write('\n'.join(k for k in kk))
        exec('subfile' + str(rep+1) + '.close()')
        # subfile.close()

    #remove the temporary files that contained the info that needed to be transposed
    os.remove(args.inVcfs+ 'VCF_Pruned/'+args.prefix+".rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt")
    os.remove(args.inVcfs+ 'VCF_Pruned/'+args.prefix+".rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt")

        
print 'Finished'



print '\nStructure params:'
print '    NUMINDS  ', len(names)
print '    NUMLOCI  ', tot_sites
print '    PLOIDY   ', args.ploidy
print '    LABEL     1'
print '    POPDATA  ', ("1" if args.popFlagLength > 0 else "0")

