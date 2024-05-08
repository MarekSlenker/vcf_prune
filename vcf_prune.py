# -*- coding: utf-8 -*-
import sys, os, argparse, transposer, random, numpy, gzip, re
from numpy.core.fromnumeric import size

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(usage ='''
VCF_prune.py [<args>]

The VCF_prune.py args are:
  --inVcfs    STR     path to vcfs (required)
  --winSize   INT     size of scaffold window (required, if --regions is not specified)
  --winDist   INT     distance between any 2 windows (required, if --regions is not specified)
  --regions   STR     path to file defining regions from which to take SNPs (1 SNP per region) format: "CHR:FROM-TO" (required, if --winSize and --winDist is not specified)

  --missing   FLOAT   amount of missing data to allow per site (between 0-1) [1]
  --minf      FLOAT   minimum allele frequency, ALT or REF (between 0-1) [0]
  --maxf      FLOAT   maximum allele frequency, ALT or REF (between 0-1) [1]
  --minSnps   INT     Minimal amount of SNPs in window. If less, window will be skipped [1]

  --prefix    STR     output prefix
  --reps      INT     number of replicates to produce (data sets) [1]
  
  --popFlagLength  INT     length of population name
  --ploidy         INT     ploidy of output Structure file [2]

  --subsample        subsample polyploid data to create psuedo-diploid data
  --gz               use if vcfs are gzipped
  --vcf              use if you want to print also pruned VCF files 
''')
parser.add_argument('--inVcfs', help = 'path to vcfs',
                    type = str, metavar = 'vcf_path',
                    required = True)
parser.add_argument('--winSize', help = 'size of scaffold window',
                    type = int, metavar = 'Window_Size', 
                    required = False)
parser.add_argument('--winDist', help = 'distance between any 2 windows',
                    type = int, metavar = 'Window_Distance',
                    required = False)
parser.add_argument('--regions', help = 'path to file defining regions from which to take SNPs (1 SNP per region)',
                    type = str, metavar = 'regions',
                    required = False)                    
parser.add_argument('--reps', help = 'Number of replicate data sets', 
                    type = int, metavar = 'Number_Replications', 
                    required = False, default = '1')
parser.add_argument('--missing', help = 'amount of missing data to allow per site (between 0-1)', 
                    type = float, metavar = 'Missing_Data',
                    required = False, default = '1')
parser.add_argument('--minf', help='minimum allele frequency, ALT or REF (between 0-1)', 
                    type = float, metavar='Minimum_ALT_REF_Frequency', 
                    required = False, default = '0')
parser.add_argument('--maxf', help = 'maximum allele frequency, ALT or REF (between 0-1)',
                    type = float, metavar = 'Maximum_ALT_REF_Frequency',
                    required = False, default = '1')
parser.add_argument('--minSnps', help = 'minimal amount of SNPs in window. If less, window will be skipped.',
                    type = int, metavar = 'Minimal amount of SNPs in window. If less, window will be skipped.',
                    required = False, default = '1')
parser.add_argument('--prefix', help = 'Vcfs retain original scaffold name but the concatenated Structure input file will be a text file with specified by output and within the VCF_Pruned directory',
                    type = str, metavar = 'Output_Prefix', 
                    required = False, default='')
parser.add_argument('--popFlagLength', help = 'length of population name', 
                    type = int, metavar = 'PopFlag_Length',
                    required = False, default = '0' )
parser.add_argument('--ploidy', help='ploidy of Structure file',
                    type=int, metavar='Ploidy',
                    required = False, default = '2')
parser.add_argument('--subsample', help = 'if true, this will subsample polyploid data to create psuedo-diploid data',
                    action="store_true")
parser.add_argument('--gz', help='are vcfs gzipped (true) or not (false)',
                    action="store_true")
parser.add_argument('--vcf', help='if used, pruned VCF files will be printed',
                    action="store_true")
parser.add_argument('--verbose', help='provides additional details',
                    action="store_true")



#output population as 2nd column (must be integer for STRUCTURE*)
#

args = parser.parse_args()

print '# ARGS:', 
for arg in sys.argv:
    print arg,  
print; print; 

if ((args.regions is None) and (args.winSize is not None) and (args.winDist is not None)):
    pass # WGS
elif ((args.regions is not None) and (args.winSize is None) and (args.winDist is None)):
    pass # RADSEQ
else:
    print '\n  ERROR: you have to provide either --winSize and --winDist, OR --regions'
    exit()

if (args.prefix != ''): args.prefix = args.prefix + "."

if args.inVcfs.endswith("/") is False:
    args.inVcfs += "/"
if os.path.exists(args.inVcfs + 'VCF_Pruned/') == False: #Create folder for output if it doesn't already exist
    os.mkdir(args.inVcfs + 'VCF_Pruned/')

if (args.regions is not None): # parse regions
    with open(args.regions,'r') as file:
        lines = [re.split(':|-|\n',line) for line in file]



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
    for str in info:   
        if "AN" in str:
            AN = float(str.split("=")[1])
        if "AC" in str:
            AC = float(str.split("=")[1])
    currentAlleles = colsToProcess[9:]
    missingData = sum(map(lambda x : x[0] == '.', currentAlleles))
                                                                        # Missing_Data            # min. ALT allele frequency     # MAX. ALT allele frequency     # min. ALT allele frequency           # MAX. ALT allele frequency
    result = colsToProcess[6] == 'PASS' and (float(missingData) / len(currentAlleles)) <= float(args.missing) and AC / AN >= float(args.minf) and AC / AN <= float(args.maxf) and (AN-AC) / AN >= float(args.minf) and (AN-AC) / AN <= float(args.maxf)    
    return result

def ProcessCurrentWindow(current_window, vcf_sites, area):
    print '\n  SNPs in ', area, ': ', len(current_window),
    if len(current_window) < args.minSnps: # SKIP
        print '  too few SNPs to select any',
        return vcf_sites
    
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

        exec('structtempfile' + str(rep+1) + '.write("\t".join(str(item) for item in genos))')
        exec('structtempfile' + str(rep+1) + '.write("""\n""")')

        if args.subsample:
            exec('subtempfile' + str(rep+1) + '.write("\t".join(str(item) for item in dgenos))')
            exec('subtempfile' + str(rep+1) + '.write("""\n""")')

        if args.vcf:
            exec('newVCF' + str(rep+1) + '.write(\'\t\'.join(current_window[rx]))')
            exec('newVCF' + str(rep+1) + '.write("""\n""")')

    vcf_sites+=1
    return vcf_sites

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

    exec('structtempfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'rep'+str(rep+1)+'.VCF_Pruned.TransposedStruct.txt",' + ' \'w\', 0)')   # VYHOD ", 0"
    exec('subtempfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'rep'+str(rep+1)+'.VCF_Pruned.TransposedStructSubSample.txt",' + ' \'w\', 0)')  # VYHOD ", 0"
    exec('structfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'StructureInput.rep'+str(rep+1)+'.VCF_Pruned.txt",' + ' \'w\', 0)')   # VYHOD ", 0"
    if (args.popFlagLength > 0):
        exec('structfile' + str(rep+1) + '.write("\t")')
    exec('structfile' + str(rep+1) + '.write("\t")')

tot_sites = 0
first_site=True
names = []

if args.subsample: #Create files if subset is true.
    for rep in range(int(args.reps)): 
        exec('subfile' + str(rep+1) + '= open("' + args.inVcfs + 'VCF_Pruned/'+args.prefix+'StructureInput.rep'+str(rep+1)+'.VCF_Pruned.Diploidized.txt",' + ' \'w\')')
        exec('subfile' + str(rep+1) + '.write("\t")')


for vcf in vcf_list:
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

    
    
    vcf_sites, passedSnps, skippedSnps = 0, 0, 0

    # PREOCESS HEADER
    # evaluate contents of each line of input file
    while True: #Cycle over lines in the VCF file
        line = src.readline()
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

            if first_site == True: #Initial setup of info for output files.
                names.extend(cols[9:])
                #Write individual name information for the temporary file that is to be transposed.  Names need to be repeated for each observed allele.
                for rep in range(int(args.reps)):
                    exec('structtempfile' + str(rep+1) + '.write("\t".join(item for item in names for i in range(args.ploidy)))')
                    exec('structtempfile' + str(rep+1) + '.write("""\n""")')

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
                    for rep in range(int(args.reps)):
                        exec('structtempfile' + str(rep+1) + '.write("""\n""")')

                    if args.subsample:  #Same task as steps immediately above, but adjusted to accomodate subsetting.  
                        for rep in range(int(args.reps)):
                            exec('subtempfile' + str(rep+1) + '.write("\t".join(item for item in names for i in range(2)))')
                            exec('subtempfile' + str(rep+1) + '.write("""\n""")')
                        popcount = 1
                        
                        if len(names)>0:
                                oldname=names[1][:args.popFlagLength]
                        for j,item in enumerate(names):
                            if oldname!=item[:args.popFlagLength]:
                                oldname=item[:args.popFlagLength]
                                popcount+=1
                            for rep in range(int(args.reps)):
                                exec('subtempfile' + str(rep+1) + '.write(str("%s\t" % str(popcount)) * 2)')
                        for rep in range(int(args.reps)):
                            exec('subtempfile' + str(rep+1) + '.write("""\n""")')
                
                first_site=False
            break # skip remaining lines


    # PROCESS REGIONS
    if (args.regions is not None):  
        
        # peek line - to find current chrom
        xpos = src.tell()
        xline = src.readline()
        src.seek(xpos)
        
        chrom = xline.replace('\n', '').split('\t')[0]  # get CHROM

        currentRegions, rowsOfRegions = [], []
        for line in lines:
            if line[0] == chrom:
                currentRegions.append(line)
    
        currentRegsBegins = numpy.array([int(x) for x in numpy.array(currentRegions)[:,1]])
        currentRegsEnds = numpy.array([int(x) for x in numpy.array(currentRegions)[:,2]])
            
        for i in range(0,len(currentRegions)):
            rowsOfRegions.append([])


        for line in src: #Cycle over lines in the VCF file
            cols = line.replace('\n', '').split('\t')  #Split each line of vcf
                
            position = int(cols[1])
            alt_base = ConvertAllele(cols[4])
            ref_base = ConvertAllele(cols[3])
            
            if ((alt_base == '-99') or (ref_base == '-99')):
                continue
            
            rr = numpy.where((position >= currentRegsBegins) == (position <= currentRegsEnds))[0] # find where position is inside the range
            if rr.size == 1:
                if TestSnpQuality(cols):
                    rowsOfRegions[rr[0]].append(cols)
                    passedSnps+=1
                else:
                    skippedSnps+=1
            #if rr.size == 0:
            #    print "!!", chrom, ": ", position, "out of regions"
            #if rr.size > 1:
            #    print "!!", chrom, ": ", position, "in more than one region"                    

        # pocess regions
        for r in range(0, len(rowsOfRegions), 1): 
            regionName = currentRegions[r][0] + ":" + currentRegions[r][1] + "-" + currentRegions[r][2]
            if rowsOfRegions[r]:
                vcf_sites = ProcessCurrentWindow(rowsOfRegions[r], vcf_sites, regionName)
            else:
                if args.verbose:
                    print '\n  SNPs in ', regionName, ': 0; region is empty/invariant',



    # PROCESS USING  winSize and winDist
    if (args.winSize is not None) and (args.winDist is not None): 
        start, end = 0, 0 - int(args.winDist)
        current_window = []

        for line in src: #Cycle over lines in the VCF file
            cols = line.replace('\n', '').split('\t')  #Split each line of vcf
                
            position = int(cols[1])
            alt_base = ConvertAllele(cols[4])
            ref_base = ConvertAllele(cols[3])
            
            if ((alt_base == '-99') or (ref_base == '-99')):
                continue

            #All lines caught by this statement are within the current window
            if position >= start and position < end: #  and site_num!=0
                if TestSnpQuality(cols):
                    current_window.append(cols)
                    passedSnps+=1
                else:
                    skippedSnps+=1

            #Here, we have moved past the current window and now need to select a site from all sites within the current window before resetting window bounds and moving on to next window.
            #if first line of vcf does not pass filter, then this statement will catch the first line that does
            elif position >= end:  # MK neries kraviny, ak prelezie end, tak vyber jeden ztych OK  # MK: END je uz zaciatok dalsieho okna

                if len(current_window) !=0: # MK ak tam nieco je, tak z toho vyber

                    vcf_sites = ProcessCurrentWindow(current_window, vcf_sites, "window")
                    current_window = [] 

                if TestSnpQuality(cols) and position >= end + int(args.winDist): # ak SNP ktory mame je odst daleko a vyhovuje, tak ho spracuj
                    start = position
                    end = position + args.winSize
                    current_window.append(cols)
                    passedSnps+=1
                else:
                    skippedSnps+=1

        if len(current_window) !=0: # MK ak tam nieco ostalo, tak z toho vyber
            vcf_sites = ProcessCurrentWindow(current_window, vcf_sites, "window")
            current_window = [] 


    if args.vcf:
        for rep in range(int(args.reps)): 
                exec('newVCF' + str(rep+1) + '.close()')

    print '\n  sites for vcf: ', vcf_sites
    print '  numb of sites, passed \"--missing\", \"--minf\", \"--maxf\"... : ', passedSnps
    print '  numb of sites, skipped due to  \"--missing\", \"--minf\", \"--maxf\"... : ', skippedSnps, '\n'



    tot_sites = tot_sites + vcf_sites






#This section transposes the temporary files created during this replicate and adds headers, so that they are formatted according to structure
for rep in range(int(args.reps)):
    exec('structtempfile' + str(rep+1) + '.close()')
    exec('subtempfile' + str(rep+1) + '.close()')

    #Transposes file
    jj=transposer.transpose(i=args.inVcfs+ 'VCF_Pruned/'+args.prefix+"rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt",d="\t",)

    #Write header names for each marker
    exec('structfile' + str(rep+1) + '.write("\t".join(str(marker) for marker in markernames' + str(rep+1) + '))')
    
    exec('structfile' + str(rep+1) + '.write("""\n""")')
    exec('structfile' + str(rep+1) + '.write("""\n""".join(j for j in jj))')
    
    exec('structfile' + str(rep+1) + '.close()')
    if args.subsample:
        kk=transposer.transpose(i=args.inVcfs+ 'VCF_Pruned/'+args.prefix+"rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt",d="\t",)
        exec('subfile' + str(rep+1) + '.write("\t".join(str(marker) for marker in markernames' + str(rep+1) + '))')
        exec('subfile' + str(rep+1) + '.write("""\n""")')
        exec('subfile' + str(rep+1) + '.write("""\n""".join(j for j in jj))')
        exec('subfile' + str(rep+1) + '.close()')

    #remove the temporary files that contained the info that needed to be transposed
    os.remove(args.inVcfs+ 'VCF_Pruned/'+args.prefix+"rep"+str(rep+1)+".VCF_Pruned.TransposedStruct.txt")
    os.remove(args.inVcfs+ 'VCF_Pruned/'+args.prefix+"rep"+str(rep+1)+".VCF_Pruned.TransposedStructSubSample.txt")

        
print 'Finished'



print '\nStructure params:'
print '    NUMINDS  ', len(names)
print '    NUMLOCI  ', tot_sites
print '    PLOIDY   ', args.ploidy
print '    LABEL     1'
print '    POPDATA  ', ("1" if args.popFlagLength > 0 else "0")
print
