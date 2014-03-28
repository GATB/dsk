import sys,re

print "Memory\tDisk\tTime\tPasses\tPartitions\tMaxVolumeOfDistinctKmersPerPart"
file=open(sys.argv[1])
result=""
max_mem=0
corrupt_log=False
for line in file:
    s=line.split()
    try:
        if "Sequentially" in line:
            result=s[s.index("memory")-3][1:]+"\t"+s[s.index("disk")-3][1:]+"\t"
            token = "partition(s)" if "partition(s)" in s else "partitions"
            nb_partitions=s[s.index(token)-1]
        if  "Aborted" in line:
            result+="crashed (out of memory)"
        if  "Too many open" in line:
            result+="crashed (too many open files)"
        if  "disk full" in line:
            result+="crashed (disk full)"
        if  "Pass" in line:
            nb_passes=re.findall('[\d]+',line)[1]
        if  "distinct" in line:
            max_mem=max(max_mem,int(int(s[s.index("distinct)")-3])*64.0/8.0/1024.0/1024.0))
        if "real" in line:
            if "crashed" not in result:
                result +=line.split()[1]
                print result+"\t"+str(nb_passes)+"\t"+str(nb_partitions)+"\t"+str(max_mem)
            else:
                print result
            result=""
            max_mem=0
    except:
        corrupt_log=True
if corrupt_log:
    print "BTW, log was corrupt"
