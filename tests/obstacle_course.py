import pysam

print "SAM File Obstacle Course"

header = {'HD': {'VN': '1.0'},
          'SQ': [{'LN': 0001, 'SN': 'gaps1'},
                 {'LN': 0002, 'SN': 'gaps2'},
                 {'LN': 0003, 'SN': 'gaps3'},
                 {'LN': 0004, 'SN': 'normal'},
                 {'LN': 0005, 'SN': 'tumour'},
                 {'LN': 0006, 'SN': 'stats'},
                 {'LN': 0007, 'SN': 'random'},
                 {'LN': 0008, 'SN': 'iupac'},
                 {'LN': 0009, 'SN': 'indels'}]}

outfile = pysam.Samfile("frankensam.sam", "wh", header=header)

# GAPS
# Create three sequences aligned to the same reference but with
# gaps in coverage.  Total depth is 3.

# This case will be used to test joint iterators that must traverse

length = 10
size = 1000
number = 1000
# outer loop: sequences
for seq in [('gaps1', 0, 0), ('gaps2', 4, 1), ('gaps3', (length * 4) + 16, 2)]:
    # inner loop: columns
    print "entering"
    for depth in range(0,3):
        for p in range(0, size, (length * 4)):
            gap = False
            number += 1
            offset = seq[1]
            position = offset + p + (depth * 4)
            qname = "read_00000_00000_" + str(number)
            read = ""
            qual = ""
            mapq = 20
            for i in range(0,length):
                read += "ACGT"

            for i in range(0,length/2):
                qual += "<92.<<<2"

            # If this region is marked for a gap continue
            if (p / (length * 4)) % 5 == 0:
                gap = True
            
            big_gap = seq[2] * 100
            if position > big_gap and position < big_gap + 100:
                gap = True
            
            if gap:
                continue
            a = pysam.AlignedRead() 
            # commit the alignment read info
            a.seq = read
            a.qual = qual
            a.qname = qname
            a.tlen = (length * 4)
            a.flag = 195
            a.rname = seq[2]
            a.mapq = 99
            a.cigar = [(0,length * 4)]
            a.pos = position
            outfile.write(a)
            print position, read, qual

# Column Tests
# Create a 100 x 100 base genome for testing statistical methods

for depth in range(0,100):
    seq = ""
    position = 0
    qual = ""
    a = pysam.AlignedRead()
    for position in range(0,100):
        base = "A"
        if depth - position <= 0:
            if position > 0 and position < 50:
                base = "C"
            else:
                base = "A"

            if depth >= 25:
                base = "C"
            if depth >= 50:
                base = "G"
            if depth >= 75:
                base = "T"
        
        seq += base
        qual += "<"
    print seq
    a.seq = seq
    a.qname = "read_" + str(depth)
    #a.qname = "asdf"
    a.flag = 195
    a.rname = 5
    a.tlen = 100
    a.mapq = 99
    a.cigar = [(0,100)]
    #a.qqual = qual
    outfile.write(a)

reffile = open("frankensam.fa", "w")
reference = ""
for i in range(0,100):

    reference += "A"

print >> reffile, ">stats"
print >> reffile, reference

reference = ""
for i in range (0,300):
    reference += "ACGT" 
print >> reffile, ">gaps1"
print >> reffile, reference

print >> reffile, ">gaps2"
print >> reffile, reference

print >> reffile, ">gaps3"
print >> reffile, reference
