#/home/linhx/bin/local/bin/python
import sys
import re
import os 

assert len(sys.argv) == 6, "USAGE: barcode_seq fastq_1.gz fastq_2.gz output_fastq1.gz output_fastq2.gz".format(sys.argv[0])
# this program may have small part duplicate reads in two sample, as two barcode can also search about 0.000175, but just a small amount, use this program more fast(can parrallle) 
# sometimes theory can't seperate, read1 with one barcode read2 with another
# this program will miss 1% reads can't find barcode
# no consider revser compliment of the barcode seq, as this proram can indentify 99% reads, if exists the situtation that need consider the rev-com, the indentify rate should not over 50%.
# trim barcode have been used
barcode_seq = sys.argv[1] 
barcode_len = 6

IN = os.popen("gzip -dc " + sys.argv[2])
IN2 = os.popen("gzip -dc " + sys.argv[3])
OT1 = os.popen("gzip  -c > "+ sys.argv[4], "w")
OT2 = os.popen("gzip  -c > "+ sys.argv[5], "w")

# only allow one mismatch
LIMIT = barcode_len  #not allow barcode have mismatch, only 6bp, can't be wrong

flag = 0
readname1 = ""
seq1 = ""
qual1 = ""
readname2 = ""
seq2 = ""
qual2 = ""
for buff in IN:
  buff2 = IN2.readline()
  buff = buff.rstrip()
  buff2 = buff2.rstrip()
  flag += 1
  if flag == 1:
    readname1 = buff
    readname2 = buff2
  if flag == 2:
    seq1 = buff
    seq2 = buff2
  if flag == 4:
    flag = 0
    qual1 = buff
    qual2 = buff2
    max_value1 = 0
    max_pos1 = -1
    max_value2 = 0
    max_pos2 = -1
    if seq1[0:7]  == barcode_seq and seq2[0:7]  == barcode_seq: 
      OT1.write("{}\n{}\n+\n{}\n".format(readname1, seq1[barcode_len:len(seq1)],qual1[barcode_len:len(seq1)]))
      OT2.write("{}\n{}\n+\n{}\n".format(readname2, seq2[barcode_len:len(seq2)],qual2[barcode_len:len(seq2)]))
    else:
 # extend search
 # only consider first 9 bp
      for i in range(0,10 - barcode_len):
        m = 0
        for j in range(0, barcode_len):
          if seq1[j+i] == barcode_seq[j] or seq1[j+i] == "N":
            m += 1
        if max_value1 < m: 
          max_value1 = m
          max_pos1 =  i
          if max_value1 == barcode_len:
            break
      
# second
      for i in range(0,10 - barcode_len):
        m = 0
        for j in range(0,barcode_len):
          if seq2[j+i] == barcode_seq[j] or seq2[j+i] == "N":
            m += 1
        if max_value2 < m: 
          max_value2 = m
          max_pos2 =  i
          if max_value2 == barcode_len:
            break
      if not( max_value1 < LIMIT and max_value2 <LIMIT):
        if max_value1 < LIMIT: #no barcode
          OT1.write("{}\n{}\n+\n{}\n".format(readname1, seq1,qual1))
        else:
          OT1.write("{}\n{}\n+\n{}\n".format(readname1, seq1[max_pos1 + barcode_len:len(seq1)],qual1[max_pos1 + barcode_len:len(seq1)]))
        if max_value2 < LIMIT:
          OT2.write("{}\n{}\n+\n{}\n".format(readname2, seq2,qual2))
        else:
          OT2.write("{}\n{}\n+\n{}\n".format(readname2, seq2[max_pos2 + barcode_len:len(seq2)],qual2[max_pos2 + barcode_len:len(seq2)]))

OT1.close()
OT2.close()
    
    
     
    
