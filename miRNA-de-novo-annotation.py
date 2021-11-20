##miRNA de novo annotation

##subject files:
####1.Results.txt (from ShortStack.pl results)
####2.Counts.txt (from ShortStack.pl results)
####3.Cluster.txt (from : ShortStack.pl results ----cat MIRNA/* > Cluster.txt)
####4.mature_extension.txt (from miRBase)
####5.hairpin_T.fa (from miRBase)
####6.short_miRNA_anno3.pl
####7.short_miRNA_anno4.pl
####8.mature_extension2.pl

##Author: QiYongli22@njfu.edu.cn


import re
import os
from sys import argv


def extract_Results(list_MIRNA):
    if list_MIRNA==['Y']:
        w1=open('shortstack_Y.fasta','w')
    else:
        w1=open('shortstack_N15_Y.fasta','w')
    with open ('Results.txt','r') as a1:
        for line in a1:
            li = line.split()
            Name = li[0]
            sequence = li[9]
            MIRNA = li[13]
            if MIRNA in list_MIRNA:
                w1.write('>' + Name + '\n' + sequence + '\n')
    w1.close()
    if list_MIRNA==['Y']:
        return 'shortstack_Y.fasta'
    else:
        return 'shortstack_N15_Y.fasta'


def RNA2DNA(file_in):
    file_out=file_in.split('.fasta')[0]+'_T.fasta'
    w1=open(file_out, 'w')
    with open(file_in,'r') as a1:
        for line in a1:
            if line.startswith('>'):
                w1.write(line)
            else:
                w1.write(line.replace('U', 'T'))
    w1.close()
    return file_out


def run_patman(file_in,database):
    file_out=file_in+'_out.txt'
    val=os.system('patman -e 4 -g 0 -D '+database+' -P '+file_in+' --singlestrand   -o '+file_out)
    return file_out


def run_perl(file_in):
    file_out=file_in.replace('out','hit')
    val=os.system('perl short_miRNA_anno3.pl '+file_in+' '+file_out)
    return file_out


def run_perl2(file_in):
    file_out=file_in.replace('out','hit')
    val=os.system('perl short_miRNA_anno4.pl '+file_in+' '+file_out)
    return file_out


def name_known(file_in,pre_name):
    a1 = open(file_in, 'r')
    list_letter = []
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append(chr(i))
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append('a'+chr(i))
    list_family = []
    dict_family_count = {}
    dict_loci2family = {}
    dict_loci2name = {}
    for line in a1:
        li = line.split()
        if int(li[5]) + int(li[6]) <= 4:
            family = li[7]
            if family not in list_family:
                i = 1
                dict_loci2family[li[1]] = family
                dict_family_count[family] = i
                list_family.append(family)
                name = pre_name+'-miR' + family + 'a'
                dict_loci2name[li[1]] = name
            else:
                dict_loci2family[li[1]] = family
                i = 1
                i = i + dict_family_count[family]
                dict_family_count[family] = i
                name = pre_name+'-miR' + family + list_letter[i - 1]  # pdl  populus deltoides
                dict_loci2name[li[1]] = name
    a2=open(file_in, 'r')
    for line in a2:
        li = line.split()
        if int(li[5]) + int(li[6]) <= 4:
            if dict_family_count[dict_loci2family[li[1]]] == 1:
                dict_loci2name[li[1]] = pre_name+'-miR' + li[7]
    w1=open('known.txt', 'w')
    for key in dict_loci2name:
        w1.write(key + '\t' + dict_loci2name[key] + '\n')
    w1.close()
    return 'known.txt'


def delect_known_loci(file_in,file_in2,file_out):
    a1=open(file_in,'r')
    a2=open(file_in2,'r')
    list=[]
    for line in a1:
        li = line.split()
        list.append(li[0])
    dict = {}
    for line in a2:
        if line.startswith('>'):
            chr = line.split()[0][1:]
        else:
            dict[chr] = line.split()[0]
    a3=open(file_in2, 'r')
    w1=open(file_out, 'w')
    for line in a3:
        if line.startswith('>'):
            chr = line.split()[0][1:]
            if chr not in list:
                w1.write(line + dict[chr] + '\n')
    w1.close()


def Cluster():
    #val = os.system('cat ../* > Cluster.txt')
    a1=open('Cluster.txt','r')
    a2=open('Cluster.txt','r')
    w1=open('shortstack.txt','w')
    dict_loci2cluster = {}
    dict_loci2hairpin_loci = {}
    dict_loci2hairpin_seq = {}
    dict_loci2hairpin_struct = {}
    dict_loci2strand = {}
    dict_loci2mature_loci = {}
    dict_loci2mature_seq = {}
    dict_loci2star_loci = {}
    dict_loci2star_seq = {}
    dict_loci2star_predict = {}
    for line in a1:
        li = line.split()
        if line.startswith('Cluster'):
            i = 0
            cluster = li[0]
            loci = li[3]
            hairpin_loci = li[6]
            strand = li[8]
            dict_loci2cluster[loci] = cluster
            dict_loci2hairpin_loci[loci] = hairpin_loci
            dict_loci2strand[loci] = strand
        else:
            i = i + 1
        if i == 2:
            hairpin_seq = li[0]
            dict_loci2hairpin_seq[loci] = hairpin_seq
        if i == 3:
            hairpin_struct = li[0]
            dict_loci2hairpin_struct[loci] = hairpin_struct
        if i == 4:
            mature_seq = ''
            hairpin_start = int(hairpin_loci.split(':')[1].split('-')[0])
            mature_start = ''
            for a in li[0]:
                if a != '.':
                    mature_seq = mature_seq + a
                    if len(mature_seq) == 1:
                        mature_start = str(hairpin_start)
                if mature_start != '' and len(mature_seq) == hairpin_start - int(mature_start) + 1:
                    mature_end = str(hairpin_start)
                    mature_loci = mature_start + '-' + mature_end
                    dict_loci2mature_loci[loci] = mature_loci
                hairpin_start = hairpin_start + 1
            dict_loci2mature_seq[loci] = mature_seq
        if i == 5:
            star_seq = ''
            hairpin_start = int(hairpin_loci.split(':')[1].split('-')[0])
            star_start = ''
            for a in li[0]:
                if a != '.':
                    star_seq = star_seq + a
                    if len(star_seq) == 1:
                        star_start = str(hairpin_start)
                if star_start != '' and len(star_seq) == hairpin_start - int(star_start) + 1:
                    star_end = str(hairpin_start)
                    star_loci = star_start + '-' + star_end
                    dict_loci2star_loci[loci] = star_loci
                hairpin_start = hairpin_start + 1
            dict_loci2star_seq[loci] = star_seq
            star_start_number = int(star_start) - int(hairpin_loci.split(':')[1].split('-')[0])
            star_end_numbe = int(star_end) - int(hairpin_loci.split(':')[1].split('-')[0]) + 1
            predict = hairpin_seq[star_start_number:star_end_numbe]
            dict_loci2star_predict[loci] = predict
        if i == 4 and strand == '-':
            mature_end = int(hairpin_loci.split(':')[1].split('-')[1]) - int(mature_start) + int(
                hairpin_loci.split(':')[1].split('-')[0])
            mature_start = mature_end - len(mature_seq) + 1
            mature_loci = str(mature_start) + '-' + str(mature_end)
            dict_loci2mature_loci[loci] = mature_loci
        if i == 5 and strand == '-':
            star_end = int(hairpin_loci.split(':')[1].split('-')[1]) - int(star_start) + int(
                hairpin_loci.split(':')[1].split('-')[0])
            star_start = star_end - len(star_seq) + 1
            star_loci = str(star_start) + '-' + str(star_end)
            dict_loci2star_loci[loci] = star_loci
    # print(dict_loci2cluster['Chr05:6806135-6806155'])
    # print(dict_loci2hairpin_loci['Chr05:6806135-6806155'])
    # print(dict_loci2hairpin_seq['Chr05:6806135-6806155'])
    # print(dict_loci2hairpin_struct['Chr05:6806135-6806155'])
    # print(dict_loci2strand['Chr05:6806135-6806155'])
    # print(dict_loci2mature_loci['Chr05:6806135-6806155'])
    # print(dict_loci2mature_seq['Chr05:6806135-6806155'])
    # print(dict_loci2star_loci['Chr05:6806135-6806155'])
    # print(dict_loci2star_seq['Chr05:6806135-6806155'])
    # print(dict_loci2star_predict['Chr05:6806135-6806155'])
    for line in a2:
        li = line.split()
        if line.startswith('Cluster'):
            loci = li[3]
            w1.write(loci + '\t' + dict_loci2cluster[loci] + '\t' + dict_loci2hairpin_loci[loci].split(':')[0] + '\t' +
                    dict_loci2hairpin_loci[loci].split(':')[1].split('-')[0] + '\t' +
                    dict_loci2hairpin_loci[loci].split(':')[1].split('-')[1] + '\t' + dict_loci2hairpin_seq[
                        loci] + '\t' + dict_loci2strand[loci] + '\t' + dict_loci2mature_loci[loci].split('-')[
                        0] + '\t' + dict_loci2mature_loci[loci].split('-')[1] + '\t' + dict_loci2mature_seq[
                        loci] + '\t' + dict_loci2star_loci[loci].split('-')[0] + '\t' +
                    dict_loci2star_loci[loci].split('-')[1] + '\t' + dict_loci2star_seq[loci] + '\t' +
                    dict_loci2star_predict[loci] + '\t' + dict_loci2hairpin_struct[loci] + '\n')

    w1.close()
    return 'shortstack.txt'


def info(list_file_in1,file_in2,file_out,reads):
    a2 = open(file_in2, 'r')
    a3 = open('shortstack.txt', 'r')
    a4 = open('Counts.txt', 'r')
    w1 = open(file_out, 'w')
    dictl_loci2name = {}
    for i in list_file_in1:
        a1 = open(i, 'r')
        for line in a1:
            li = line.split()
            dictl_loci2name[li[0]] = li[1]
    for line in a2:
        li = line.split()
        if line.startswith('>'):
            loci = li[0][1:]
            name = li[0][1:]
            dictl_loci2name[loci] = name
    for key in dictl_loci2name:
        original_loci = key
        name = dictl_loci2name[key]
        a3.seek(0)
        for line2 in a3:
            li2 = line2.split()
            if li2[0] == original_loci:
                Chr = li2[2]
                hairpin_start = li2[3]
                hairpin_end = li2[4]
                hairpin_seq = li2[5]
                strand = li2[6]
                mature_start = li2[7]
                mature_end = li2[8]
                mature_seq_inter = li2[9]
                mature_seq = ''
                for i in mature_seq_inter:
                    if i == 'a':
                        if strand == '-':
                            i = 'U'
                        else:
                            i = 'A'
                    if i == 't':
                        if strand == '-':
                            i = 'A'
                        else:
                            i = 'U'
                    if i == 'g':
                        if strand == '-':
                            i = 'C'
                        else:
                            i = 'G'
                    if i == 'c':
                        if strand == '-':
                            i = 'G'
                        else:
                            i = 'C'
                    mature_seq = mature_seq + i
                star_start = li2[10]
                star_end = li2[11]
                star_seq = li2[13]
                struct = li2[14]
                if li2[12].startswith('X'):
                    star_detected = 'NO'
                else:
                    star_detected = 'YES'
                a4.seek(0)
        for line3 in a4:
            li3 = line3.split()
            if li3[0] == original_loci:
                counts = str(int(li3[2]) / reads * 10000000)
        w1.write(
            original_loci + '\t' + name + '\t' + Chr + '\t' + hairpin_start + '\t' + hairpin_end + '\t' + hairpin_seq + '\t' + struct + '\t' + strand + '\t' + mature_start + '\t' + mature_end + '\t' + mature_seq + '\t' + star_start + '\t' + star_end + '\t' + star_seq + '\t' + star_detected + '\t' + counts + '\n')
    w1.close()


def overlap(interval1,interval2):
    s=max(interval1[0],interval2[0])
    e=min(interval1[1],interval2[1])
    x=e-s+1
    return x


def miRNA_overlap(file_in,file_out):
    a1 = open(file_in, 'r')
    a2 = open(file_in, 'r')
    w1 = open(file_out, 'w')
    for line in a1:
        li = line.split()
        Chr = li[2]
        hairpin_start = li[3]
        hairpin_end = li[4]
        mature_start = li[8]
        mature_end = li[9]
        star_start = li[11]
        star_end = li[12]
        interval_hairpin = (int(hairpin_start), int(hairpin_end))
        interval_mature = (int(mature_start), int(mature_end))
        interval_star = (int(star_start), int(star_end))
        a2.seek(0)
        for line2 in a2:
            li2 = line2.split()
            Chr2 = li2[2]
            hairpin_start2 = li2[3]
            hairpin_end2 = li2[4]
            mature_start2 = li2[8]
            mature_end2 = li2[9]
            star_start2 = li2[11]
            star_end2 = li2[12]
            interval_hairpin2 = (int(hairpin_start2), int(hairpin_end2))
            interval_mature2 = (int(mature_start2), int(mature_end2))
            interval_star2 = (int(star_start2), int(star_end2))
            if Chr == Chr2:
                if overlap(interval_hairpin, interval_hairpin2) >= 1:
                    if overlap(interval_mature, interval_mature2) >= 1 or overlap(interval_mature, interval_star2) >= 1:
                        if line != line2:
                            w1.write(line + line2 + '\n')
    w1.close()
    return file_out


def overlap_delect(file_in1,file_in2,file_out):
    a1 = open(file_in1, 'r')
    a2 = open(file_in2, 'r')
    w1 = open(file_out, 'w')
    number = 1
    dict_name2tptm = {}
    list_remove = []
    for line in a1:
        if number == 1:
            li = line.split()
            name1 = li[1]
            dict_name2tptm[name1] = li[15]
        if number == 2:
            li = line.split()
            name2 = li[1]
            dict_name2tptm[name2] = li[15]
            if float(dict_name2tptm[name2]) < float(dict_name2tptm[name1]):
                list_remove.append(name2)
            else:
                list_remove.append(name1)
        if number == 3:
            number = 0
        number = number + 1
    for line in a2:
        li = line.split()
        if li[1] not in list_remove:
            w1.write(line)
    w1.close()


def extract_novel_T():
    a1=open('info2.txt','r')
    w1=open('novel_T.fasta','w')
    for line in a1:
        li = line.split()
        if li[1] == li[0]:
            w1.write('>' + li[1] + '\n' + li[10].replace('U', 'T') + '\n')
    w1.close()
    return 'novel_T.fasta'


def hairpin_patman_name(pre_name):
    a1 = open('info2.txt', 'r')
    a2 = open('novel_T.fasta_hit.txt', 'r')
    a3 = open('info2.txt', 'r')
    w1 = open('known2.txt', 'w')
    list_letter = []
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append(chr(i))
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append('a'+chr(i))
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append('b'+chr(i))
    patt = r'miR\d+'
    dict_family2num = {}
    for line in a1:
        li = line.split()
        if li[1] != li[0]:
            pattern = re.compile(patt)
            result = pattern.findall(li[1])
            family = result[0]
            if family not in dict_family2num:
                dict_family2num[family] = 1
            else:
                dict_family2num[family] = dict_family2num[family] + 1
    list_novel = []
    dict_loci2name = {}
    for line in a2:
        li = line.split('\t')
        pattern = re.compile(patt)
        result = pattern.findall(li[0].split()[-2])
        family = result[0]
        locus = li[1]
        list_novel.append(li[1])
        if family in dict_family2num:
            for key in dict_loci2name:
                if dict_loci2name[key].split(pre_name+'-miR')[1].endswith('.1') and dict_loci2name[key].split(pre_name+'-miR')[1][-3].isdigit():
                    if dict_loci2name[key].split(pre_name+'-miR')[1][:-2] == family[3:] :
                        dict_loci2name[key] = pre_name+'-' + family + 'a.1'
            name = pre_name+'-' + family + list_letter[dict_family2num[family]] + '.1'
            dict_loci2name[locus] = name
            dict_family2num[family] = dict_family2num[family] + 1
        else:
            name = pre_name+'-' + family + '.1'
            dict_loci2name[locus] = name
            dict_family2num[family] = 1
    for line in a3:
        li = line.split()
        if li[1] in list_novel:
            w1.write(li[0] + '\t' + dict_loci2name[li[1]] + '\n')
    w1.close()
    return 'known2.txt'


def extract_novel_star_T():
    a1=open('info4.txt','r')
    w1=open('novel_star_T.fasta','w')
    for line in a1:
        li = line.split()
        if li[1] == li[0]:
            w1.write('>' + li[1] + '\n' + li[13].replace('U', 'T') + '\n')
    w1.close()
    return 'novel_star_T.fasta'


def perl_extension(file_in):
    if '_T' in file_in:
        file_out=file_in.replace('_T','_extension')
    else:
        file_out=file_in.split('.fasta')[0]+'_extension.fasta'
    val=os.system('perl mature_extension2.pl '+file_in+' '+file_out)


def run_patman2():
    os.system('patman -e 4 -g 0 -D shortstack_Y_filter2_extension.fasta -P shortstack_Y_filter2_T.fasta --singlestrand -o shortstack_Y_filter2.fasta_out.txt')
    os.system('patman -e 4 -g 0 -D novel_star_extension.fasta -P shortstack_Y_filter2_T.fasta --singlestrand -o shortstack_Y_filter2-star.fasta_out.txt')


def novel_Cluster():
    a1 = open('shortstack_Y_filter2.fasta_out.txt', 'r')
    a2 = open('shortstack_Y_filter2-star.fasta_out.txt', 'r')
    a3 = open('shortstack_Y_filter2.fasta', 'r')
    w1 = open('novel_cluster.txt', 'w')
    dict_sub2query = {}
    for line in a1:
        li = line.split()
        if abs(int(li[2]) - 26) + int(li[5]) <= 4:
            sub = li[0]
            if sub not in dict_sub2query:
                dict_sub2query[sub] = li[1]
            else:
                dict_sub2query[sub] = dict_sub2query[sub] + '\t' + li[1]
    for line in a2:
        li = line.split()
        if abs(int(li[2]) - 26) + int(li[5]) <= 4:
            sub = li[0]
            if sub not in dict_sub2query:
                dict_sub2query[sub] = li[1]
            else:
                dict_sub2query[sub] = dict_sub2query[sub] + '\t' + li[1]
    for line in a3:
        if line.startswith('>'):
            li = line.split()
            if li[0][1:] not in dict_sub2query:
                dict_sub2query[li[0][1:]] = li[0][1:]
    for key in dict_sub2query:
        w1.write(key + '\t' + dict_sub2query[key] + '\n')
    w1.close()


def sort():
    a1 = open('novel_cluster.txt', 'r')
    w1 = open('novel_cluster_sorted.txt', 'w')
    patt = r'\d+'
    pattern = re.compile(patt)
    w1.write(''.join(sorted(a1, key=lambda s: (int(pattern.findall(s.split()[0].split(':')[0])[0]), int(s.split()[0].split(':')[1].split('-')[0])))))
    w1.close()


def novel_family():
    a1 = open('novel_cluster_sorted.txt', 'r')
    w1 = open('novel_cluster_family.txt', 'w')
    number = 1
    dict_name2family = {}
    for line in a1:
        li = line.split()
        list_pass = ''
        for i in li:
            if i in dict_name2family:
                list_pass = 'IN'
                fimily_known = dict_name2family[i]
        if list_pass == 'IN':
            for a in li:
                dict_name2family[a] = fimily_known
        if list_pass != 'IN':
            for a in li:
                dict_name2family[a] = 'miRN' + str(number)
            number = number + 1
    key_number = 1
    while key_number < number:
        for o in dict_name2family:
            if dict_name2family[o] == 'miRN' + str(key_number):
                w1.write(o + '\t')
        w1.write('\n')
        key_number = key_number + 1
    w1.close()


def novel_delect():
    global number
    global dict_name2family
    number = 1
    dict_name2family = {}
    with open('novel_cluster_family.txt') as a1:
        for line in a1:
            li=line.split()
            list_pass = ''
            for i in li:
                i=i.split('-star')[0]
                if i in dict_name2family:
                    list_pass = 'IN'
                    fimily_known = dict_name2family[i]
            if list_pass == 'IN':
                for a in li:
                    a = a.split('-star')[0]
                    dict_name2family[a] = fimily_known
            if list_pass != 'IN':
                for a in li:
                    a=a.split('-star')[0]
                    dict_name2family[a] = 'miRN' + str(number)
                number = number + 1


def novel_write():
    w1=open('novel_cluster_family2.txt','w')
    key_number = 1
    while key_number < number:
        for o in dict_name2family:
            if dict_name2family[o] == 'miRN' + str(key_number):
                w1.write(o + '\t')
        w1.write('\n')
        key_number = key_number + 1
    w1.close()


def novel_family_number():
    w1=open('novel_cluster_family3.txt','w')
    with open('novel_cluster_family2.txt','r') as a1:
        number=1
        for line in a1:
            w1.write('miRN'+str(number)+'\t'+line)
            number=number+1
    w1.close()


def novel_name(pre_name):
    a1 = open('novel_cluster_family3.txt', 'r')
    w1 = open('novel_cluster_name.txt', 'w')
    list_letter = []
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append(chr(i))
    for i in range(ord('a'), ord('a') + 26):
        list_letter.append('a'+chr(i))
    dict_loci2name = {}
    for line in a1:
        li = line.split()
        if len(li) == 2:
            for i in li:
                if i.startswith('miR'):
                    family = i
                else:
                    name = pre_name+'-' + family
                    loci = i
                    w1.write(loci + '\t' + name + '\n')
        else:
            number = 0
            for i in li:
                if i.startswith('miR'):
                    family = i
                else:
                    name = pre_name+'-' + family + list_letter[number]
                    loci = i
                    w1.write(loci + '\t' + name + '\n')
                    number = number + 1
    w1.close()


def info2(pre_name):
    a0 = open('novel_cluster_name.txt', 'r')
    a1 = open('info4.txt', 'r')
    w1 = open('info5.txt', 'w')
    dict_loci2name = {}
    for line in a0:
        li = line.split()
        dict_loci2name[li[0]] = li[1]
    for line in a1:
        li = line.split()
        counts = 0
        for i in li:
            if counts != 1 and counts != 15:
                w1.write(i + '\t')
            if counts == 1:
                if i.startswith(pre_name) == False:
                    w1.write(dict_loci2name[i] + '\t')
                else:
                    w1.write(i + '\t')
            if counts == 15:
                w1.write(i + '\n')
            counts = counts + 1
    w1.close()


def count(string):
    count_left=0
    count_right=0
    for i in string:
        if i =='(':
            count_left=count_left+1
        if i ==')':
            count_right=count_right+1
    if count_left>count_right:
        left_or_right='('
    if count_right>count_left:
        left_or_right=')'
    return left_or_right


def add_3p5p():
    a1 = open('info5.txt', 'r')
    w1 = open('info6.txt', 'w')
    for line in a1:
        li = line.split()
        name = li[1]
        mature_start = li[8]
        mature_end = li[9]
        strand = li[7]
        hairpin_start = li[3]
        hairpin_end = li[4]
        struct = li[6]
        if strand == '+':
            string = struct[int(mature_start) - int(hairpin_start):int(mature_end) - int(hairpin_start)]
            if count(string) == '(':
                five_or_three = '5'
            if count(string) == ')':
                five_or_three = '3'
        if strand == '-':
            string = struct[::-1]
            string = string[int(mature_start) - int(hairpin_start):int(mature_end) - int(hairpin_start)]
            if count(string) == '(':
                five_or_three = '5'
            if count(string) == ')':
                five_or_three = '3'
        w1.write(line[:-1] + '\t' + five_or_three + '\n')
    w1.close()


def run_master(pre_name,reads):
    reads=int(reads)
    ###1
    #name_known(run_perl(run_patman(RNA2DNA(extract_Results('Results.txt',['N15','Y'])))),pre_name)
    delect_known_loci(name_known(run_perl(run_patman(RNA2DNA(extract_Results(['N15','Y'])),'mature_extension.fa')),pre_name),extract_Results(['Y']),'shortstack_Y_filter.fasta')

    ###2
    Cluster()
    info(['known.txt'],'shortstack_Y_filter.fasta','info1.txt',reads)
    overlap_delect(miRNA_overlap('info1.txt','overlap'),'info1.txt','info2.txt')

    ###3
    run_perl2(run_patman(extract_novel_T(),'hairpin_T.fa'))
    delect_known_loci(hairpin_patman_name(pre_name),'shortstack_Y_filter.fasta','shortstack_Y_filter2.fasta')
    info(['known.txt','known2.txt'],'shortstack_Y_filter2.fasta','info3.txt',reads)
    overlap_delect(miRNA_overlap('info3.txt','overlap2'),'info3.txt','info4.txt')

    ###4
    perl_extension(extract_novel_star_T())
    perl_extension('shortstack_Y_filter2.fasta')
    RNA2DNA('shortstack_Y_filter2.fasta')
    run_patman2()
    novel_Cluster()
    sort()
    novel_family()
    novel_delect()
    novel_write()
    novel_family_number()
    novel_name(pre_name)
    info2(pre_name)
    add_3p5p()

script,pre_name,reads=argv
run_master(pre_name,reads)







