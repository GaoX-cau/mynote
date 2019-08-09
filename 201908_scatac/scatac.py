#!/usr/bin/env python3


import sys
import os
import time
import collections
import pysam
import numpy as np


class fragment(object):
    """A fragment object that has the following attributes:
    Attributes:
        chrom: chromsome name
        start: start position
        end: end position
        mapq: mapping quality
        is_proper_pair: whether properly paired
        is_single: whether it is a single read
        is_secondary: whether it is a secondary alignment read
        flen: fragment length
    """
    def __init__(self, qname, chrom, pos, flen, mapq, is_single, is_secondary, is_proper_pair):
        """Return a qc object"""
        self.qname = qname
        self.chrom = chrom
        self.pos = pos
        self.flen = flen
        self.mapq = mapq
        self.is_single = is_single
        self.is_secondary = is_secondary
        self.is_proper_pair = is_proper_pair


class qc(object):
    """A quality control object that has the following attributes:
    Attributes:
        total: total number of sequenced fragments.
        mapped: number of mappable fragments.
        chrM: number of fragments mapped to chrM.
        paired: number of fragments are paired.
        single: number of fragments are from single read.
        proper_paired: number of paired reads are properly paired.
        usable: number of usable fragments.
        uniq: number of unique fragments.
        isize: average insert size distribution.
    """
    def __init__(self):
        """Return a qc object"""
        self.id = 0
        self.total = 0
        self.mapped = 0
        self.single = 0
        self.secondary = 0
        self.paired = 0
        self.proper_paired = 0
        self.proper_flen = 0
        self.usable = 0
        self.uniq = 0
        self.chrM = 0
        self.final = 0


def readGenomeSizeFromTxt(fname):
    """
    Read genome information.
    
    Args:
    -----
        fname: a txt file contains genome information
    Returns:
    -----
        A dictionary contains SQ as key and SL as value, otherwise None
    """
    # first check if fname is a bam file
    res = dict()
    with open(fname) as fin:
        for line in fin:
            elems = line.split()
            chr_name = elems[0]
            chr_len = int(elems[1])
            res[chr_name] = chr_len
    return res


def getBinsFromGenomeSize(genome_dict, bin_size):
    """Create a dictionary contains all bins of the same size across the genome
    
    Attributes:
        binSize: bin size (i.e. 5000)
    
        genomeDict: a dictionary contains chromosome sizes
    
    Return:
        A dictionary contains all bins and its index (start from 1)
    """
    bin_dict = collections.OrderedDict()
    i = 1
    for _chrom in genome_dict:
    	for _start in range(1, genome_dict[_chrom], bin_size):
            _end = min(_start + bin_size - 1, genome_dict[_chrom])
            _binId = (_chrom , _start, _end)
            bin_dict[_binId] = i  ## set bin index
            i += 1
    return bin_dict


def group_reads_by_barcode(input_bam):
    """ Group reads based on the barcodes
    
    Args:
        input_bam: a bam file
    Returns:
        Generator that contains reads sharing the same barcode
    """
    if not os.path.exists(input_bam): 
        print(("Error @group_reads_by_barcode: " + input_bam + " does not exist!"));
 
    read_group_list = []; 
    pre_barcode = ""
    samfile = pysam.AlignmentFile(input_bam, "rb")
    for cur_read in samfile:
        cur_barcode = cur_read.qname.split(":")[0]
        if cur_barcode == pre_barcode:
            read_group_list.append(cur_read)
        else:
            if pre_barcode != "":
                # return read group
                yield (x for x in read_group_list)
            read_group_list = [cur_read] # add the first read
            pre_barcode = cur_barcode
    # reads from the last barcode
    yield (x for x in read_group_list)


def pairReadsByName(read_list):
    """ Pair reads based on read names
    
    Args:
        read_list: a list of reads that share the same barcode
    Returns:
        Generator contains read pairs from the same fragment
        and a boolen variable indicates whether it is supplementary alignment
    """
    # pair up 
    for read1 in read_list:
        # read until read1 is not a supplementary alignment
        while(read1.is_supplementary):
            yield (read1, None, False, True)
            try:
                #print "Warning: skipping ", read1.qname;
                read1 = next(read_list)
            except:
            	break        
        try:
            read2 = next(read_list)
        except:
        	break
        while(read2.is_supplementary):
            yield (read2, None, False, True)
            try:
                #print "Warning: skipping ", read2.qname;
                read2 = next(read_list)
            except:
            	break
        if(read1.qname != read2.qname):
            while (read1.qname != read2.qname):
                yield(read1, None, False, False)
                read1 = read2
                try:
                    read2 = next(read_list)
                    while(read2.is_supplementary):
                        try:
                            #print "Warning: skipping ", read2.qname;
                            read2 = next(read_list)
                        except:
                        	break
                except:
                    break
        yield (read1, read2, True, False)


def readPairToFragment(read1, read2, is_secondary):
    """ convert read pairs to fragments
    
    Args:
        read1: R1 read
        read2: R2 read
    Returns:
        Generator that contains a fragment object
    """
    try:
        read1.qname
        read2.qname
        if read1.qname != read2.qname:
            sys.exit('read_pair_to_fragment: read1 and read2 name does not match!');  
    except ValueError as e:
        sys.exit('read_pair_to_fragment: can not get read1 or read2 name!');   
    barcode = read1.qname.split(":")[0]
    mapq = min(read1.mapq, read2.mapq)
    try:
        chrom1 = read1.reference_name
        start1 = read1.reference_start
        strand1 = "-" if read1.is_reverse else "+"
        chrom2 = read2.reference_name
        start2 = read2.reference_start
        strand2 = "-" if read1.is_reverse else "+"
        # it is possible that flen1 is None  
        flen1 = read1.reference_length if read1.reference_length != None else 0
        flen2 = read2.reference_length if read2.reference_length != None else 0
        end1   = start1 + flen1
        end2   = start2 + flen2
        start = min(start1, end1, start2, end2)
        end = max(start1, end1, start2, end2)
        return fragment(read1.qname, chrom1, start, abs(start - end), mapq, False, is_secondary, read1.is_proper_pair)
    except ValueError as e:
        return fragment(read1.qname, None, None, None, mapq, False, is_secondary, False)


def isBamQuerynameSorted(fname):
    """Check if a given bam file is sorted based on the queryname.
    
    Args:
        fname: bam file name
    Returns:
        True if fname is a bam file sorted by read name, otherwise False
    """
    samfile = pysam.AlignmentFile(fname, "rb");
    flag = False
    header = samfile.header
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "queryname"):
                flag = True
    samfile.close()
    return flag


def getBarcodesFromBam(input_bam):
    """Identify unique barcodes from the bam file
    
    Args:
        input_bam: a bam file
    Returns:
        A dictionary contains all barcodes, otherwise None
    """
    
    barcode_dict = collections.OrderedDict()
    samfile = pysam.AlignmentFile(input_bam, "rb")
    i = 1
    for _read in samfile:
        barcode = _read.qname.split(":")[0].upper()
        if barcode not in barcode_dict:
            barcode_dict[barcode] = i
            i += 1
    samfile.close()
    return barcode_dict


def readToFragment(read1, is_secondary):
    """ convert a single read to fragment
    
    Args:
        read1: a single read
    Returns:
        Generator contains a fragment object
    """
    try:
        read1.qname
    except ValueError as e:
        sys.exit('readto_fragment: can not get read1 name!');   

    barcode = read1.qname.split(":")[0]
    mapq = read1.mapq
    try:
        chrom1 = read1.reference_name
        start1 = read1.reference_start
        strand1 = "-" if read1.is_reverse else "+"   
        # it is possible that flen1 is None  
        flen1 = read1.reference_length if read1.reference_length != None else 0
        end1   = start1 + flen1
        start = min(start1, end1)
        end = max(start1, end1)
        #(qname, chrom, pos, flen, mapq, is_single, is_secondary, is_proper_pair):
        return fragment(read1.qname, chrom1, start, abs(start - end), mapq, True, is_secondary, False)
    except ValueError as e:
        return fragment(read1.qname, None, None, None, mapq, True, is_secondary, False)


def getBinsFromGenomeSize(genome_dict, bin_size):
    """Create a dictionary contains all bins of the same size across the genome
    
    Attributes:
        binSize: bin size (i.e. 5000)
    
        genomeDict: a dictionary contains chromosome sizes
    
    Return:
        A dictionary contains all bins and its index (start from 1)
    """
    bin_dict = collections.OrderedDict()
    i = 1
    for _chrom in genome_dict:
    	for _start in range(1, genome_dict[_chrom], bin_size):
            _end = min(_start + bin_size - 1, genome_dict[_chrom])
            _binId = (_chrom , _start, _end)
            bin_dict[_binId] = i
            i = i + 1
    return bin_dict


def preProcess(input_file, genome_name, bin_size=5000, min_mapq=10, min_flen=0, max_flen=1000):
    """
    Pre-processing to create a snap file from a bam or bed that contains alignments or a bed file that contains fragments.
    Args:
    --------
    input_file: 
        a bam file contains fragments 
        
    output: 
        idx
        idy
    
    Optional
    --------
    min_mapq: 
            min mappability [10]. fragments with mappability less than 10 will be filtered
               
    min_flen: 
            min fragment size [0]. fragments of length shorter than min_flen will be filtered
    max_flen: 
            max fragment size [1000]. fragments of length bigger than min_flen will be filtered
    min_cov:
            min coverage per barcode. barcodes with sequencing fragments less than min_cov will be filtered before processed  
    max_num:
            max number of barcodes to be stored. Based on the coverage, top max_barcode barcodes are selected and stored.
    
    barcode_file: 
            a txt file that contains selected barcodes to create count matrix [None]
    keep_chrm:
            a boolen variable indicates whether to keep reads mapped to chrM [True]
    keep_single:
            a boolen variable indicates whether to keep single-end reads [False]
    keep_secondary:
            a boolen variable indicates whether to keep secondary alingments [False]
    keep_discordant:
            a boolen variable indicates whether to keep discordant alingments [False]
            
    tmp_folder: 
            where to store the temporary files [None];
    
    qc_file: 
            a boolen variable indicates whether to create a master qc file [True];
    verbose: 
            a boolen variable indicates whether to output the progress [True];
    """
    if not os.path.exists(input_file):
        print(('error: ' + input_file + ' does not exist!'))
        sys.exit(1)
    
    if not isBamQuerynameSorted(input_file):
        print(f'error: {input_file} is not a sorted by read name!')
        sys.exit(1)

    barcode_dict = getBarcodesFromBam(input_file)
    genome_dict = readGenomeSizeFromTxt(genome_name)
    bin_dict = getBinsFromGenomeSize(genome_dict, bin_size)

    check_point = 0
    start_time = time.time()
    qc_dict = collections.defaultdict(qc)
    num_barcode = len(barcode_dict)

    idxList = []  # barcode index list
    idyList = []  # bin index list
    countList = []  # number of count

    for read_group in group_reads_by_barcode(input_file):
        frag_list = []
        for (read1, read2, is_paired, is_secondary) in pairReadsByName(read_group):   
            if is_paired:
                frag = readPairToFragment(read1, read2, is_secondary)
            else: # supplementary alignments or unmated reads;
                frag = readToFragment(read1, is_secondary)
            # extract the barcode
            barcode = frag.qname.split(":")[0].upper()
            # only for printing the progress
            check_point += 1
            if check_point%100000 == 0:
                print(f'{check_point} tags, {time.time() - start_time:.2f} seconds')
                #print(("%d\ttags, %s seconds " % (check_point, time.time() - start_time)));  
            # skip this barcode if it is not in the barcode list 
            if barcode not in barcode_dict: 
                break
            # total number of sequencing fragments (exclude supplementary alignments)
            if frag.is_secondary == False:
                qc_dict[barcode].total += 1

            ## 1. Filter non-uniquely mapped fragments
            if frag.mapq < min_mapq:
                continue;           
            qc_dict[barcode].mapped += 1

            # 2. if it is single read, keep it only if keep_unpaired is true
            #    if it is paired, keep it only if it is approperly paired
            if frag.is_single == True: 
                if frag.is_secondary:
                    qc_dict[barcode].secondary += 1
                    if not keep_secondary:
                        continue
                else:
                    qc_dict[barcode].single += 1
                    if not keep_single:
                        continue                    
            else: # paired-end reads
                qc_dict[barcode].paired += 1
                if frag.is_proper_pair:
                    qc_dict[barcode].proper_paired += 1
                else:
                    if not keep_discordant:
                        continue    
            # 3. check fragment size
            if frag.flen > min_flen and frag.flen < max_flen:
                qc_dict[barcode].proper_flen += 1
            else:
                continue
            # 4. combine single and paired as fragments
            frag_list.append((frag.chrom, frag.pos, frag.pos+frag.flen, barcode))
        
        # 5. remove duplicate fragments
        qc_dict[barcode].usable = len(frag_list)
        frag_list_uniq = set(frag_list) # remove duplicated fragments
        qc_dict[barcode].uniq = len(frag_list_uniq)

        # 6. cell by bins

        bins = collections.defaultdict(lambda : 0)
        for item in frag_list_uniq:
            bin_chr = item[0]
            for bin_pos in set([int(item[1]/bin_size) * bin_size + 1, int(item[2]/bin_size) * bin_size + 1]):
                bins[(bin_chr, bin_pos, bin_pos + bin_size - 1)] += 1
    
        for key in bins:
            if key in bin_dict:
                idyList.append(bin_dict[key])
                countList.append(bins[key])
                idxList.append(barcode_dict[barcode])

        del frag_list, frag_list_uniq
    
    return idxList, idyList, countList, qc_dict


if __name__ == '__main__':

    input_file = sys.argv[1]
    genome_name = sys.argv[2]
    idxList, idyList, countList, qc_dict = preProcess(input_file, genome_name)

    barIdx = dict(list(zip(list(set(idyList)),range(len(set(idyList))))))

    x = np.zeros((len(set(idxList)),len(set(idyList))), dtype='uint8')

    for i in range(len(idyList)):
        x[idxList[i]-1][barIdx[idyList[i]]] = 1


    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca.fit(x)
    X_pca = pca.transform(x)
    plt.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.2, linewidths=0)
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    from sklearn.manifold import TSNE
    X_tsne = TSNE(n_components=2).fit_transform(X)
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], alpha=0.2, linewidths=0)