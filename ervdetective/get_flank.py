# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/10 14:25

"""
from Bio import SeqIO
from public_function import reversecomp_seq,get_all_path

class GetFlanks(object):

    def __init__(self, host_seq, blast_out_dir,
                 flanks_lenth, out_path):
        """
        Get flanks sequence near the blast hit sites
        :param host_seq: path of input genome or sequence
        :param blast_out_dir: path of blast_result
        :param flanks_lenth: the lenth of flanks_lenth
        :param out_path: path flanks sequence
        """
        super(GetFlanks, self).__init__()

        self.host_path = host_seq
        self.blast_out_dir = blast_out_dir # Support for multiple probes
        self.out_path = out_path
        self.flank_length= flanks_lenth


    def run(self):
        blast_file_list = get_all_path(self.blast_out_dir)

        EVEs = {}  # to combine the blast results of protein probes

        for each_blast_path in blast_file_list:
            matchlenth_threshold = 0
            each_blast_path = each_blast_path.strip()

            if each_blast_path.upper().find("RT") != -1:
                matchlenth_threshold = 50   # alignment length of aa

            elif each_blast_path.upper().find("GAG") != -1:
                matchlenth_threshold = 150  # alignment length of aa

            elif each_blast_path.upper().find("ENV") != -1:
                matchlenth_threshold = 200  # alignment length of aa


            blast_file = open(each_blast_path,"r",encoding="utf-8")

            for line in blast_file:

                blast_line_list = line.strip().split("\t")

                if int(blast_line_list[3]) > matchlenth_threshold:
                    if int(blast_line_list[-4]) > int(blast_line_list[-3]):

                        # On the complementary strand of the input genome sequence
                        strand = "complement"

                        """
                        the starting position is relabeled according to the 
                        numerical value from small to large
                        """
                        start = int(blast_line_list[-3])
                        end = int(blast_line_list[-4])
                    else:
                        strand = "same"  # On the same chain of the input sequence
                        start = int(blast_line_list[-4])
                        end = int(blast_line_list[-3])

                    new_line_list = [start, end, strand,
                                     float(blast_line_list[2]),
                                     int(blast_line_list[3]),
                                     float(blast_line_list[-2])]

                    # start, end, strand, pident, lenth, e-value

                    if not EVEs.__contains__(blast_line_list[1]):
                        EVEs[blast_line_list[1]] = []

                    # Add RT protein comment
                    EVEs[blast_line_list[1]].append(new_line_list)

            blast_file.close()

        ref_EVEs = {}    # processing hit event of continuous blast sites

        for subject_name in EVEs.keys():

            hit_list = EVEs[subject_name]

            # Sort the nested list according to the starting position
            hit_list_sort = sorted(hit_list,
                                   key=lambda x: x[0])

            # Sort the nested list according to the direction of the chain
            hit_list_sort = sorted(hit_list_sort,
                                   key=lambda x: x[2])
            # print(subject_name,hit_list_sort)

            # the first hit in the list as the starting value
            last_start = hit_list_sort[0][0]
            last_end = hit_list_sort[0][1]
            last_strand = hit_list_sort[0][2]

            fragment_distance = 10000  # distance threshold for blast consecutive hit events

            if not ref_EVEs.__contains__(subject_name):
                ref_EVEs[subject_name] = []

            """
            A more ingenious method is used here. After exiting the continuous 
            hit area, the range of the last continuous hit area is recorded. 
            """
            for each_hit in hit_list_sort:

                # If two blast hits on the same chain
                if each_hit[2] == last_strand:

                    # Judging the starting point distance
                    if each_hit[0] <= last_start + fragment_distance:

                        # Take the smaller value as the starting point
                        last_start = last_start

                        # Take the larger value as the end point
                        last_end = max(last_end, each_hit[1])

                    # if out of range
                    else:
                        # Record the result of the previous
                        ref_EVEs[subject_name].append(
                            [last_start, last_end, last_strand])

                        last_start = each_hit[0]
                        last_end = each_hit[1]
                        last_strand = each_hit[2]

                # If two blast hits on different chains
                else:
                    # Record the last result
                    ref_EVEs[subject_name].append(
                        [last_start, last_end, last_strand])

                    last_start = each_hit[0]
                    last_end = each_hit[1]
                    last_strand = each_hit[2]

            # Supplement the last result
            ref_EVEs[subject_name].append([last_start, last_end, last_strand])

        host_seq_dict = SeqIO.to_dict(SeqIO.parse(self.host_path, "fasta"))
        """
        Note: the key value in host_seq_dict only contains the part before
        the space in the sequence name 
        """

        # self.flank_length = 15000
        seqs = {}
        sum_lsit = []

        for subject_name in ref_EVEs.keys():
            for a in ref_EVEs[subject_name]:
                sub_id = subject_name
                s_start = a[0]
                s_end = a[1]
                strand = a[2]

                start = int(s_start) - int(self.flank_length)
                end = int(s_end) + int(self.flank_length)

                new_seq = ""

                if start <= 1 and end <= len(host_seq_dict[sub_id].seq[:]):
                    start = 1
                    end = end
                    new_seq = host_seq_dict[sub_id].seq[start - 1:end]

                elif start >= 1 and end <= len(host_seq_dict[sub_id].seq[:]):
                    start = start
                    end = end
                    new_seq = host_seq_dict[sub_id].seq[start - 1: end]

                elif start >= 1 and end >= len(host_seq_dict[sub_id].seq[:]):
                    start = start
                    end = len(host_seq_dict[sub_id].seq[:])
                    new_seq = host_seq_dict[sub_id].seq[start - 1: end]

                elif start <= 1 and end >= len(host_seq_dict[sub_id].seq[:]):
                    start = 1
                    end = len(host_seq_dict[sub_id].seq[:])
                    new_seq = host_seq_dict[sub_id].seq[start - 1:end]

                else:
                    print("error")

                name = (">" + sub_id + "|" + str(start) + "|" + str(end) + "|"
                        + strand)
                seq = (str(new_seq))

                seqs[name] = seq

                sum_lsit.append(name)

        #  Remove areas with full containment relations (leave wide areas) in same chain
        use_list = []

        for line1 in sum_lsit:
            qury_name, qury_start, qury_end, qury_direction = line1.split("|",
                                                                          3)
            for line2 in sum_lsit:

                if line1 != line2:

                    refer_name, refer_start, refer_end, refer_direction = line2.split(
                        "|", 3)
                    if (qury_name == refer_name and
                            qury_direction == refer_direction):
                        if qury_start >= refer_start and qury_end <= refer_end:
                            break

            else:
                use_list.append(line1)



        ganjing = set()

        finnaly_use = []

        for each_seqname in use_list:
            each_seqname = each_seqname.strip()

            if not each_seqname in ganjing:
                finnaly_use.append(each_seqname)
                ganjing.add(each_seqname)

        # get flanking sequence
        with open(self.out_path, "w",encoding="utf-8") as erv_flanking:
            for line in finnaly_use:

                if line.split("|")[-1] == "complement":
                    erv_flanking.write(line + "\n"
                                        + reversecomp_seq(
                        seqs[line].strip()) + "\n")   #  Reverse complement the sequence on the negative strand
                else:
                    erv_flanking.write(line + "\n" + seqs[line].strip() + "\n")


