# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/10 18:40

"""
import sys
import os
import platform

def reversecomp_seq(input_str):
    """
    Reverse complement of sequence
    :param input_str: str
    :return: str
    """

    trantab1 = str.maketrans('ACGTRYMKVBHDSWNacgtrymkvbhdswn',
                             'TGCAYRKMBVDHWSNtgcayrkmbvdhwsn')

    trantab2 = str.maketrans('ACGURYMKVBHDSWNacgurymkvbhdswn',
                             'UGCAYRKMBVDHWSNugcayrkmbvdhwsn')

    line = input_str.strip()

    line = line[::-1]  # reverse sequence

    if line.find("T") != -1 or line.find("t") != -1:
        line = line.translate(trantab1)
    else:
        line = line.translate(trantab2)

    return line



def make_dir(input_dir):
    """
    :param input_dir: directory which need to create
    :return:
    """
    if os.path.isdir(input_dir) == False:
        os.makedirs(input_dir)


def get_all_path(open_dir_path):
    """
    Gets all the path of files in a folder
    :param open_dir_path: folder
    :return: path of files
    """

    rootdir = open_dir_path

    path_list = []

    lists = os.listdir(rootdir)

    for i in range(0, len(lists)):
        com_path = os.path.join(rootdir, lists[i])

        if platform.system().lower() == "windows":
            com_path = com_path.replace("\\", "/")
        # print(com_path)

        if os.path.isfile(com_path):
            path_list.append(com_path)


    return path_list



def read_gff3(input_path):
    """
    :param input: gff3 file
    :return: seqs and annotation dic
    """
    gff3_file = open(input_path,"r",encoding="utf-8")
    line_num = 0
    seq_nums = []
    i = 0
    seqs = {}  # sequence name

    annotations = {}

    for line in gff3_file:
        line_num += 1
        line = line.strip()
        
        if line_num >= 2 and line != "":
            if line.startswith("##seq"):
                seq_num = line.split(" ")[3]
                # print(str(seq_num))
                seq_nums.append(seq_num)

            elif line.startswith("###"):
                pass

            elif line.startswith("#"):
                seqs[seq_nums[i]] = line[1:]
                i += 1

            else:
                anno_list = line.split("\t")

                seq_name = anno_list[0]

                if seq_name not in annotations:
                    annotations[seq_name] = [line]

                else:
                    annotations[seq_name].append(line)



    gff3_file.close()

    return seqs, annotations



def best_domain(input_list,n):
    """
    After sorting the nested list, a list with the smallest element in it
    :param input_list: Enter a nested list
    :param n: number of column used for sort
    :return: the best domain
    """
    hit_domain = "not found"

    if len(input_list) >= 1:
        hit_list_sort = sorted(input_list, key=lambda x: x[n])

        hit_list = hit_list_sort[0]

        hit_domain = (hit_list[0] + "," + str(hit_list[1]) + "," + hit_list[2]
                            + "," + str(hit_list[3]) + "," + str(hit_list[4]))
        # name of domain, c-Evalue, ORF, Start and end point, lenth of domain


    return hit_domain


def domain_locate_genome(seq_coordinates, seq_stand, domain_aa_range,orf):
    """
    domain relocation on the original genome

    :param seq_coordinates: The coordinates of the sequence on the original
    genome, for example: "2000..4902"

    :param seq_stand: The direction of the sequence, "+" means the same as the
    input genome, "-" means the complementary strand

    :param domain_coordinates: The coordinates (amino acids) of the domain
    predicted by HMMER in the sequence, the format is for example: "100..202"

    :param orf: When the domain is predicted, which set of coding frames of
    the sequence is used, such as "ORF1"

    :return: The location on the original genome (for example: "1900..2202")
    and the chain where the domain is located
    """
    seq_coordinates_start = int(seq_coordinates.split("..")[0])
    seq_coordinates_end = int(seq_coordinates.split("..")[1])

    seq_lenth = seq_coordinates_end - seq_coordinates_start + 1

    domain_aa_range = (domain_aa_range.split("(")[-1]).split(")")[0]

    domain_in_strand = "ERROR"
    site_range = "(..)"


    if seq_stand == "same":
        if orf == "ORF1":
            domain_nt_start = 3 *(int(domain_aa_range.split("..")[0]) - 1)  + 1
            domain_nt_end = 3 * int(domain_aa_range.split("..")[1])

            genome_start = seq_coordinates_start + domain_nt_start - 1
            genome_end = seq_coordinates_start + domain_nt_end - 1

            domain_in_strand = "same"
            site_range = "(" + str(genome_start) + ".." + str(genome_end) + ")"

        elif orf == "ORF2":
            domain_nt_start = 3 * (int(domain_aa_range.split("..")[0]) - 1) + 1 + 1
            domain_nt_end = 3 * int(domain_aa_range.split("..")[1]) + 1
            genome_start = seq_coordinates_start + domain_nt_start - 1
            genome_end = seq_coordinates_start + domain_nt_end - 1

            domain_in_strand = "same"
            site_range = "(" + str(genome_start) + ".." + str(genome_end) + ")"

        elif orf == "ORF3":
            domain_nt_start = 3 * (
                        int(domain_aa_range.split("..")[0]) - 1) + 1 + 2
            domain_nt_end = 3 * int(domain_aa_range.split("..")[1]) + 2
            genome_start = seq_coordinates_start + domain_nt_start - 1
            genome_end = seq_coordinates_start + domain_nt_end - 1

            domain_in_strand = "same"
            site_range = "(" + str(genome_start) + ".." + str(genome_end) + ")"


    elif seq_stand == "complement":
        if orf == "ORF1":
            domain_nt_start = 3 * (
                    int(domain_aa_range.split("..")[0]) - 1) + 1
            domain_nt_end = 3 * int(domain_aa_range.split("..")[1])

            domain_nt_restart = seq_lenth - domain_nt_end + 1
            domain_nt_reend = seq_lenth - domain_nt_start + 1

            genome_start = seq_coordinates_start + domain_nt_restart - 1
            genome_end = seq_coordinates_start + domain_nt_reend - 1

            domain_in_strand = "complement"
            site_range = ("complement(" + str(genome_start)
                          + ".." + str(genome_end) + ")")

        elif orf == "ORF2":
            domain_nt_start = 3 * (
                    int(domain_aa_range.split("..")[0]) - 1) + 1 + 1
            domain_nt_end = 3 * int(domain_aa_range.split("..")[1]) + 1

            domain_nt_restart = seq_lenth - domain_nt_end + 1
            domain_nt_reend = seq_lenth - domain_nt_start + 1

            genome_start = seq_coordinates_start + domain_nt_restart - 1
            genome_end = seq_coordinates_start + domain_nt_reend - 1

            domain_in_strand = "complement"
            site_range = ("complement(" + str(genome_start)
                         + ".." + str(genome_end) + ")")

        elif orf == "ORF3":
            domain_nt_start = 3 * (
                    int(domain_aa_range.split("..")[0]) - 1) + 1 + 2
            domain_nt_end = 3 * int(domain_aa_range.split("..")[1]) + 2

            domain_nt_restart = seq_lenth - domain_nt_end + 1
            domain_nt_reend = seq_lenth - domain_nt_start + 1

            genome_start = seq_coordinates_start + domain_nt_restart - 1
            genome_end = seq_coordinates_start + domain_nt_reend - 1

            domain_in_strand = "complement"
            site_range = ("complement(" + str(genome_start)
                          + ".." + str(genome_end) + ")")



    return (domain_in_strand,site_range)



def hmmer_contain(input_path):

    domain = []
    line_num = 0
    with open(input_path,"r", encoding="utf-8") as input_file:
        for line in input_file:
            line_num += 1
            if line_num <= 5:
                if not line.startswith("#"):
                    domain.append(line.strip())

    if domain == []:
        return "False"

    else:
        return "True"


