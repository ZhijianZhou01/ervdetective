# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/10 20:46

"""

import os

from public_function import read_gff3

class SimplyGFF(object):

    def __init__(self, input_gff3,
                 simplify_gff3,
                 duplicate_removal_gff3,
                 mintsd):
        super(SimplyGFF, self).__init__()
        
        self.gff3_path = input_gff3
        self.simplify_gff3_path = simplify_gff3
        self.duplicate_removal_path = duplicate_removal_gff3
        self.mintsd = int(mintsd)

    def run(self):

        record_dict = {}

        seqs, annotations = read_gff3(self.gff3_path)


        for key in annotations.keys():
            annos = annotations[key]

            name = ""

            for line in annos:
                anno_list = line.strip().split("\t")

                if anno_list[2] == "repeat_region":
                    region = anno_list[3] + ".." + anno_list[4]
                    name = anno_list[0] + "|" + region
                    if not record_dict.__contains__(name):
                        record_dict[name] = []

                elif anno_list[2] == "LTR_retrotransposon":
                    similarity = anno_list[-1].split(";")[-2].replace("ltr_similarity=", "")
                    record_dict[name].append(similarity)

                elif anno_list[2] == "long_terminal_repeat":
                    ltr_region = anno_list[3] + ".." + anno_list[4]
                    record_dict[name].append(ltr_region)


        # print(record_dict)

        out_simplify_file = open(self.simplify_gff3_path, "w",encoding="utf-8")
        out_simplify_file.write(
            "#name of flank seq" + "\t" + "repeat_region in flank seq"
            + "\t" + "similarity of pairLTRs" + "\t" + "LTR1 in flank seq"
            + "\t" + "LTR2 in flank seq" + "\t" * 2
            + "seq_name_reset(map to the genome)" + "\t"
            + "LTR1 in genome" + "\t" + "LTR2 in genome" + "\n")

        for k, y in record_dict.items():

            qury_seq_name = ""

            if seqs == {}:
                LTRharvest_out_dir = os.path.dirname(self.gff3_path)

                LTRharvest_out_file = LTRharvest_out_dir + "/" + "LTRharvest_out.fas"

                with open(LTRharvest_out_file, "r", encoding='utf-8') as f:

                    for line in f:
                        if line.startswith(">"):
                            qury_seq_name = line.split("(dbseq")[0].strip().strip(">")

            else:
                qury_seq_name = seqs[k.split("|")[0]]  #  name of query sequence

            repeat_region_range = k.split("|")[1].strip()

            ltr1_range = y[1]  # range of ltr1
            ltr2_range = y[2]  # range of ltr2

            # print(y)

            """
            Recalculate the start and end positions of the endogenous 
            viral element on the original genome for renaming the sequence name 
            """

            # The start position of the flank sequence on the input genome
            start_nt = int(qury_seq_name.split("|")[1])

            # The end position of the flank sequence on the input genome
            end_nt = int(qury_seq_name.split("|")[2])
            seq_lenth = end_nt - start_nt + 1

            # start and end position in candidate pair-LTR ervs
            erv_start = int(repeat_region_range.split("..")[0])
            erv_end = int(repeat_region_range.split("..")[1])

            ltr1_start = int(repeat_region_range.split("..")[0])
            # print(ltr1_start)

            ltr1_end = int(ltr1_range.split("..")[1])

            ltr2_start = int(ltr2_range.split("..")[0])
            ltr2_end = int(repeat_region_range.split("..")[1])

            # whether it is on the complementary chain
            if qury_seq_name.split("|")[-1] == "complement":

                erv_start = seq_lenth - int(repeat_region_range.split("..")[1]) + 1
                erv_end = seq_lenth - int(repeat_region_range.split("..")[0]) + 1


                ltr1_start = seq_lenth - int(ltr1_range.split("..")[1]) + 1

                ltr1_end = seq_lenth - int(repeat_region_range.split("..")[0]) + 1

                ltr2_start = seq_lenth - int(repeat_region_range.split("..")[1]) + 1


                ltr2_end = seq_lenth - int(ltr2_range.split("..")[0]) + 1


            erv_start_site = start_nt + erv_start - 1
            erv_end_site = start_nt + erv_end - 1

            ltr1_start_site = start_nt + ltr1_start - 1
            ltr1_end_site = start_nt + ltr1_end - 1

            ltr2_start_site = start_nt + ltr2_start - 1
            ltr2_end_site = start_nt + ltr2_end - 1

            seq_name_reset = (
                        qury_seq_name.split("|")[0] + "|"
                        + str(erv_start_site) + "|"
                        + str(erv_end_site) + "|"
                        + qury_seq_name.split("|")[-1])


            out_simplify_file.write(
                qury_seq_name + "\t" + k.split("|")[1] + "\t" + y[0] + "\t"
                + str(int(repeat_region_range.split("..")[0])) + ".."
                + str(int(ltr1_range.split("..")[1])) + "\t"
                + str(int(ltr2_range.split("..")[0])) + ".."
                + str(int(repeat_region_range.split("..")[1]))

                + "\t" * 2 + seq_name_reset + "\t"
                + str(ltr1_start_site) + ".." + str(ltr1_end_site) + "\t"
                + str(ltr2_start_site) + ".." + str(ltr2_end_site) + "\n")


        out_simplify_file.close()

        out_duplicate_removal = open(self.duplicate_removal_path, "w",encoding="utf-8")

        simplify_gff3 = open(self.simplify_gff3_path, "r",encoding="utf-8")

        lines_seen = set()
        for line in simplify_gff3:
            if line.startswith("#"):
                out_duplicate_removal.write(line)
            else:
                seq_name = line.split("\t")[-3]

                if seq_name not in lines_seen:
                    out_duplicate_removal.write(line)

                    lines_seen.add(seq_name)

        simplify_gff3.close()
        out_duplicate_removal.close()



