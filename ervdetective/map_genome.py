# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/11 11:14

"""
import os

from public_function import domain_locate_genome

class PairLTRervMap(object):
    def __init__(self,hmmer_simplify_path,gff3_final_path, out_file):
        super(PairLTRervMap,self).__init__()
        
        self.hmmer_simplify_path = hmmer_simplify_path
        self.gff3_final_path = gff3_final_path
        self.domain_in_genome = out_file


    def run(self):

        HMMER_simplify_file = open(self.hmmer_simplify_path, "r",
                                   encoding="utf-8")
        domain_in_genome_file = open(self.domain_in_genome, "w",encoding="utf-8")

        domain_in_genome_file.write(
            "##Source of ERVs with paired_LTR "
            + "[host|start|end|chain]"
            + "\t" + "5'-LTR in genome" + "\t" + "3'-LTR in genome"
            + "\t" + "LTR_similarity" + "\t" + "GAG in genome" + "\t"
            + "DUT in genome" + "\t" + "AP in genome" + "\t"
            + "RT in genome" + "\t" + "RNaseH in genome" + "\t"
            + "INT in genome" + "\t" + "ENV in genome" + "\n")


        gff3_final_file = open(self.gff3_final_path, "r",encoding="utf-8")

        ltr_dic = {}

        for line in gff3_final_file:
            line = line.strip()
            if not line.startswith("#"):
                line_list = line.split("\t")
                pair_ltr_name = line_list[-3]
                ltrs = (line_list[-2] + "\t" + line_list[-1]
                        + "\t" + line_list[2])
                ltr_dic[pair_ltr_name] = ltrs

        for line in HMMER_simplify_file:
            line = line.strip()
            if not line.startswith("##"):
                line_list = line.split("\t")

                seq_name = line_list[0]

                ltrs = ltr_dic[seq_name]

                ltr1, ltr2, ltr_similarary = ltrs.split("\t")

                new_domain_list = []

                same_strand_count = 0
                complement_strand_count = 0

                for n in range(1, len(line_list)):
                    new_domain = "not found"
                    domain_record = line_list[n]
                    if domain_record != "not found":
                        (domain_name, c_Evalues, orfs,
                         domain_aa_site,
                         domain_nts_lenth) = domain_record.split(",", 4)

                        seq_coordinate = (str(seq_name.split("|")[1]) + ".."
                                          + str(seq_name.split("|")[2]))

                        seq_stand = seq_name.split("|")[-1]

                        re_strand, re_domain_location = domain_locate_genome(
                            seq_coordinate, seq_stand, domain_aa_site, orfs)

                        new_domain = (domain_name + "," + c_Evalues + ","
                                      + re_domain_location)

                        if re_strand == "same":
                            same_strand_count += 1
                        elif re_strand == "complement":
                            complement_strand_count += 1

                    new_domain_list.append(new_domain)

                if same_strand_count >= complement_strand_count:
                    finnaly_strand = "same"

                    ltr1_start = int(ltr1.split("..")[0])
                    ltr2_start = int(ltr2.split("..")[0])

                    if ltr1_start < ltr2_start:
                        five_ltr = ltr1  # Take the smallest value as 5’-LTR
                        three_ltr = ltr2

                    else:
                        five_ltr = ltr2  # Take the smallest value as 5’-LTR
                        three_ltr = ltr1


                else:
                    finnaly_strand = "complement"

                    ltr1_start = int(ltr1.split("..")[0])
                    ltr2_start = int(ltr2.split("..")[0])

                    if ltr1_start < ltr2_start:
                        # The complementary chain takes the larger value as 5’-LTR
                        five_ltr = "complement(" + ltr2 + ")"
                        three_ltr = "complement(" + ltr1 + ")"

                    else:
                        five_ltr = "complement(" + ltr1 + ")"
                        three_ltr = "complement(" + ltr2 + ")"


                new_line = (seq_name.replace(seq_name.split("|")[-1],
                                             "").strip("|")
                            + "|" + finnaly_strand + "\t" + five_ltr
                            + "\t" + three_ltr
                            + "\t" + ltr_similarary + "\t"
                            + "\t".join(new_domain_list) + "\n")

                domain_in_genome_file.write(new_line)


        gff3_final_file.close()
        domain_in_genome_file.close()




class NoPairLTRervMap(object):
    def __init__(self, hmmer_simplify_path, out_map,
                 out_map_duplicate_removal):
        super(NoPairLTRervMap,self).__init__()
        self.hmmer_simplify_path = hmmer_simplify_path
        self.domain_in_genome = out_map
        self.out_map_duplicate_removal = out_map_duplicate_removal


    def run(self):

        HMMER_simplify_file = open(self.hmmer_simplify_path, "r",encoding="utf-8")

        domain_in_genome_file = open(self.domain_in_genome, "w",encoding="utf-8")
        domain_in_genome_file.write(
            "##flank sequence used for HMMER"
            + "\t"
            + "potential ERVs without paired-LTRs "
            + "[host|start|end|chain]"
            + "\t" + "GAG in genome" + "\t"
            + "DUT in genome" + "\t" + "AP in genome" + "\t"
            + "RT in genome" + "\t" + "RNaseH in genome" + "\t"
            + "INT in genome" + "\t" + "ENV in genome" + "\n")


        domain_duplicate_removal = open(self.out_map_duplicate_removal, "w",
                                        encoding="utf-8")
        domain_duplicate_removal.write(
            "##potential ERVs without paired-LTRs "
            + "[host|start|end|chain]"
            + "\t" + "GAG in genome" + "\t"
            + "DUT in genome" + "\t" + "AP in genome" + "\t"
            + "RT in genome" + "\t" + "RNaseH in genome" + "\t"
            + "INT in genome" + "\t" + "ENV in genome" + "\n")

        for line in HMMER_simplify_file:
            line = line.strip()
            if not line.startswith("##"):
                line_list = line.split("\t")

                seq_name = line_list[0]

                new_domain_list = []

                same_strand_count = 0
                complement_strand_count = 0

                max_site_list = []
                min_site_list = []

                for n in range(1, len(line_list)):
                    new_domain = "not found"
                    domain_record = line_list[n]
                    if domain_record != "not found":
                        (domain_name, c_Evalues, orfs,
                         domain_aa_site,
                         domain_nts_lenth) = domain_record.split(",", 4)

                        seq_coordinate = (str(seq_name.split("|")[1]) + ".."
                                          + str(seq_name.split("|")[2]))

                        seq_stand = seq_name.split("|")[-1]

                        re_strand, re_domain_location = domain_locate_genome(
                            seq_coordinate, seq_stand, domain_aa_site, orfs)

                        new_domain = (domain_name + "," + c_Evalues + ","
                                      + re_domain_location)

                        domain_range = \
                        (re_domain_location.split("(")[-1]).split(")")[0]
                        domain_start = int(domain_range.split("..")[0])
                        domain_end = int(domain_range.split("..")[1])

                        max_site_list.append(domain_end)
                        min_site_list.append(domain_start)

                        if re_strand == "same":
                            same_strand_count += 1
                        elif re_strand == "complement":
                            complement_strand_count += 1

                    new_domain_list.append(new_domain)

                if same_strand_count >= complement_strand_count:
                    finnaly_strand = "same"

                else:
                    finnaly_strand = "complement"

                max_site_list.sort()
                min_site_list.sort()

                max_site = max_site_list[-1]
                min_site = min_site_list[0]

                new_line = (seq_name + "\t"
                            + seq_name.split("|")[0] + "|"
                            + str(min_site) + "|" + str(max_site)
                            + "|" + finnaly_strand + "\t"
                            + "\t".join(new_domain_list) + "\n")

                domain_in_genome_file.write(new_line)

        domain_in_genome_file.close()


        """
        Remove potential regions of erv that overlap with 
        the scope after relocation 
        """
        domain_in_genome_file = open(self.domain_in_genome, "r",
                                     encoding="utf-8")

        old_seq_name = []

        name_dic = {}
        use_list = []

        for line in domain_in_genome_file:
            if not line.startswith("##"):
                seq_name = line.split("\t")[1]
                old_seq_name.append(seq_name.strip())

                name_dic[seq_name] = (line.split("\t")[1] + "\t" + 
                                      line.split("\t")[2] + "\t" +
                                      line.split("\t")[3] + "\t" +
                                      line.split("\t")[4] + "\t" +
                                      line.split("\t")[5] + "\t" +
                                      line.split("\t")[6] + "\t" +
                                      line.split("\t")[7] + "\t" +
                                      line.split("\t")[8]).strip()

                                                 

        for line1 in old_seq_name:

            qury_name, qury_start, qury_end, qury_direction = line1.split("|",
                                                                          3)

            for line2 in old_seq_name:
                if line1 != line2:

                    (refer_name, refer_start,
                     refer_end, refer_direction) = line2.split("|", 3)
                    if (qury_name == refer_name
                        and qury_direction) == refer_direction:
                        if qury_start >= refer_start and qury_end <= refer_end:
                            break

            else:
                use_list.append(line1)

        """
        Annotated results in the same area but on a different chain, 
        one is kept at random 
        """
        ganjing = set()

        finnaly_use = []

        for each_seqname in use_list:
            labels = each_seqname.replace(each_seqname.split("|")[-1], "")
            if not labels in ganjing:
                finnaly_use.append(each_seqname)
                ganjing.add(labels)

        for line in finnaly_use:
            domain_duplicate_removal.write(name_dic[line] + "\n")

        domain_duplicate_removal.close()

        domain_in_genome_file.close()

        os.remove(self.domain_in_genome)

