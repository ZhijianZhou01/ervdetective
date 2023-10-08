# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/10 23:03

"""

from public_function import best_domain

class SimplyHMMER(object):
    def __init__(self, hmmer_file, out_simplify, header_str):
        super(SimplyHMMER,self).__init__()

        self.hmmer_result_path = hmmer_file
        self.hmmer_simplify_path = out_simplify
        self.header_str = header_str

    def run(self):

        hmmer_result_file = open(self.hmmer_result_path, "r",encoding="utf-8")

        HMMER_simplify_file = open(self.hmmer_simplify_path, "w",
                                   encoding="utf-8")

        HMMER_simplify_file.write(self.header_str)

        sum_dict = {}

        for line in hmmer_result_file:

            if not line.startswith("#"):
                line_list = [i for i in line.strip().split(' ') if i != ""]

                # print(line_list)

                new_line = []

                domain = line_list[0]  # record of predicted domain

                ERV_candidate_seqname = (line_list[3].split("#")[0])

                orf = line_list[3].split("#")[-1]

                c_Evalue = float(line_list[11])

                # print(c_Evalue)

                start_aa = (int(line_list[17]))

                end_aa = (int(line_list[18]))

                domain_aa_lenth = end_aa - start_aa + 1

                new_line.append(domain)
                new_line.append(c_Evalue)
                new_line.append(orf)
                new_line.append(
                    "AA(" + str(start_aa) + ".." + str(end_aa) + ")")
                new_line.append(domain_aa_lenth)

                if not sum_dict.__contains__(ERV_candidate_seqname):
                    sum_dict[ERV_candidate_seqname] = []
                sum_dict[ERV_candidate_seqname].append(new_line)

        for key in sum_dict:

            line_flage = "False"

            # print(key, sum_dict[key])
            hit_list = sum_dict[key]
            GAG_list = []
            DUT_list = []
            AP_list = []
            RT_lsit = []
            RNaseH_list = []
            INT_list = []
            ENV_list = []

            for each_hit in hit_list:

                # print(each_hit)
                match_domain = each_hit[0]

                if match_domain.split("_")[0] == "GAG":
                    GAG_list.append(each_hit)

                elif match_domain.split("_")[0] == "DUT":
                    DUT_list.append(each_hit)

                elif match_domain.split("_")[0] == "AP":
                    AP_list.append(each_hit)

                elif match_domain.split("_")[0] == "RT":
                    RT_lsit.append(each_hit)

                elif match_domain.split("_")[0] == "RNaseH":
                    RNaseH_list.append(each_hit)

                elif match_domain.split("_")[0] == "INT":
                    INT_list.append(each_hit)

                elif match_domain.split("_")[0] == "ENV":
                    ENV_list.append(each_hit)

            hit_GAG = best_domain(GAG_list, 1)
            if hit_GAG != "not found":
                line_flage = "True"

            hit_DUT = best_domain(DUT_list, 1)
            if hit_DUT != "not found":
                line_flage = "True"

            hit_AP = best_domain(AP_list, 1)
            if hit_AP != "not found":
                line_flage = "True"

            hit_RT = best_domain(RT_lsit, 1)
            if hit_RT != "not found":
                line_flage = "True"

            hit_RNaseH = best_domain(RNaseH_list, 1)
            if hit_RNaseH != "not found":
                line_flage = "True"

            hit_INT = best_domain(INT_list, 1)
            if hit_INT != "not found":
                line_flage = "True"

            hit_ENV = best_domain(ENV_list, 1)
            if hit_ENV != "not found":
                line_flage = "True"

            if line_flage == "True":
                HMMER_simplify_file.write(
                    key + "\t" + hit_GAG + "\t" + hit_DUT + "\t" + hit_AP
                    + "\t" + hit_RT + "\t" + hit_RNaseH + "\t" + hit_INT
                    + "\t" + hit_ENV + "\n")

        hmmer_result_file.close()
        HMMER_simplify_file.close()
