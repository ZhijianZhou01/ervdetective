# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/11 11:14

"""
import os

from public_function import domain_locate_genome

def domain_start_end(domian):


    if domian.find("complement") == -1:  # 正链上

        domian_start = int((domian.split("(")[-1]).split("..")[0])

        domian_end = int((domian.split("(")[-1]).split("..")[-1].replace(")", ""))

    else: # 负链上
        domian_start = int(
            (domian.split("(")[-1]).split("..")[-1].replace(")", ""))

        domian_end = int((domian.split("(")[-1]).split("..")[0])


    return [domian_start,domian_end]




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
            "##Source of ERVs with paired_LTRs "
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

        gff3_final_file.close()

        domain_in_genome_dic = {}

        name_list = []

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
                    domain_record = line_list[n]
                    if domain_record != "not found" and domain_record.find("relatively short") == -1:
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

                    else:
                        new_domain = domain_record

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


                name_ = seq_name.replace(seq_name.split("|")[-1],"").strip("|") + "|" + finnaly_strand

                new_line = ( name_ + "\t" + five_ltr
                            + "\t" + three_ltr
                            + "\t" + ltr_similarary + "\t"
                            + "\t".join(new_domain_list))

                name_list.append(name_)
                domain_in_genome_dic[name_] = new_line




        use_list = []
        for line1 in name_list:

            qury_name, qury_start, qury_end, qury_direction = line1.split("|",
                                                                          3)

            for line2 in name_list:
                if line1 != line2:

                    (refer_name, refer_start,
                     refer_end, refer_direction) = line2.split("|", 3)
                    if (qury_name == refer_name
                        and qury_direction) == refer_direction:
                        if qury_start >= refer_start and qury_end <= refer_end:
                            break

            else:
                use_list.append(line1)




        woqu_ = {}

        finnaly_use = []

        for each_name in use_list:
            erv_contain_list = domain_in_genome_dic[each_name].split("\t")  

            xxx_list = []

            bat_seq_index = erv_contain_list[0].split("|")[0]  

            ltrs_similary = float(erv_contain_list[3]) 

            xxx_list.append(bat_seq_index)

            for n in range(4,len(erv_contain_list)):

                domain_ = erv_contain_list[n]

                if domain_ != "not found" and domain_.find("relatively short") == -1:
                    xxx_list.append(domain_.split(",")[-1]) 

                else:

                    xxx_list.append(domain_)

            label_ = "\t".join(xxx_list) 

            if not label_ in woqu_: 
                finnaly_use.append(each_name) 
                woqu_[label_] = [each_name,ltrs_similary]

            else:
                if ltrs_similary > woqu_[label_][1]: 
                    finnaly_use.remove(woqu_[label_][0])
                    finnaly_use.append(each_name) 
                    woqu_[label_] = [each_name, ltrs_similary] 


        for erv_name in finnaly_use:

            flage = 0

            poss_erv = domain_in_genome_dic[erv_name]

            poss_erv_contain_list = poss_erv.split("\t")

            erv_name = poss_erv_contain_list[0]

            chain_direction = erv_name.split("|")[-1]

            GAG_ = poss_erv_contain_list[4].strip()

            PRO_ = poss_erv_contain_list[6].strip()

            RT_ = poss_erv_contain_list[7].strip()

            RH_ = poss_erv_contain_list[8].strip()

            INT_ = poss_erv_contain_list[9].strip()

            ENV_ = poss_erv_contain_list[-1].strip()

            if GAG_.find(",") != -1 and RT_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)

                RT_site = domain_start_end(RT_)

                if chain_direction == "same" and RT_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and RT_site[0] > Gag_site[
                    0]:

                    flage = 1

            if GAG_.find(",") != -1 and RH_.find(",") != -1:
                Gag_site = domain_start_end(GAG_)
                RH_site = domain_start_end(RH_)

                if chain_direction == "same" and RH_site[0] < Gag_site[0]:
                    flage = 1

                elif chain_direction == "complement" and RH_site[0] > Gag_site[
                    0]:
                    flage = 1



            if GAG_.find(",") != -1 and INT_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)

                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and INT_site[0] > Gag_site[0]:

                    flage = 1

            if GAG_.find(",") != -1 and ENV_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)
                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and ENV_site[0] > Gag_site[
                    0]:

                    flage = 1

            if RT_.find(",") != -1 and ENV_.find(",") != -1:
                RT_site = domain_start_end(RT_)
                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < RT_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > RT_site[
                    0]:

                    flage = 1


            if RT_.find(",") != -1 and RH_.find(",") != -1:
                RT_site = domain_start_end(RT_)
                RH_site = domain_start_end(RH_)

                if chain_direction == "same" and RT_site[0] > RH_site[0]: 
                    flage = 1

                elif chain_direction == "complement" and RT_site[0] < RH_site[0]:
                    flage = 1

                if abs(RT_site[0] - RH_site[0]) > 5000: 
                    flage = 1


            if RT_.find(",") != -1 and INT_.find(",") != -1:


                RT_site = domain_start_end(RT_)
                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < RT_site[0]:  
                    flage = 1

                elif chain_direction == "complement" and INT_site[0] > RT_site[0]:
                    flage = 1


            if RH_.find(",") != -1 and INT_.find(",") != -1:
                RH_site =  domain_start_end(RH_)
                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < RH_site[0]:  
                    flage = 1

                elif chain_direction == "complement" and INT_site[0] > RH_site[0]:
                    flage = 1



            if PRO_.find(",") != -1 and RT_.find(",") != -1:
                PRO_site = domain_start_end(PRO_)
                RT_site = domain_start_end(RT_)

                if chain_direction == "same" and RT_site[0] < PRO_site[0]:
                    flage = 1

                elif chain_direction == "complement" and RT_site[0] > PRO_site[0]:
                    flage = 1


            if PRO_.find(",") != -1 and ENV_.find(",") != -1:
                PRO_site = domain_start_end(PRO_)

                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < PRO_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > PRO_site[
                    0]:

                    flage = 1

            if INT_.find(",") != -1 and ENV_.find(",") != -1:
                INT_site = domain_start_end(INT_)

                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < INT_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > INT_site[
                    0]:

                    flage = 1

            if flage == 0:


                domain_in_genome_file.write(domain_in_genome_dic[erv_name] + "\n")


        domain_in_genome_file.close()




class NoPairLTReveMap(object):
    def __init__(self, pairltr_ervs_map_genome,hmmer_simplify_path, out_map,
                 out_map_duplicate_removal):
        super(NoPairLTReveMap,self).__init__()

        self.pairltr_ervs_map_genome = pairltr_ervs_map_genome
        self.hmmer_simplify_path = hmmer_simplify_path
        self.domain_in_genome = out_map
        self.out_map_duplicate_removal = out_map_duplicate_removal


    def run(self):

        HMMER_simplify_file = open(self.hmmer_simplify_path, "r",encoding="utf-8")

        domain_in_genome_file = open(self.domain_in_genome, "w",encoding="utf-8")
        domain_in_genome_file.write(
            "##flank sequence used for HMMER"
            + "\t"
            + "potential ERVLEs without paired-LTRs "
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

                    domain_record = line_list[n]

                    if domain_record != "not found" and domain_record.find("relatively short") == -1:
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

                    else:
                        new_domain = domain_record

                    new_domain_list.append(new_domain)

                if same_strand_count >= complement_strand_count:
                    finnaly_strand = "same"

                else:
                    finnaly_strand = "complement"

                max_site_list.sort()
                min_site_list.sort()

                max_site = max_site_list[-1]
                min_site = min_site_list[0]

                if abs(max_site - min_site) <= 12000: 

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

        ervles_dic = {}
        use_list = []

        for line in domain_in_genome_file:
            if not line.startswith("##"):
                seq_name = line.split("\t")[1]
                old_seq_name.append(seq_name.strip())

                ervles_dic[seq_name] = (line.split("\t")[1] + "\t" + 
                                      line.split("\t")[2] + "\t" +
                                      line.split("\t")[3] + "\t" +
                                      line.split("\t")[4] + "\t" +
                                      line.split("\t")[5] + "\t" +
                                      line.split("\t")[6] + "\t" +
                                      line.split("\t")[7] + "\t" +
                                      line.split("\t")[8]).strip()


        domain_in_genome_file.close()

        os.remove(self.domain_in_genome)



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


        ganjing = set()

        finnaly_use1 = []

        for each_seqname in use_list:
            if not each_seqname in ganjing:
                finnaly_use1.append(each_seqname)
                ganjing.add(each_seqname)


  
        woqu_ = set()

        finnaly_use2 = []

        for each_seqname in finnaly_use1:
            erv_contain_list = ervles_dic[each_seqname].split("\t") 

            xxx_list = []

            bat_seq_index = erv_contain_list[0].split("|")[0] 

            xxx_list.append(bat_seq_index)

            for n in range(1,len(erv_contain_list)):

                domain_ = erv_contain_list[n]

                if domain_ != "not found" and domain_.find("relatively short") == -1:
                    xxx_list.append(domain_.split(",")[-1]) 

                else:

                    xxx_list.append(domain_)

            label_ = "\t".join(xxx_list)

            if not label_ in woqu_:
                finnaly_use2.append(each_seqname)
                woqu_.add(label_)



        ERV_anno_domian_dic = {}

        with open(self.pairltr_ervs_map_genome,"r",encoding="utf-8") as pairltr_ervs_anno:

            for line in pairltr_ervs_anno:

                if not line.startswith("##"):

                    line = line.strip()

                    line_list = line.split("\t")

                    key = line_list[0].split("|")[0] + "_" + \
                          line_list[0].split("|")[-1]

                    if key not in ERV_anno_domian_dic:
                        ERV_anno_domian_dic[key] = []

                    for n in range(4, len(line_list)):
                        domain = line_list[n]
                        if domain.find("relatively short") == -1 and domain != "not found":
                            ERV_anno_domian_dic[key].append(domain.split(",")[-1])



        finnaly_use3 = []

        for ervle_name in finnaly_use2:

            flage = 0

            line = ervles_dic[ervle_name]

            line = line.strip()

            line_list = line.split("\t")

            key = line_list[0].split("|")[0] + "_" + line_list[0].split("|")[-1]


            for n in range(1, len(line_list)):
                domain = line_list[n]

                if domain.find("relatively short") == -1 and domain != "not found":

                    local_ = domain.split(",")[-1]

                    if key in ERV_anno_domian_dic:

                        ERV_anno_domian_list = ERV_anno_domian_dic[key]

                        if local_ in ERV_anno_domian_list:

                            flage = 1

                            break

            if flage == 0:
                finnaly_use3.append(ervle_name)


        domain_duplicate_removal = open(self.out_map_duplicate_removal, "w",
                                        encoding="utf-8")

        domain_duplicate_removal.write(
            "##potential ERVLEs without paired-LTRs "
            + "[host|start|end|chain]"
            + "\t" + "GAG in genome" + "\t"
            + "DUT in genome" + "\t" + "AP in genome" + "\t"
            + "RT in genome" + "\t" + "RNaseH in genome" + "\t"
            + "INT in genome" + "\t" + "ENV in genome" + "\n")



        for line in finnaly_use3:

            flage = 0

            poss_erv = ervles_dic[line]

            poss_erv_contain_list = poss_erv.split("\t")

            erv_name = poss_erv_contain_list[0]

            chain_direction = erv_name.split("|")[-1]

            GAG_ = poss_erv_contain_list[1].strip()

            PRO_ = poss_erv_contain_list[3].strip()

            RT_ = poss_erv_contain_list[4].strip()

            RH_ = poss_erv_contain_list[5].strip()

            INT_ = poss_erv_contain_list[6].strip()

            ENV_ = poss_erv_contain_list[-1].strip()


            if GAG_.find(",") != -1 and RT_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)

                RT_site = domain_start_end(RT_)

                if chain_direction == "same" and RT_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and RT_site[0] > Gag_site[0]:

                    flage = 1

            if GAG_.find(",") != -1 and RH_.find(",") != -1:
                Gag_site = domain_start_end(GAG_)
                RH_site = domain_start_end(RH_)

                if chain_direction == "same" and RH_site[0] < Gag_site[0]:
                    flage = 1

                elif chain_direction == "complement" and RH_site[0] > Gag_site[0]:
                    flage = 1


            if GAG_.find(",") != -1 and INT_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)

                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and INT_site[0] > Gag_site[0]:

                    flage = 1



            if GAG_.find(",") != -1 and ENV_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)
                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and ENV_site[0] > Gag_site[0]:

                    flage = 1


                if abs(Gag_site[0] - ENV_site[0]) > 12000:  
                    flage = 1




            if RT_.find(",") != -1 and ENV_.find(",") != -1:
                RT_site = domain_start_end(RT_)
                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < RT_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > RT_site[0]:

                    flage = 1

                if abs(RT_site[0] - ENV_site[0]) > 6000:  
                    flage = 1



            if PRO_.find(",") != -1 and RT_.find(",") != -1:
                PRO_site = domain_start_end(PRO_)
                RT_site = domain_start_end(RT_)

                if chain_direction == "same" and RT_site[0] < PRO_site[0]:
                    flage = 1

                elif chain_direction == "complement" and RT_site[0] > PRO_site[0]:
                    flage = 1



            if RT_.find(",") != -1 and INT_.find(",") != -1:


                RT_site = domain_start_end(RT_)
                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < RT_site[
                    0]:  
                    flage = 1

                elif chain_direction == "complement" and INT_site[0] > RT_site[0]:
                    flage = 1


            if RH_.find(",") != -1 and INT_.find(",") != -1:
                RH_site =  domain_start_end(RH_)
                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < RH_site[0]: 
                    flage = 1

                elif chain_direction == "complement" and INT_site[0] > RH_site[0]:
                    flage = 1



            if RT_.find(",") != -1 and RH_.find(",") != -1:
                RT_site = domain_start_end(RT_)
                RH_site = domain_start_end(RH_)

                if chain_direction == "same" and RT_site[0] > RH_site[0]: 
                    flage = 1

                elif chain_direction == "complement" and RT_site[0] < RH_site[0]:
                    flage = 1


                if abs(RT_site[0] - RH_site[0]) > 5000: 
                    flage = 1


            if PRO_.find(",") != -1 and ENV_.find(",") != -1:
                PRO_site = domain_start_end(PRO_)

                ENV_site = domain_start_end(ENV_)


                if chain_direction == "same" and ENV_site[0] < PRO_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > PRO_site[0]:

                    flage = 1


            if INT_.find(",") != -1 and ENV_.find(",") != -1:
                INT_site = domain_start_end(INT_)

                ENV_site = domain_start_end(ENV_)


                if chain_direction == "same" and ENV_site[0] < INT_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > INT_site[0]:

                    flage = 1


            if flage == 0:

                domain_duplicate_removal.write(ervles_dic[erv_name] + "\n")

        domain_duplicate_removal.close()
        





class NoPairLTReveMap2(object):

    """
    if not ERV with paired LTRs in genome

    """

    def __init__(self, hmmer_simplify_path, out_map,out_map_duplicate_removal):
        super(NoPairLTReveMap2,self).__init__()

        self.hmmer_simplify_path = hmmer_simplify_path
        self.domain_in_genome = out_map
        self.out_map_duplicate_removal = out_map_duplicate_removal


    def run(self):

        HMMER_simplify_file = open(self.hmmer_simplify_path, "r",encoding="utf-8")

        domain_in_genome_file = open(self.domain_in_genome, "w",encoding="utf-8")
        domain_in_genome_file.write(
            "##flank sequence used for HMMER"
            + "\t"
            + "potential ERVLEs without paired-LTRs "
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

                    domain_record = line_list[n]

                    if domain_record != "not found" and domain_record.find("relatively short") == -1:
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

                    else:
                        new_domain = domain_record

                    new_domain_list.append(new_domain)

                if same_strand_count >= complement_strand_count:
                    finnaly_strand = "same"

                else:
                    finnaly_strand = "complement"

                max_site_list.sort()
                min_site_list.sort()

                max_site = max_site_list[-1]
                min_site = min_site_list[0]

                if abs(max_site - min_site) <= 12000:

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

        ervles_dic = {}
        use_list = []

        for line in domain_in_genome_file:
            if not line.startswith("##"):
                seq_name = line.split("\t")[1]
                old_seq_name.append(seq_name.strip())

                ervles_dic[seq_name] = (line.split("\t")[1] + "\t" +
                                      line.split("\t")[2] + "\t" +
                                      line.split("\t")[3] + "\t" +
                                      line.split("\t")[4] + "\t" +
                                      line.split("\t")[5] + "\t" +
                                      line.split("\t")[6] + "\t" +
                                      line.split("\t")[7] + "\t" +
                                      line.split("\t")[8]).strip()


        domain_in_genome_file.close()

        os.remove(self.domain_in_genome)


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


        ganjing = set()

        finnaly_use1 = []

        for each_seqname in use_list:
            if not each_seqname in ganjing:
                finnaly_use1.append(each_seqname)
                ganjing.add(each_seqname)


        woqu_ = set()

        finnaly_use2 = []

        for each_seqname in finnaly_use1:
            erv_contain_list = ervles_dic[each_seqname].split("\t") 

            xxx_list = []

            bat_seq_index = erv_contain_list[0].split("|")[0]

            xxx_list.append(bat_seq_index)

            for n in range(1,len(erv_contain_list)):

                domain_ = erv_contain_list[n]

                if domain_ != "not found" and domain_.find("relatively short") == -1:
                    xxx_list.append(domain_.split(",")[-1]) 

                else:

                    xxx_list.append(domain_)

            label_ = "\t".join(xxx_list)

            if not label_ in woqu_:
                finnaly_use2.append(each_seqname)
                woqu_.add(label_)




        domain_duplicate_removal = open(self.out_map_duplicate_removal, "w",
                                        encoding="utf-8")

        domain_duplicate_removal.write(
            "##potential ERVLEs without paired-LTRs "
            + "[host|start|end|chain]"
            + "\t" + "GAG in genome" + "\t"
            + "DUT in genome" + "\t" + "AP in genome" + "\t"
            + "RT in genome" + "\t" + "RNaseH in genome" + "\t"
            + "INT in genome" + "\t" + "ENV in genome" + "\n")



        for line in finnaly_use2:

            flage = 0

            poss_erv = ervles_dic[line]

            poss_erv_contain_list = poss_erv.split("\t")

            erv_name = poss_erv_contain_list[0]

            chain_direction = erv_name.split("|")[-1]

            GAG_ = poss_erv_contain_list[1].strip()

            PRO_ = poss_erv_contain_list[3].strip()

            RT_ = poss_erv_contain_list[4].strip()

            RH_ = poss_erv_contain_list[5].strip()

            INT_ = poss_erv_contain_list[6].strip()

            ENV_ = poss_erv_contain_list[-1].strip()


            if GAG_.find(",") != -1 and RT_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)

                RT_site = domain_start_end(RT_)

                if chain_direction == "same" and RT_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and RT_site[0] > Gag_site[0]:

                    flage = 1

            if GAG_.find(",") != -1 and RH_.find(",") != -1:
                Gag_site = domain_start_end(GAG_)
                RH_site = domain_start_end(RH_)

                if chain_direction == "same" and RH_site[0] < Gag_site[0]:
                    flage = 1

                elif chain_direction == "complement" and RH_site[0] > Gag_site[0]:
                    flage = 1


            if GAG_.find(",") != -1 and INT_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)

                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and INT_site[0] > Gag_site[0]:

                    flage = 1



            if GAG_.find(",") != -1 and ENV_.find(",") != -1:

                Gag_site = domain_start_end(GAG_)
                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < Gag_site[0]:
                    flage = 1


                elif chain_direction == "complement" and ENV_site[0] > Gag_site[0]:

                    flage = 1


                if abs(Gag_site[0] - ENV_site[0]) > 12000:
                    flage = 1




            if RT_.find(",") != -1 and ENV_.find(",") != -1:
                RT_site = domain_start_end(RT_)
                ENV_site = domain_start_end(ENV_)

                if chain_direction == "same" and ENV_site[0] < RT_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > RT_site[0]:

                    flage = 1

                if abs(RT_site[0] - ENV_site[0]) > 6000:
                    flage = 1



            if PRO_.find(",") != -1 and RT_.find(",") != -1:
                PRO_site = domain_start_end(PRO_)
                RT_site = domain_start_end(RT_)

                if chain_direction == "same" and RT_site[0] < PRO_site[0]:
                    flage = 1

                elif chain_direction == "complement" and RT_site[0] > PRO_site[0]:
                    flage = 1



            if RT_.find(",") != -1 and INT_.find(",") != -1:

                RT_site = domain_start_end(RT_)
                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < RT_site[
                    0]:
                    flage = 1

                elif chain_direction == "complement" and INT_site[0] > RT_site[0]:
                    flage = 1


            if RH_.find(",") != -1 and INT_.find(",") != -1:
                RH_site =  domain_start_end(RH_)
                INT_site = domain_start_end(INT_)

                if chain_direction == "same" and INT_site[0] < RH_site[0]:
                    flage = 1

                elif chain_direction == "complement" and INT_site[0] > RH_site[0]:
                    flage = 1



            if RT_.find(",") != -1 and RH_.find(",") != -1:
                RT_site = domain_start_end(RT_)
                RH_site = domain_start_end(RH_)

                if chain_direction == "same" and RT_site[0] > RH_site[0]: 
                    flage = 1

                elif chain_direction == "complement" and RT_site[0] < RH_site[0]:
                    flage = 1


                if abs(RT_site[0] - RH_site[0]) > 5000:
                    flage = 1


            if PRO_.find(",") != -1 and ENV_.find(",") != -1:
                PRO_site = domain_start_end(PRO_)

                ENV_site = domain_start_end(ENV_)


                if chain_direction == "same" and ENV_site[0] < PRO_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > PRO_site[0]:

                    flage = 1


            if INT_.find(",") != -1 and ENV_.find(",") != -1:
                INT_site = domain_start_end(INT_)

                ENV_site = domain_start_end(ENV_)


                if chain_direction == "same" and ENV_site[0] < INT_site[0]:

                    flage = 1

                elif chain_direction == "complement" and ENV_site[0] > INT_site[0]:

                    flage = 1


            if flage == 0:

                domain_duplicate_removal.write(ervles_dic[erv_name] + "\n")

        domain_duplicate_removal.close()

