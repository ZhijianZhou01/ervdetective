# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/11 22:22

"""

class PairLTRErvAnno(object):
    def __init__(self, pairltr_ervs_map_genome, pairltr_ervs_anno, out_refix):

        """
        annotate domain  on the extracted pair of LTR_ervs
        :param pairltr_ervs_map_genome:
        :param pairltr_ervs_anno:
        :param out_refix:
        """
        super(PairLTRErvAnno,self).__init__()

        self.pairltr_ervs_map_genome = pairltr_ervs_map_genome
        self.pairltr_ervs_anno = pairltr_ervs_anno
        self.out_refix = out_refix
    
    def run(self):            

        domain_in_genome_file = open(self.pairltr_ervs_map_genome, "r",
                                     encoding="utf-8")
       
        out_domain_annos = open(self.pairltr_ervs_anno,"w",encoding="utf-8")

        out_domain_annos.write("##Name of ERVs with paired-LTRs" + "\t"
                               + "5'-LTR region" + "\t"
                               + "3'-LTR region" + "\t"
                               + "GAG" + "\t"
                               + "DUT" + "\t"
                               + "AP" + "\t"
                               + "RT" + "\t"
                               + "RNaseH" + "\t"
                               + "INT"+ "\t"
                               + "ENV" + "\n")


        for line in domain_in_genome_file:

            if not line.startswith("##"):

                line = line.strip()
                line_list = line.split("\t")

                seq_name = line_list[0]
                strand = seq_name.split("|")[-1]

                flage = 0

                line_record_list = []

                if strand == "same":
                    start_site = int(seq_name.split("|")[1])
                    # end_site =int(seq_name.split("|")[2])
                    for n in range(1, len(line_list)):
                        if n == 1 or n == 2:
                            sites = line_list[n]

                            LTR_start = int(
                                sites.split("..")[0]) - start_site + 1

                            LTR_end = int(sites.split("..")[1]) - start_site + 1

                            line_record_list.append(str(LTR_start) + ".."
                                                    + str(LTR_end))

                        elif n == 3:
                            pass

                        else:
                            if line_list[n] != "not found" and line_list[n].find("relatively short") == -1:
                                xx = (line_list[n].split(",")[-1]).split("(")[0]
                                if xx == "":
                                    sites = \
                                    (line_list[n].split("(")[1]).split(")")[0]
                                    domain_start = int(
                                        sites.split("..")[0]) - start_site + 1
                                    domain_end = int(
                                        sites.split("..")[-1]) - start_site + 1
                                    line_record_list.append(
                                        str(domain_start) + ".."
                                        + str(domain_end))
                                elif xx == "complement":
                                    sites = \
                                    (line_list[n].split("(")[1]).split(")")[0]
                                    domain_start = int(
                                        sites.split("..")[0]) - start_site + 1
                                    domain_end = int(
                                        sites.split("..")[-1]) - start_site + 1
                                    line_record_list.append("complement" + "("
                                                            + str(
                                        domain_start) + ".."
                                                            + str(
                                        domain_end) + ")")
                                    flage = 1

                            else:
                                line_record_list.append(line_list[n])


                elif strand == "complement":
                    start_site = int(seq_name.split("|")[2])
                    for n in range(1, len(line_list)):

                        if n == 1 or n == 2:

                            ltr_sites = (line_list[n].split("(")[1]).split(")")[
                                0]

                            LTR_start = start_site - int(
                                ltr_sites.split("..")[1]) + 1

                            LTR_end = start_site - int(
                                ltr_sites.split("..")[0]) + 1

                            line_record_list.append(str(LTR_start) + ".."
                                                    + str(LTR_end))
                        elif n == 3:
                            pass

                        else:
                            if line_list[n] != "not found" and line_list[n].find("relatively short") == -1:
                                xx = (line_list[n].split(",")[-1]).split("(")[0]
                                if xx == "complement":
                                    sites = \
                                    (line_list[n].split("(")[1]).split(")")[0]
                                    domain_start = start_site - int(
                                        sites.split("..")[1]) + 1
                                    domain_end = start_site - int(
                                        sites.split("..")[0]) + 1

                                    line_record_list.append(
                                        str(domain_start) + ".."
                                        + str(domain_end))
                                elif xx == "":
                                    sites = \
                                    (line_list[n].split("(")[1]).split(")")[0]
                                    domain_start = start_site - int(
                                        sites.split("..")[1]) + 1
                                    domain_end = start_site - int(
                                        sites.split("..")[0]) + 1

                                    line_record_list.append("complement" + "("
                                                            + str(
                                        domain_start) + ".."
                                                            + str(
                                        domain_end) + ")")
                                    flage = 1

                            else:
                                line_record_list.append(line_list[n])

                new_line = "\t".join(line_record_list)

                out_domain_annos.write(seq_name + "\t"
                    + new_line + "\n")

                if flage == 1:
                    out_domain_annos.write(
                        "##Warn: the annotation and chain direction of '" +
                        seq_name
                        + "' needs to be checked in file of 'domain_annotation_in_genome_final.txt'!" + "\n")


        domain_in_genome_file.close()
        out_domain_annos.close()



class PotentialEveAnno(object):
    def __init__(self, map_duplicate_removal,
                 domain_in_potential_erv,
                 out_refix):
        self.out_map_duplicate_removal = map_duplicate_removal
        self.domain_in_potential_erv = domain_in_potential_erv
        self.out_refix = out_refix

        super(PotentialEveAnno, self).__init__()

    def run(self):

        domain_duplicate_removal = open(self.out_map_duplicate_removal, "r",
                                        encoding="utf-8")
        out_domain_annos = open(self.domain_in_potential_erv, "w",
                                encoding="utf-8")

        out_domain_annos.write("##Name of potential ERVLEs without paired-LTRs" + "\t"
                               + "GAG" + "\t"
                               + "DUT" + "\t"
                               + "AP" + "\t"
                               + "RT" + "\t"
                               + "RNaseH" + "\t"
                               + "INT" + "\t"
                               + "ENV" + "\n")


        for line in domain_duplicate_removal:

            if not line.startswith("##"):

                line = line.strip()
                line_list = line.split("\t")

                seq_name = line_list[0]
                strand = seq_name.split("|")[-1]

                flage = 0

                line_record_list = []
                if strand == "same":
                    start_site = int(seq_name.split("|")[1])
                    # end_site =int(seq_name.split("|")[2])
                    for n in range(1, len(line_list)):
                        if line_list[n] != "not found" and line_list[n].find("relatively short") == -1:
                            xx = (line_list[n].split(",")[-1]).split("(")[0]
                            if xx == "":
                                sites = (line_list[n].split("(")[1]).split(")")[
                                    0]
                                domain_start = int(
                                    sites.split("..")[0]) - start_site + 1
                                domain_end = int(
                                    sites.split("..")[-1]) - start_site + 1
                                line_record_list.append(str(domain_start) + ".."
                                                        + str(domain_end))
                            elif xx == "complement":
                                sites = (line_list[n].split("(")[1]).split(")")[
                                    0]
                                domain_start = int(
                                    sites.split("..")[0]) - start_site + 1
                                domain_end = int(
                                    sites.split("..")[-1]) - start_site + 1
                                line_record_list.append("complement" + "("
                                                        + str(
                                    domain_start) + ".."
                                                        + str(domain_end) + ")")
                                flage = 1
                        else:
                            line_record_list.append(line_list[n])


                elif strand == "complement":
                    start_site = int(seq_name.split("|")[2])
                    for n in range(1, len(line_list)):

                        if line_list[n] != "not found" and line_list[n].find("relatively short") == -1:

                            xx = (line_list[n].split(",")[-1]).split("(")[0]
                            if xx == "complement":
                                sites = (line_list[n].split("(")[1]).split(")")[
                                    0]
                                domain_start = start_site - int(
                                    sites.split("..")[1]) + 1
                                domain_end = start_site - int(
                                    sites.split("..")[0]) + 1

                                line_record_list.append(str(domain_start) + ".."
                                                        + str(domain_end))
                            elif xx == "":
                                sites = (line_list[n].split("(")[1]).split(")")[
                                    0]
                                domain_start = start_site - int(
                                    sites.split("..")[1]) + 1
                                domain_end = start_site - int(
                                    sites.split("..")[0]) + 1

                                line_record_list.append("complement" + "("
                                                        + str(
                                    domain_start) + ".."
                                                        + str(domain_end) + ")")
                                flage = 1

                        else:
                            line_record_list.append(line_list[n])

                new_line = "\t".join(line_record_list)

                out_domain_annos.write(seq_name + "\t"
                    + new_line + "\n")

                if flage == 1:
                    out_domain_annos.write(
                        "##Warn: the annotation and chain direction of '" +
                        seq_name
                        + "' needs to be checked in file of 'domain_annotation_in_genome_final.txt'!" + "\n")


        domain_duplicate_removal.close()
        out_domain_annos.close()


