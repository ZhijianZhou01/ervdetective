# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/11 17:47

"""
from Bio import SeqIO
from public_function import reversecomp_seq


class GetPairLTRervs(object):
    def __init__(self, domain_map_genome, host_path, doamin_seq_dir,
                  erv_seq_path, pairLTRs_path):
        """
        Extract pairLTR ERVs and domains
        :param domain_map_genome:
        :param host_path:
        :param doamin_seq_dir:
        :param erv_seq_path:
        :param pairLTRs_path:
        """
        super(GetPairLTRervs,self).__init__()

        self.domain_map_genome = domain_map_genome
        self.doamin_seq_dir = doamin_seq_dir
        self.host_path = host_path
        self.erv_seq_path = erv_seq_path
        self.pairLTRs_path = pairLTRs_path
    
    def run(self):

        domain_in_genome_file = open(self.domain_map_genome,"r",encoding="utf-8")

        host_seq_dict = SeqIO.to_dict(SeqIO.parse(self.host_path, "fasta")) # 宿主基因组

        erv_seq = open(self.erv_seq_path,"w",encoding="utf-8")

        out_pair_LTRs = open(self.pairLTRs_path,"w",encoding="utf-8")
        
        for line in domain_in_genome_file:
            if not line.startswith("##"):
        
                line_list = line.split("\t")
                seq_name = line_list[0]
                chain_strand = seq_name.split("|")[-1]
        
                erv_start = int(seq_name.split("|")[1])
                erv_end = int(seq_name.split("|")[2])
        
                key = seq_name.split("|")[0]
        
                five_ltr = line_list[1]
                three_ltr = line_list[2]
        
                if chain_strand == "same":
        
                    ervs_seq_contain = str(host_seq_dict[key].seq[erv_start-1 :erv_end])
        
                    five_ltr_start= int(five_ltr.split("..")[0])
                    five_ltr_end = int(five_ltr.split("..")[1])
                    five_ltr_seq = str(host_seq_dict[key].seq[five_ltr_start - 1: five_ltr_end])
        
                    three_ltr_start = int(three_ltr.split("..")[0])
                    three_ltr_end = int(three_ltr.split("..")[1])
                    three_ltr_seq = str(host_seq_dict[key].seq[three_ltr_start - 1: three_ltr_end])
        
        
        
                else:
                    ervs_seq_contain = reversecomp_seq(str(host_seq_dict[key].seq[erv_start - 1:erv_end]))
        
                    five_ltr_site = (five_ltr.split("(")[-1]).split(")")[0]
                    five_ltr_start = int(five_ltr_site.split("..")[0])
                    five_ltr_end = int(five_ltr_site.split("..")[1])
                    five_ltr_seq = reversecomp_seq(str(
                        host_seq_dict[key].seq[five_ltr_start - 1: five_ltr_end]))
        
                    three_ltr_site = (three_ltr.split("(")[-1]).split(")")[0]
                    three_ltr_start = int(three_ltr_site.split("..")[0])
                    three_ltr_end = int(three_ltr_site.split("..")[1])
                    three_ltr_seq = reversecomp_seq(str(
                        host_seq_dict[key].seq[three_ltr_start - 1: three_ltr_end]))
        
        
                erv_seq.write(">" + seq_name + "\n" + ervs_seq_contain + "\n")
                out_pair_LTRs.write(">" + seq_name + "_five_ltr" + "\n"
                                    + five_ltr_seq + "\n"
                                    + ">" + seq_name + "_three_ltr" + "\n"
                                    + three_ltr_seq + "\n")


                for n in range(4,len(line_list)):
                    domain_record = line_list[n].strip()
                    if domain_record != "not found" and domain_record.find("relatively short") == -1:
                        domain_name = (domain_record.split(",")[0]).split("_")[0]
        
                        domain_range = (domain_record.split("(")[-1]).strip(")")
        
                        domain_start = int(domain_range.split("..")[0])
                        domain_end = int(domain_range.split("..")[1])
        
                        if (domain_record.split(",")[-1]).split("(")[0]=="complement":
                            domain_seq = reversecomp_seq(str(host_seq_dict[key].seq[domain_start-1 :domain_end]))
        
        
                        else:
                            domain_seq = str(host_seq_dict[key].seq[domain_start-1 :domain_end])
        
        
                        out_domain = open(self.doamin_seq_dir + r"/" + domain_name + ".fasta", "a",
                            encoding="utf-8")
        
                        out_domain.write(">" + seq_name + "_" + domain_name + "\n"
                              +  domain_seq + "\n")

        domain_in_genome_file.close()
        erv_seq.close()
        out_pair_LTRs.close()


class GetPotentialEves(object):
    def __init__(self, map_duplicate_removal, host_path,
                 doamin_seq_dir, erv_seq_path):
        """
        Extract potential ERVs and domains on the flank sequence
        :param map_duplicate_removal:
        :param host_path:
        :param doamin_seq_dir:
        :param erv_seq_path:
        """
        super(GetPotentialEves, self).__init__()

        self.map_duplicate_removal = map_duplicate_removal
        self.host_path = host_path
        self.doamin_seq_dir = doamin_seq_dir
        self.erv_seq_path = erv_seq_path

    def run(self):
        domain_duplicate_removal = open(self.map_duplicate_removal, "r",
                                        encoding="utf-8")

        host_seq_dict = SeqIO.to_dict(SeqIO.parse(self.host_path, "fasta"))

        erv_seq = open(self.erv_seq_path, "w",encoding="utf-8")

        for line in domain_duplicate_removal:
            if not line.startswith("##"):
                line_list = line.split("\t")
                seq_name = line_list[0]
                chain_strand = seq_name.split("|")[-1]

                erv_start = int(seq_name.split("|")[1])
                erv_end = int(seq_name.split("|")[2])

                key = seq_name.split("|")[0]

                if chain_strand == "same":

                    ervs_seq_contain = str(
                        host_seq_dict[key].seq[erv_start - 1:erv_end])

                else:
                    ervs_seq_contain = reversecomp_seq(
                        str(host_seq_dict[key].seq[erv_start - 1:erv_end]))

                erv_seq.write(">" + seq_name + "\n" + ervs_seq_contain + "\n")


                for n in range(1, len(line_list)):
                    domain_record = line_list[n].strip()
                    if domain_record != "not found" and domain_record.find("relatively short") == -1:

                        domain_name = (domain_record.split(",")[0]).split("_")[
                            0]

                        domain_range = (domain_record.split("(")[-1]).strip(")")

                        domain_start = int(domain_range.split("..")[0])
                        domain_end = int(domain_range.split("..")[1])

                        if (domain_record.split(",")[-1]).split("(")[
                            0] == "complement":
                            domain_seq = reversecomp_seq(str(
                                host_seq_dict[key].seq[
                                domain_start - 1:domain_end]))
                        else:
                            domain_seq = str(host_seq_dict[key].seq[
                                             domain_start - 1:domain_end])

                        out_domain = open(self.doamin_seq_dir + r"/"
                                          + domain_name + ".fasta", "a",
                                          encoding="utf-8")

                        out_domain.write(
                            ">" + seq_name + "_" + domain_name + "\n"
                            + domain_seq + "\n")

        domain_duplicate_removal.close()
        erv_seq.close()


