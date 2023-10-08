# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/10 20:36

"""

from Bio.Seq import Seq
from Bio import SeqIO
from public_function import make_dir


class GetSeqForHmmer(object):
    def __init__(self, blast_flank_seq, simplify_gff3,
                 duplicate_removal_gff3,
                 seq_with_ltrs_nt, seq_with_ltrs_aa,
                 flank_without_ltrs_nt, flank_without_ltrs_aa):
        """
        Make sequence for HMMER analysis
        :param blast_flank_seq: flank seq from balst
        :param simplify_gff3: simplified gff3 file
        :param duplicate_removal_gff3: simplified gff3 with duplicate removal
        """
        super(GetSeqForHmmer, self).__init__()

        self.blast_flank_seq = blast_flank_seq
        self.simplify_gff3_path = simplify_gff3
        self.duplicate_removal_path = duplicate_removal_gff3

        self.seq_with_ltrs_nt = seq_with_ltrs_nt
        self.seq_with_ltrs_aa = seq_with_ltrs_aa

        self.flank_without_ltrs_nt = flank_without_ltrs_nt
        self.flank_without_ltrs_aa = flank_without_ltrs_aa

    def run(self):


        flank_seq_dict = SeqIO.to_dict(SeqIO.parse(self.blast_flank_seq, "fasta"))

        duplicate_removal_file = open(self.duplicate_removal_path,
                                      "r",encoding="utf-8")
        

        # 1. Extracted ERVs with paired LTR
        pairLTR_eves_file = open(self.seq_with_ltrs_nt, "w",
                                 encoding="utf-8")
        pairLTR_evesAA_file = open(self.seq_with_ltrs_aa, "w",
                                   encoding="utf-8")

        for line in duplicate_removal_file:
            if not line.startswith("#"):
                ltr_record_list = line.strip().split("\t")
                ltr1_start = int(ltr_record_list[3].split("..")[0])
                ltr2_end = int(ltr_record_list[4].split("..")[1])

                qury_erv_start = ltr1_start
                qury_erv_end = ltr2_end
        
                erv_seq_contain = Seq(
                    str(flank_seq_dict[ltr_record_list[0]][
                        qury_erv_start - 1: qury_erv_end].seq))

                # the full-length nucleotide sequence of pairLTR ervs
                pairLTR_eves_file.write(">" + ltr_record_list[-3] + "\n"
                                        + str(erv_seq_contain) + "\n")
        
                # the aa sequence(6 reading frames) of pairLTR ervs
        
                pairLTR_evesAA_file.write(">" + ltr_record_list[-3]
                                          + "#ORF1" + "\n"
                                          + str(
                    erv_seq_contain.translate(table=1)) + "\n")
        
                pairLTR_evesAA_file.write(">" + ltr_record_list[-3]
                                          + "#ORF2" + "\n"
                                          + str(
                    erv_seq_contain[1:].translate(table=1))
                                          + "\n")
        
                pairLTR_evesAA_file.write(">" + ltr_record_list[-3]
                                          + "#ORF3" + "\n"
                                          + str(
                    erv_seq_contain[2:].translate(table=1))
                                          + "\n")
        
                pairLTR_evesAA_file.write(">" + ltr_record_list[-3]
                                          + "#ORF4" + "\n"
                                          + str(
                    erv_seq_contain.reverse_complement().translate(table=1))
                                          + "\n")
        
                pairLTR_evesAA_file.write(">" + ltr_record_list[-3]
                                          + "#ORF5" + "\n"
                                          + str(
                    erv_seq_contain.reverse_complement()[1:].translate(table=1))
                                          + "\n")
        
                pairLTR_evesAA_file.write(">" + ltr_record_list[-3]
                                          + "#ORF6" + "\n"
                                          + str(
                    erv_seq_contain.reverse_complement()[2:].translate(table=1))
                                          + "\n")
        
        duplicate_removal_file.close()
        pairLTR_eves_file.close()
        pairLTR_evesAA_file.close()


        # 2. Extract flank sequences without paired LTR
        # using the harvest_gff3_simplify.txt
        gff3_file = open(self.simplify_gff3_path, "r",encoding="utf-8")
        
        no_pairLTR_seq = open(self.flank_without_ltrs_nt, "w",encoding="utf-8")
        no_pairLTR_seqAA = open(self.flank_without_ltrs_aa, "w",encoding="utf-8")
        
        ltr_envs_record = set()
        
        for line in gff3_file:
            ltr_envs_record.add(line.split("\t")[0].strip())
        
        for key in flank_seq_dict:
            if key.strip() not in ltr_envs_record:
                seq_contain = Seq(str(flank_seq_dict[key].seq))
        
                no_pairLTR_seq.write(">" + key + "\n" + str(seq_contain) + "\n")
        

                no_pairLTR_seqAA.write(">" + key + "#ORF1" + "\n"
                                       + str(seq_contain.translate(table=1))
                                       + "\n")
        
                no_pairLTR_seqAA.write(">" + key + "#ORF2" + "\n"
                                       + str(seq_contain[1:].translate(table=1))
                                       + "\n")
        
                no_pairLTR_seqAA.write(">" + key + "#ORF3" + "\n"
                                       + str(seq_contain[2:].translate(table=1))
                                       + "\n")
        
                no_pairLTR_seqAA.write(">" + key + "#ORF4" + "\n"
                                       + str(
                    seq_contain.reverse_complement().translate(table=1))
                                       + "\n")
        
                no_pairLTR_seqAA.write(">" + key + "#ORF5" + "\n"
                                       + str(
                    seq_contain.reverse_complement()[1:].translate(table=1))
                                       + "\n")
        
                no_pairLTR_seqAA.write(">" + key + "#ORF6" + "\n"
                                       + str(
                    seq_contain.reverse_complement()[2:].translate(table=1))
                                       + "\n")
        
        no_pairLTR_seq.close()
        no_pairLTR_seqAA.close()



def GetSeqForHmmer2(input_file, outnt, outAA):
    flank_seq_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    no_pairLTR_seq = open(outnt, "w", encoding="utf-8")
    no_pairLTR_seqAA = open(outAA, "w", encoding="utf-8")

    for key in flank_seq_dict:
        seq_contain = Seq(str(flank_seq_dict[key].seq))

        no_pairLTR_seq.write(">" + key + "\n" + str(seq_contain) + "\n")

        no_pairLTR_seqAA.write(">" + key + "#ORF1" + "\n"
                               + str(seq_contain.translate(table=1))
                               + "\n")

        no_pairLTR_seqAA.write(">" + key + "#ORF2" + "\n"
                               + str(seq_contain[1:].translate(table=1))
                               + "\n")

        no_pairLTR_seqAA.write(">" + key + "#ORF3" + "\n"
                               + str(seq_contain[2:].translate(table=1))
                               + "\n")

        no_pairLTR_seqAA.write(">" + key + "#ORF4" + "\n"
                               + str(
            seq_contain.reverse_complement().translate(table=1))
                               + "\n")

        no_pairLTR_seqAA.write(">" + key + "#ORF5" + "\n"
                               + str(
            seq_contain.reverse_complement()[1:].translate(table=1))
                               + "\n")

        no_pairLTR_seqAA.write(">" + key + "#ORF6" + "\n"
                               + str(
            seq_contain.reverse_complement()[2:].translate(table=1))
                               + "\n")

    no_pairLTR_seq.close()
    no_pairLTR_seqAA.close()