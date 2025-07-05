# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2021/11/9 15:19

"""

import sys
import os
import platform
# import psutil

app_dir = os.path.dirname(os.path.realpath(__file__))
# print(app_dir)


if platform.system().lower() == "windows":
    app_dir = app_dir.replace("\\", "/")

# sys.path.append(app_dir) # has been in __init__


import time
import subprocess
import zipfile
import argparse

from datetime import datetime
from get_flank import GetFlanks
from simplify_gff3 import SimplyGFF
from seqs_for_hmmer import GetSeqForHmmer,GetSeqForHmmer2
from simplify_hmmer import SimplyHMMER
from map_genome import PairLTRervMap, NoPairLTReveMap
from get_ervs import GetPairLTRervs, GetPotentialEves
from annotion_ervs import PairLTRErvAnno, PotentialEveAnno
from public_function import make_dir,hmmer_contain,get_all_path



#  ***** External program *****


makeblastdb_execute = "makeblastdb"     # from blast
tblastn_execute = "tblastn"             # from blast
gt_execute = "gt"                       # from genometools
hmmpress_execute = "hmmpress"           # from hmmer
hmmscan_execute = "hmmscan"             # from hmmer



if platform.system().lower() == "windows":
    makeblastdb_execute = "makeblastdb.exe"
    tblastn_execute = "tblastn.exe"
    gt_execute = "gt.exe"
    hmmpress_execute = "hmmpress.exe"
    hmmscan_execute = "hmmscan.exe"




#  ***** Built-in database *****

ERVprofiles_zip = app_dir + "/Internal_database/ERVprofiles/ERVprofile.zip"
ERVprofiles_path = app_dir + "/Internal_database/ERVprofiles/ERV.profile"
probe_dir =  app_dir + "/Internal_database/probe"



example_use = r'''
----------------☆ Example of use ☆-----------------

ervdetective -i myotis_lucifugus.fna -p myotis_lucifugus -n 10 -o output
  
----------------------☆  End  ☆---------------------

'''

"""
def calculate_memory():
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    memory = info.uss / 1024 / 1024
    return memory
"""


# ###  *** Step 0, Processing input parameters  **

def parameter():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        prog="ervdetective",
        description="",
        epilog=example_use)

    parser.add_argument(
        "-i", dest="host",
        help="The file-path of host genome sequence, the suffix is generally *.fna, *.fas, *.fasta.",
        default="")

    parser.add_argument(
        "-n", dest="thread",
        help="Specify the number of threads used, default: 1.",
        type=int,
        default=1)

    parser.add_argument(
        "-p", dest="prefix",
        help="The prefix of output file, default character: 'host'.",
        type=str,
        default="host")

    parser.add_argument(
        "-o", dest="output",
        help="The path of output folder to store all the results.",
        default="")

    parser.add_argument(
        "-eb", dest="eblast",
        help="Specify threshold of e-value for BLAST search, default: 1e-5.",
        type=float,
        default=0.00001)

    parser.add_argument(
        "-f", dest="flank",
        help="The length of extended flank sequence on either side of "
             "the blast hit-site, default: 15000.",
        type=int,
        default=15000)

    parser.add_argument(
        "-l1", dest="minltr",
        help="Specify minimum length of LTR, default: 100.",
        type=int,
        default=100)

    parser.add_argument(
        "-l2", dest="maxltr",
        help="Specify maximum length of LTR, default: 1000.",
        type=int,
        default=1000)

    parser.add_argument(
        "-s", dest="ltrsimilar",
        help="Specify threshold(%%) of the similarity of paired LTRs, default: 80.",
        type=float,
        default=80)

    parser.add_argument(
        "-d1", dest="mindistltr",
        help="The minimum interval of paired-LTRs start-positions, default: 1000.",
        type=int,
        default=1000)

    parser.add_argument(
        "-d2", dest="maxdistltr",
        help="The maximum interval of paired-LTRs start-positions, default: 15000.",
        type=int,
        default=15000)

    parser.add_argument(
        "-t1", dest="mintsd",
        help="The minimum length for each TSD site, default: 4.",
        type=int,
        default=4)

    parser.add_argument(
        "-t2", dest="maxtsd",
        help="The maximum length for each TSD site, default: 6.",
        type=int,
        default=6)

    parser.add_argument(
        "-motif", dest="motif",
        help="Specify start-motif (2 nucleotides) and end-motif (2 nucleotides), default string: TGCA.",
        type=str,
        default="TGCA")

    parser.add_argument(
        "-mis", dest="mismotif",
        help="The maximum number of mismatches nucleotides in motif, default: 1.",
        type=int,
        default=1)

    parser.add_argument(
        "-ed", dest="ehmmer",
        help="Specify threshold of e-value using for HMMER search, default: 1e-6.",
        type=float,
        default=0.000001)

    parser.add_argument(
        "--gag", dest="GAG_length",
        help="The threshold of length of GAG protein in HMMER search, default: 250 aa.",
        type=int,
        default=250)

    parser.add_argument(
        "--pro", dest="PRO_length",
        help="The threshold of length of PRO protein in HMMER search, default: 50 aa.",
        type=int,
        default=50)

    parser.add_argument(
        "--rt", dest="RT_length",
        help="The threshold of length of RT protein in HMMER search, default: 150 aa.",
        type=int,
        default=150)

    parser.add_argument(
        "--rh", dest="RNaseH_length",
        help="The threshold of length of RNaseH protein in HMMER search, default: 65 aa.",
        type=int,
        default=65)

    parser.add_argument(
        "--int", dest="INT_length",
        help="The threshold of length of INT protein in HMMER search, default: 150 aa.",
        type=int,
        default=150)

    parser.add_argument(
        "--env", dest="ENV_length",
        help="The threshold of length of ENV protein in HMMER search, default: 250 aa.",
        type=int,
        default=250)



    theargs = parser.parse_args(sys.argv[1:])

    return theargs


def starts():

    print("\n" + "-------------------------------------------------")

    print("  Name: Endogenous retroviruses detective (ERVdetective)")

    print("  Description: An efficient pipeline for identification and annotation of endogenous retroviruses.")

    print("  Version: 1.0.9 (2024-05-01)")

    print("  Author: Zhi-Jian Zhou")

    print("-------------------------------------------------" + "\n")


    myargs = parameter()
    
    print(myargs)
    
    print("\n")

    time_start = datetime.today().now()

    # start_memory = calculate_memory()

    parameter_dic = {}

    parameter_dic["host_path"] = myargs.host

    parameter_dic["out_prefix"] = myargs.prefix

    parameter_dic["out_dir"] = myargs.output

    parameter_dic["thread"] = str(myargs.thread)

    parameter_dic["blastE"] = str(myargs.eblast)

    parameter_dic["flanklen"] = myargs.flank

    parameter_dic["minlength"] = str(myargs.minltr)

    parameter_dic["maxlength"] = str(myargs.maxltr)

    parameter_dic["LTRsimilar"] = str(myargs.ltrsimilar)

    parameter_dic["mindistltr"] = str(myargs.mindistltr)

    parameter_dic["maxdistltr"] = str(myargs.maxdistltr)

    parameter_dic["mintsd"] = str(myargs.mintsd)

    parameter_dic["maxtsd"] = str(myargs.maxtsd)

    parameter_dic["motif_str"] = myargs.motif

    parameter_dic["motifmis"] = str(myargs.mismotif)

    parameter_dic["HmmerE"] = str(myargs.ehmmer)

    parameter_dic["GAG_length"] = int(myargs.GAG_length)

    parameter_dic["PRO_length"] = int(myargs.PRO_length)

    parameter_dic["RT_length"] = int(myargs.RT_length)

    parameter_dic["RNaseH_length"] = int(myargs.RNaseH_length)

    parameter_dic["INT_length"] = int(myargs.INT_length)

    parameter_dic["ENV_length"] = int(myargs.ENV_length)



    if platform.system().lower() == "windows":
        parameter_dic["host_path"] = parameter_dic["host_path"].replace("\\", "/")
        parameter_dic["out_dir"] = parameter_dic["out_dir"].replace("\\", "/")


    run_log_dir = parameter_dic["out_dir"] + "/Run_log"

    make_dir(parameter_dic["out_dir"])
    make_dir(run_log_dir)


    # When the input and output are not right
    if parameter_dic["host_path"] == "":
        print("Error: the path of input host-genome was not specified! Please use "
              "'ervdetective -h' to view help document.")
        sys.exit()

    if parameter_dic["out_dir"] == "":
        print("Error: the path of out folder was not specified! Please use "
              "'ervdetective -h' to view help document.")
        sys.exit()


    ###  *** Step 1, run tblastn ***

    print(">>> Blast search for ERVs loci....")

    blast_result_path = run_log_dir + "/Record1_blast_result"

    host_db_dir = blast_result_path + "/seqdb"

    blast_out_dir = blast_result_path + "/blast_result"

    make_dir(blast_result_path)

    make_dir(host_db_dir)

    make_dir(blast_out_dir)


    makeblastdb_command = (makeblastdb_execute + " -in "
                            + parameter_dic["host_path"]
                            + " -dbtype nucl "
                            + "-out " + host_db_dir + "/"
                            + parameter_dic["out_prefix"]
                            + ".")


    makeblastdb_process = subprocess.Popen(makeblastdb_command,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT,
                                           universal_newlines=True,
                                           shell=True)

    while True:

        output = makeblastdb_process.stdout.readline()
        if makeblastdb_process.poll() is not None:
            break

        elif output:
            print(output)

    makeblastdb_process.terminate()


    # blast process
    probe_file_list = get_all_path(probe_dir)

    tblastn_commd_list = []
    for each_path in probe_file_list:
        each_path = each_path.strip()
        file_name = str(each_path.split("/")[-1]).strip()
        print("\n" + ">>> Use the " +  file_name + " as a probe...")
        blast_result = (blast_out_dir + "/"
                        + parameter_dic["out_prefix"]
                        + "_" + file_name + ".txt")

        each_tblastn_commd = (tblastn_execute + " -query " + each_path
                              + " -db " + host_db_dir + "/"
                              + parameter_dic["out_prefix"] + "."
                              + " -out " + blast_result
                              + " -outfmt 6"
                              + " -evalue " + parameter_dic["blastE"]
                              + " -num_threads " + parameter_dic["thread"])

        tblastn_commd_list.append(each_tblastn_commd)


    tblastn_commd = " && ".join(tblastn_commd_list)

    tblastn_process = subprocess.Popen(tblastn_commd,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT,
                                           universal_newlines=True,
                                           shell=True)

    while True:

        output = tblastn_process.stdout.readline()
        if tblastn_process.poll() is not None:
            break

        elif output:
            print(output)

    tblastn_process.terminate()


    # Determine whether the blast result file is empty

    blast_file_list = get_all_path(blast_out_dir)

    blast_file_flage = "False"
    for each_file in blast_file_list:
        file_size = os.path.getsize(each_file)
        if file_size > 0:
            blast_file_flage = "True"


    if blast_file_flage == "False":
        print("\n"
              + "*** Note: no homologues of ERVs were found in the input sequence!"
              + "\n")

        sys.exit()


    ###  ***Step 2, get flank sequence on both sides of these hit sites***

    flank_seq_dir = run_log_dir + "/Record2_flank_seqs"

    make_dir(flank_seq_dir)

    flank_seq_file = flank_seq_dir + "/blast_flank_seq.fasta"

    flank_get = GetFlanks(parameter_dic["host_path"], blast_out_dir,
                          parameter_dic["flanklen"], flank_seq_file)
    flank_get.run()



    ###  ***Step 3, identify LTR by LTRharvest***

    print("\n" + ">>> Search pair-LTR....")

    identy_pairLTR = run_log_dir + "/Record3_identy_pairedLTR"

    flankseq_index_dir = identy_pairLTR + "/flank_seq_index"

    ltrharvest_result_dir = identy_pairLTR + "/LTRharvest_result"

    flankseq_index_file = flankseq_index_dir + "/FlankSeq"

    ltrharvest_gff3_file = ltrharvest_result_dir + "/flank_seq_harvest.gff3"
    ltrharvest_out_seq = ltrharvest_result_dir + "/LTRharvest_out.fas"
    harvest_log_file = ltrharvest_result_dir + "/harvest_log.txt"

    make_dir(identy_pairLTR)
    make_dir(flankseq_index_dir)
    make_dir(ltrharvest_result_dir)

    make_index_commd = (gt_execute + " suffixerator"
                          + " -db " + flank_seq_file
                          + " -indexname " + flankseq_index_file
                          + " -tis -suf -lcp -des -ssp -sds -dna")

    ltrharvest_commd = (gt_execute + " ltrharvest"
                      + " -index " + flankseq_index_file
                      + " -minlenltr " + parameter_dic["minlength"]
                      + " -maxlenltr " + parameter_dic["maxlength"]
                      + " -similar " + parameter_dic["LTRsimilar"]
                      + " -seed 20"
                      + " -mindistltr " + parameter_dic["mindistltr"]
                      + " -maxdistltr " + parameter_dic["maxdistltr"]
                      + " -mintsd " + parameter_dic["mintsd"]
                      + " -maxtsd " + parameter_dic["maxtsd"]
                      + " -motif " + parameter_dic["motif_str"]
                      + " -motifmis " + parameter_dic["motifmis"]
                      + " -gff3 " + ltrharvest_gff3_file
                      + " -out " + ltrharvest_out_seq
                      + " > " + harvest_log_file)



    find_ltr_commd = make_index_commd + " && " + ltrharvest_commd


    find_ltr_process = subprocess.Popen(find_ltr_commd, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               universal_newlines=True, shell=True)

    while True:

        output = find_ltr_process.stdout.readline()
        if find_ltr_process.poll() is not None:
            break

        elif output:
            print(output)

    find_ltr_process.terminate()



    ###  ***Step4, Prepare the sequence for HMMER ***

    print("\n" + ">>> Prepare the sequence for HMMER...")

    seqs_for_hmmer_dir = run_log_dir + "/Record4_get_Seqs"

    make_dir(seqs_for_hmmer_dir)

    seq_with_ltrs = seqs_for_hmmer_dir + "/seq_with_paired_ltrs"

    make_dir(seq_with_ltrs)

    seq_with_ltrs_nt = (seq_with_ltrs
                                 + "/seq_with_paired_ltrs.fas")

    seq_with_ltrs_aa = (seq_with_ltrs
                                 + "/seq_with_paired_ltrs_aa.fas")

    flank_without_ltrs = seqs_for_hmmer_dir + "/flank_without_paired_ltrs"

    make_dir(flank_without_ltrs)

    flank_without_ltrs_nt = (flank_without_ltrs
                              + "/flank_without_ltrs.fasta")

    flank_without_ltrs_aa = (flank_without_ltrs
                              + "/flank_without_ltrs_AA.fasta")

    simplify_gff3_dir = identy_pairLTR + "/Simplify_gff3"
    simplify_gff3_path = simplify_gff3_dir + "/harvest_gff3_simplify.txt"

    gff_duplicate_removal_path = (simplify_gff3_dir
                                  + "/harvest_gff3_duplicate_removal.txt")

    LTRharvest_out_size = os.path.getsize(ltrharvest_out_seq)


    if LTRharvest_out_size != 0:

        make_dir(simplify_gff3_dir)

        simplify_gff3_process = SimplyGFF(ltrharvest_gff3_file,
                                           simplify_gff3_path,
                                          gff_duplicate_removal_path,
                                          parameter_dic["mintsd"])
        simplify_gff3_process.run()

        # get sequences for HMMER
        seq_for_hmmer = GetSeqForHmmer(flank_seq_file, simplify_gff3_path,
                                       gff_duplicate_removal_path,
                                       seq_with_ltrs_nt,seq_with_ltrs_aa,
                                       flank_without_ltrs_nt,flank_without_ltrs_aa)
        seq_for_hmmer.run()

    else:
        print("\n" +
            "*** Note: no pair-LTR ERV was found, and the HMMER search domain was started!")

        make_dir(flank_without_ltrs)

        GetSeqForHmmer2(flank_seq_file, flank_without_ltrs_nt,
                        flank_without_ltrs_aa)


    ###  ***Step 5 annotation domain***

    # init ERV.profile database

    if os.path.exists(ERVprofiles_path) == False:
        with zipfile.ZipFile(ERVprofiles_zip, "r") as zip_file:
            zip_file.extractall(os.path.dirname(ERVprofiles_zip))

    h3p_path = ERVprofiles_path + ".h3p"

    if os.path.exists(h3p_path) == False:
        hmmpress_commd = (hmmpress_execute + " " + ERVprofiles_path)

        hmmpress_process = subprocess.Popen(hmmpress_commd,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT,
                                            universal_newlines=True,
                                            shell=True)

        while True:

            output = hmmpress_process.stdout.readline()
            if hmmpress_process.poll() is not None:
                break

            elif output:
                print(output)

        hmmpress_process.terminate()



    ##  **Step 5.1 annotation of ervs with paired-LTRs (seq_with_ltrs_aa) **

    pairltr_ervs_map_genome = ""

    if LTRharvest_out_size != 0: 

        hmmer_result_dir = run_log_dir + "/Record5_HMMER_analysis"
        make_dir(hmmer_result_dir)

        # * HMMER analysis *

        pairLTR_ERVs_hmmer = hmmer_result_dir + "/seqs_with_paired_LTRs"

        make_dir(pairLTR_ERVs_hmmer)

        seqs_paired_LTRs_hmmer = (pairLTR_ERVs_hmmer + "/"
                                 + parameter_dic["out_prefix"]
                                 + "_seq_with_paired_LTR_hmmer.txt")

        seqs_paired_LTRs_hmmer_log = (pairLTR_ERVs_hmmer + "/"
                                  + parameter_dic["out_prefix"]
                                  + "_seq_with_paired_LTR_hmmer_log.txt")


        hmmer_commd1 = (hmmscan_execute + " --domtblout "
                        + seqs_paired_LTRs_hmmer
                        + " -E " + parameter_dic["HmmerE"]
                        + " --cpu " + parameter_dic["thread"] + " "
                        + ERVprofiles_path + " "
                        + seq_with_ltrs_aa
                        + " > " + seqs_paired_LTRs_hmmer_log)

        print("\n" + ">>> Start HMMER search for ERVs with paired-LTRs...")

        hmmer1_process = subprocess.Popen(hmmer_commd1,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT,
                                          universal_newlines=True,
                                          shell=True)

        while True:

            output = hmmer1_process.stdout.readline()
            if hmmer1_process.poll() is not None:
                break

            elif output:
                print(output)

        hmmer1_process.terminate()

        if hmmer_contain(seqs_paired_LTRs_hmmer) == "True":
            print("\n" + ">>> Annotate domain...." + "\n")
            #  * simplify the result of hmmer *
            seqs_paired_LTRs_hmmer_simplify = (pairLTR_ERVs_hmmer + "/"
                                           + parameter_dic["out_prefix"]
                                           + "_seq_with_paired_LTR_hmmer_simplify.txt")

            hmmer_head_pairedLTR = (
                        "##seq name for HMMER" + "\t" + "GAG domain"
                        + "\t" + "DUT domain" + "\t"
                        + "AP domain" + "\t" + "RT domain"
                        + "\t" + "RNaseH domain" + "\t"
                        + "INT domain" + "\t" + "ENV domain"
                        + "\n")

            simp_hmmer1 = SimplyHMMER(seqs_paired_LTRs_hmmer,
                                      seqs_paired_LTRs_hmmer_simplify,
                                      hmmer_head_pairedLTR,parameter_dic)
            simp_hmmer1.run()

            # * domian map to genome *
            annotation_dir = (parameter_dic["out_dir"]
                              + "/Annotatio_and_extract_domain")

            make_dir(annotation_dir)
            pairLTR_ervs_annotation_dir = annotation_dir + "/ERVs_with_paired_LTRs"
            make_dir(pairLTR_ervs_annotation_dir)

            pairltr_ervs_map_genome = (pairLTR_ervs_annotation_dir
                                       + "/" + parameter_dic["out_prefix"]
                                       + "_domain_annotation_in_genome_final.txt")

            domain_in_genome_map1 = PairLTRervMap(seqs_paired_LTRs_hmmer_simplify,
                                                  gff_duplicate_removal_path,
                                                  pairltr_ervs_map_genome)
            domain_in_genome_map1.run()

            #  *  extract the domain of ervs with paired-LTRs  *

            doamin_seq_dir1 = (pairLTR_ervs_annotation_dir + "/"
                               + parameter_dic["out_prefix"]
                               + "_domain_seqs")
            make_dir(doamin_seq_dir1)

            pairLTRerv_seq_path = (pairLTR_ervs_annotation_dir + "/"
                                   + parameter_dic["out_prefix"]
                                   + "_extracted_ERVs_with_paired_LTRs.fas")

            pairLTRs_path = (pairLTR_ervs_annotation_dir + "/"
                             + parameter_dic["out_prefix"]
                             + "_extracted_pairedLTR.fas")

            getpairLTRervs_process = GetPairLTRervs(pairltr_ervs_map_genome,
                                                    parameter_dic["host_path"],
                                                    doamin_seq_dir1,
                                                    pairLTRerv_seq_path,
                                                    pairLTRs_path)

            getpairLTRervs_process.run()

            #  * annotate domain  on the extracted pair of LTR_ervs *

            pairltr_ervs_anno_path = (pairLTR_ervs_annotation_dir
                                 + "/domain_annotation_in_ERVs_with_paired_LTRs.txt")

            pairltr_ervs_annotation = PairLTRErvAnno(pairltr_ervs_map_genome,
                                                     pairltr_ervs_anno_path,
                                                     parameter_dic["out_prefix"])

            pairltr_ervs_annotation.run()

            time_end = time.time()

            print("\n" + ">>> The result was saved in: "
                  + annotation_dir + "\n")

        else:
            print("\n" + "*** Note: no domain was found in the identified potential ERVs with paired-LTRs!")



    ##  ** 5.2 annotation of potential ERVLEs whitout paired-LTRs
    no_pairLTR_flankseqaa_size = os.path.getsize(flank_without_ltrs_aa)

    if no_pairLTR_flankseqaa_size != 0:

        hmmer_result_dir = run_log_dir + "/Record5_HMMER_analysis"
        make_dir(hmmer_result_dir)

        # * HMMER analysis *
        no_LTR_flank_hmmer = hmmer_result_dir + "/flank_without_paired_LTRs"

        make_dir(no_LTR_flank_hmmer)

        no_LTR_flank_hmmer_su = (no_LTR_flank_hmmer + "/"
                                 + parameter_dic["out_prefix"]
                                 + "_flank_without_paired_LTR_hmmer.txt")

        no_LTR_flank_hmmer_log = (no_LTR_flank_hmmer + "/"
                                  + parameter_dic["out_prefix"]
                                  + "_flank_without_paired_LTR_hmmer_log.txt")

        hmmer_commd2 = (hmmscan_execute + " --domtblout "
                        + no_LTR_flank_hmmer_su
                        + " -E " + parameter_dic["HmmerE"]
                        + " --cpu " + parameter_dic["thread"] + " "
                        + ERVprofiles_path + " " + flank_without_ltrs_aa
                        + " > " + no_LTR_flank_hmmer_log)

        print("\n" + ">>> Start HMMER search for ERVLEs without paired-LTRs...")
        hmmer2_process = subprocess.Popen(hmmer_commd2,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT,
                                          universal_newlines=True,
                                          shell=True)

        while True:

            output = hmmer2_process.stdout.readline()
            if hmmer2_process.poll() is not None:
                break

            elif output:
                print(output)


        hmmer2_process.terminate()

        if hmmer_contain(no_LTR_flank_hmmer_su) == "True":
            print("\n" + ">>> Annotate domain...." + "\n")

            # * simplify the result of hmmer *
            no_LTR_flank_hmmer_simplify = (no_LTR_flank_hmmer
                                           + "/" + parameter_dic["out_prefix"]
                                           + "_flank_seq_without_paired_LTR_hmmer_simplify.txt")

            hmmer_head_noLTRs = ("##seq name for HMMER" + "\t"
                                       + "GAG domain"
                                       + "\t" + "DUT domain"
                                       + "\t"
                                       + "AP domain" + "\t"
                                       + "RT domain"
                                       + "\t" + "RNaseH domain"
                                       + "\t"
                                       + "INT domain" + "\t"
                                       + "ENV domain"
                                       + "\n")

            simp_hmmer2 = SimplyHMMER(no_LTR_flank_hmmer_su,
                                      no_LTR_flank_hmmer_simplify,
                                      hmmer_head_noLTRs,parameter_dic)
            simp_hmmer2.run()

            # * domian map to genome *
            annotation_dir = (parameter_dic["out_dir"]
                              + "/Annotatio_and_extract_domain")
            make_dir(annotation_dir)

            noLTR_flank_annotation_dir = (annotation_dir
                                          + "/ERVLEs_without_paired_LTRs")

            make_dir(noLTR_flank_annotation_dir)

            noLTR_flank_map_genome = (noLTR_flank_annotation_dir + "/"
                                      + parameter_dic["out_prefix"]
                                      + "_domain_annotation_in_genome.txt")

            out_map_duplicate_removal = (noLTR_flank_annotation_dir + "/"
                                         + parameter_dic["out_prefix"]
                                         + "_domain_annotation_in_genome_final.txt")

            if LTRharvest_out_size != 0:

                domain_in_genome_map2 = NoPairLTReveMap(pairltr_ervs_map_genome,
                                                        no_LTR_flank_hmmer_simplify,
                                                        noLTR_flank_map_genome,
                                                        out_map_duplicate_removal)
                domain_in_genome_map2.run()

            elif LTRharvest_out_size == 0:
                domain_in_genome_map2 = NoPairLTReveMap2(no_LTR_flank_hmmer_simplify,
                                                        noLTR_flank_map_genome,
                                                        out_map_duplicate_removal)

                domain_in_genome_map2.run()


            #  *  extract the potential ERVs ervs and domain in frank sequence *

            doamin_seq_dir2 = (noLTR_flank_annotation_dir + "/"
                               + parameter_dic["out_prefix"] + "_domain_seqs")

            make_dir(doamin_seq_dir2)
            potential_erv_path = (noLTR_flank_annotation_dir + "/"
                                  + parameter_dic["out_prefix"]
                                  + "_extracted_potential_ERVLEs.fas")

            get_potential_erv = GetPotentialEves(out_map_duplicate_removal,
                                                 parameter_dic["host_path"],
                                                 doamin_seq_dir2,
                                                 potential_erv_path)

            get_potential_erv.run()


            #  * annotate domain on the potential ervs *

            domain_in_potential_erv = (noLTR_flank_annotation_dir
                                    + "/domain_annotation_in_potential_ERVLEs.txt")

            potential_ervs_annotation = PotentialEveAnno(out_map_duplicate_removal,
                                                         domain_in_potential_erv,
                                                         parameter_dic["out_prefix"])
            potential_ervs_annotation.run()


            print("\n" + ">>> The result was saved in: " + annotation_dir + "\n")

        else:

            print("\n" + "*** Note: no reliable domain of retroviridae was found in the "
                         "these potential ERVLEs without paired-LTRs!")



    print("\n" + ">>> Analysis have finshed.")


    duration = datetime.today().now() - time_start

    # end_memory = calculate_memory()
    #
    # used_memory = end_memory - start_memory
    #
    # print(">>> " + f"Occupied {used_memory}MB memory in total" + "\n")

    print(">>> " + "Take " + str(duration) + " seconds in total." + "\n")

    sys.exit()


if __name__ == "__main__":
    starts()

