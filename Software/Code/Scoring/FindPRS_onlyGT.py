# -*- coding: utf-8 -*-
import csv, sys, os

PRSMODELPATH = "PATH_to_my_PRSmodels_CanRisk" # Replace with the actual path to your PRS models

# Dictionary de recherche basé sur la combinaison (chrom, pos, ref, alt)
#Dictionnary to search by combination (chrom, pos, ref, alt)
search_by_variant = {}

def load_search_dict(search_file):
    """Charge le fichier de recherche PRS dans un dictionnaire pour accès rapide"""
    with open(search_file, "r") as prsfile:
        # Ignorer les 6 lignes de métadonnées
        # Skip the first 6 metadata lines
        for _ in range(6):
            next(prsfile)
        reader = csv.DictReader(prsfile, delimiter=",")
        for row in reader:
            chrom = row["Chromosome"]
            pos = row["Position"]
            ref = row["Reference_Allele"]
            alt = row["Effect_Allele"]
            search_by_variant[(chrom, pos, ref, alt)] = row

def process_vcf(vcf_file, output_file):
    """Traite le fichier VCF : filtre les variants présents dans le dictionnaire
    et ne conserve que le champ GT dans les colonnes FORMAT et échantillons"""
    with open(vcf_file, "r") as vcffile, open(output_file, "w") as output:
        for line in vcffile:
            if line.startswith("#"):
                output.write(line)
                continue

            cols = line.strip().split("\t")
            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]

            # Vérifie si la variante est dans le dictionnaire
            # Check if the variant is in the dictionary
            if (chrom, pos, ref, alt) in search_by_variant:
                # Réécriture de la ligne en ne conservant que le champ GT
                # Rewrite the line keeping only the GT field
                format_fields = cols[8].split(":")
                try:
                    gt_index = format_fields.index("GT")
                except ValueError:
                    continue  # sauter si GT est absent # skip if GT is absent

                # Remplacer le champ FORMAT par uniquement "GT"
                # Replace the FORMAT field with only "GT"
                cols[8] = "GT"

                # Réduction des colonnes d'échantillons aux seules valeurs de GT
                # Reduce sample columns to only GT values
                for i in range(9, len(cols)):
                    sample_fields = cols[i].split(":")
                    cols[i] = sample_fields[gt_index]

                output.write("chr" + "\t".join(cols) + "\n")

if __name__ == "__main__":
    
    PRS_search_file = sys.argv[1]
    vcf_file    = sys.argv[2]
    output_file = sys.argv[3]
    
    '''
    #Debugger
    print("affichage de l'output", output_file)
    print("affichage du vcf", vcf_file)
    print("affichage du search", PRS_search_file)
    '''
    load_search_dict(PRS_search_file)
    process_vcf(vcf_file, output_file)
