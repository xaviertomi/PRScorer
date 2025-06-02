#!/bin/bash

# ./run_pipeline.sh "curl -sL https://... | bash" TestPCA1

# Validate inputs
# Check if the download command is provided
if [ -z "$1" ]; then
    echo "Error: Please provide a download command."
    exit 1
fi

# Assignation des arguments
# arguments assignment
cmd=$1 #commande curl
DL_COMMAND="${cmd}"
name=$2  

# Constants
PWD="AS4P3RS3C4R3DP5SSW0RD" # Write the password in hard in the unzip cmd (line 40) for more stability # Ecrire le mot de passe en dur dans la commande unzip (ligne 40) pour plus de stabilité
PIPELINE_DIR="/storage/home/x_tomisinecvandoren/Ksein/FinalPipelineV3"
Post_imputation_path="${PIPELINE_DIR}/Data/${name}/VCF/post_imputation"

# Creation du dossier pour téléchargement des résultats
#Create the folder for the downloaded results
mkdir -p "${Post_imputation_path}"
cd "${Post_imputation_path}" || exit 1

# Exécution de la commande passée
# Execute the provided download command
echo "Exécuion de la commande: ${DL_COMMAND}"
eval "${DL_COMMAND}"

# Création du dossier unzip
# Create the unzip folder
mkdir -p "$Post_imputation_path/unzip"

# Unzip de tout les fichier .zip avec le mot de passe
# Unzipping all .zip files with the password
for z in "$Post_imputation_path"/*.zip; do
    unzip -P MYSECUREPWD $z -d $Post_imputation_path/unzip
done

# Fusion des fichiers .dose.vcf.gz en en un seul fichier .vcf contenant toutes les chromosomes
# Merging all .dose.vcf.gz files into a single .vcf file containing all chromosomes
bcftools concat -O v -o "${Post_imputation_path}/${name}_chr_all.vcf" "$Post_imputation_path"/unzip/chr*.dose.vcf

# Nettoyage des fichiers si le fichier .vcf finale existe
# Cleaning up files if the final .vcf file exists, depending on how much data must be kept
echo '''
if [ -f "${Post_imputation_path}/${name}_chr_all.vcf" ]; then
    rm -f "$Post_imputation_path"/unzip/*
    rm -f "$Post_imputation_path"/*.zip
fi
'''
# Retour au dossier de travail
# Return to the working directory
cd ${PIPELINE_DIR}/Software/Code || exit 1

# Exécution du code de la PCA
# Execution of the PCA script
bash "${PIPELINE_DIR}/Software/Code/Ethnic_PCA/1-determine-ancestry-by-PCA.sh" "${name}" &

#Execution du code de scoring
sh "${PIPELINE_DIR}/Software/Code/Scoring/Score_module.sh" "${name}"