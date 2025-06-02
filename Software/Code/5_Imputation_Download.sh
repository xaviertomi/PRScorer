#!/bin/bash

# ./run_pipeline.sh "curl -sL https://... | bash" TestPCA1

# Validate inputs
if [ -z "$1" ]; then
    echo "Error: Please provide a download command."
    exit 1
fi

# Assignation des arguments
cmd=$1 #commande curl
DL_COMMAND="${cmd}"
name=$2  

# Constants
PWD="AS4P3RS3C4R3DP5SSW0RD" #Ecrire le mot de passe en hard dans l'unzip
PIPELINE_DIR="/storage/home/x_tomisinecvandoren/Ksein/FinalPipelineV3"
Post_imputation_path="${PIPELINE_DIR}/Data/${name}/VCF/post_imputation"

# Creation du dossier pour téléchargement des résultats

mkdir -p "${Post_imputation_path}"
cd "${Post_imputation_path}" || exit 1

# Exécution de la commande passée
echo "Exécuion de la commande: ${DL_COMMAND}"
eval "${DL_COMMAND}"

# Création du dossier unzip
mkdir -p "$Post_imputation_path/unzip"

# Unzip de tout les fichier .zip avec le mot de passe
for z in "$Post_imputation_path"/*.zip; do
    unzip -P AS4P3RS3C4R3DP5SSW0RD $z -d $Post_imputation_path/unzip
done

# Fusion des fichiers .dose.vcf.gz en en un seul fichier .vcf contenant toutes les chromosomes
bcftools concat -O v -o "${Post_imputation_path}/${name}_chr_all.vcf" "$Post_imputation_path"/unzip/chr*.dose.vcf

# Nettoyage des fichiers si le fichier .vcf finale existe
echo '''
if [ -f "${Post_imputation_path}/${name}_chr_all.vcf" ]; then
    rm -f "$Post_imputation_path"/unzip/*
    rm -f "$Post_imputation_path"/*.zip
fi
'''
# Retour au dossier de travail
cd ${PIPELINE_DIR}/Software/Code || exit 1

# Exécution du code de la PCA
bash "${PIPELINE_DIR}/Software/Code/Ethnic_PCA/1-determine-ancestry-by-PCA.sh" "${name}"

#Execution du code de scoring
#sh "${PIPELINE_DIR}/Software/Code/Scoring/Score_module.sh" "${name}"