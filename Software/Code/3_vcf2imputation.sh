name=$1
token=$2
PIPELINE_DIR="/storage/home/x_tomisinecvandoren/Ksein/FinalPipelineV3"
VCFLocation="${PIPELINE_DIR}/Data/"$name"/VCF/"
echo "Step 3: Starting VCF conversion for Impution..." >> ${PIPELINE_DIR}/Software/Code/Progress.txt

#create per_chr folder if not exist 
per_chr="${PIPELINE_DIR}/Data/"$name"/VCF/per_chr"
mkdir -p $per_chr

#split chr 
for i in $(seq 1 22) X Y; 
do
    echo $i
    bcftools filter -t chr$i "${VCFLocation}${name}_no_monomorph.vcf" -o "${per_chr}/${i}.vcf.gz" -O z

    bcftools sort "${per_chr}/${i}.vcf.gz" -O v > $per_chr/$i.sorted.vcf
    #rm $i.vcf.gz
    sed -i "s/^chr$i/$i/" $per_chr/$i.sorted.vcf
    bgzip $per_chr/$i.sorted.vcf
done

echo "Step 3 done: VCF conversion for Impution done" >> ${PIPELINE_DIR}/Software/Code/Progress.txt
echo "Step 4: Sending the file to the API... - Please check your email in the next hour" >> ${PIPELINE_DIR}/Software/Code/Progress.txt
#----------- 2 options to send the file to the Imputaion server -----------
#send to API
python3 ${PIPELINE_DIR}/Software/Code/4_Imputation_APY.py $token $name

#OR use imputationbot (a permit to download 50+ files at once is necessary)
#in this case, 4_Imputation_APY.py and 5_Imputation_Download is not used, plus the code "I have received an imputation email" can be deleted
: '
# bash imputationbot impute \
    --files ${PIPELINE_DIR}/Data/${name}/VCF/per_chr/*.sorted.vcf.gz \
    --refpanel 1000g-phase-3-v5 \
    --population super_pop \
    --autoDownload \
    --password AS4P3RS3C4R3DP5SSW0RD \
    --name ${name} \
    --build hg19

path="${PIPELINE_DIR}/Data/${name}/VCF/post_imputation"
mkdir -p ${path}
mv "/storage/home/x_tomisinecvandoren/Ksein/FinalPipeline3/Software/Code/imputationV100/job-*-${name}/local/*.zip" "${path}"

# Création du dossier unzip
# Create unzip folder
mkdir -p "$path/unzip"

# Unzip de tout les fichier .zip avec le mot de passe
# Unzipping all .zip files with the password
for z in "$path"/*.zip; do
    unzip -P AS4P3RS3C4R3DP5SSW0RD $z -d $path/unzip
done

# Fusion des fichiers .dose.vcf.gz en en un seul fichier .vcf contenant toutes les chromosomes
# Merging .dose.vcf.gz files into a single .vcf file containing all chromosomes
bcftools concat -O v -o "${path}/${name}_chr_all.vcf" "$path"/unzip/chr*.dose.vcf.gz

# Exécution du code de la PCA
# Execution of the PCA code
bash "${PIPELINE_DIR}/Software/Code/Ethnic_PCA/1-determine-ancestry-by-PCA.sh" "${name}" &

#Execution du code de scoring
# Execution of the scoring code
#sh "${PIPELINE_DIR}/Software/Code/Scoring/Score_module.sh" "${name}"
' #FIN commentaire
