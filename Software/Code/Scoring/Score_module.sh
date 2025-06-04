#récup name
name=$1
PIPELINE_DIR="My_PIPELINE_DIR" # Replace with your actual pipeline directory path
vcfLocation="${PIPELINE_DIR}/Data/${name}/VCF/post_imputation/${name}_chr_all.vcf"
PRSModelPath="PaTh_to_my_CanRisk_PRSmodels" # Replace with the actual path to your PRS models
EthnicModel=("BCAC_313_PRS.prs" "BCAC_307_PRS-UKB_southAsian.prs" "BCAC_307_PRS-UKB_eastAsian.prs" "BCAC_307_PRS-UKB_african.prs")

ScorePath="${PIPELINE_DIR}/Data/${name}/PRS"
mkdir -p $ScorePath

mkdir -p "${PIPELINE_DIR}/Results/${sam}"

#récupération des du nom des échantillons
#recovery of sample names
samples=$(bcftools query -l $vcfLocation)

#Filtre les variants PRS pour chaque ethnie
#Filter PRS variants for each ethnicity
for sam in $samples;
    do
    ResultFile="${PIPELINE_DIR}/Results/${sam}/${sam}_Score.txt"
    echo "${sam}" > $ResultFile
    for PRSModel in "${EthnicModel[@]}";do

        #Debugger echos
        #echo "vcf2prs -g "${ScorePath}/${name}_${PRSModel}.vcf" -s "${sam}" "${PRSModelPath}/${EthnicModel}" >> $ResultFile"
        #echo "VOICI LE SCORE PATH: ${ScorePath}/${name}_${PRSModel}.vcf"
    
        echo "Processing sample: ${sam} with model: ${PRSModel}"
        #Ecris le modèle PRS dans le fichier de résultats
        #wrte the PRS model to the result file
        echo -e "----" >> $ResultFile
        echo "${PRSModel}" >> $ResultFile

        #python .py ~prs ~vcf ~output
        python3 ${PIPELINE_DIR}/Software/Code/Scoring/FindPRS_onlyGT.py "${PRSModelPath}/${EthnicModel}" "${vcfLocation}" "${ScorePath}/${name}_${PRSModel}.vcf"

        #vcf2prs -g ~vcf -s ~sample ~prs
        vcf2prs -g "${ScorePath}/${name}_${PRSModel}.vcf" -s "${sam}" "${PRSModelPath}/${EthnicModel}" >> $ResultFile
        
        done
    done

