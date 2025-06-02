name=$1
PIPELINE_DIR="MY_PIPELINE_DIR" # Replace with your actual pipeline directory path
GTCLocation="${PIPELINE_DIR}/Data/"$name"/GTC"
BeadChipType=$2
BeadChipFile="${PIPELINE_DIR}/Software/Cluster_Manifest/"$BeadChipType
 
bpm=$BeadChipFile"/Hg19/*.bpm"
csv=$BeadChipFile"/Hg19/*.csv"
egt=$BeadChipFile"/Cluster/*.egt"
#token pour API imputation Munich
token="MY_TOKEN" # Replace with your actual token (https://imputationserver.helmholtz-munich.de/index.html#!pages/profile)

echo "Step 2: Starting GTC to VCF conversion..." >> ${PIPELINE_DIR}/Software/Code/Progress.txt

#ref hg19
ref="/storage/data/software/highlander/reference/hg19/ucsc.hg19.fasta"
VCFLocation="${PIPELINE_DIR}/Data/"$name"/VCF"
mkdir -p $VCFLocation

#GTC->VCF
bcftools +gtc2vcf $GTCLocation/*.gtc \
    --bpm $bpm\
    --egt $egt\
    -o $VCFLocation/$name.vcf\
    -O v \
    --fasta-ref $ref \
    --no-version \
    --do-not-check-bpm \


#suppression des monomorphes (0/0 ou 1/1)
#delete monomorphic variants (0/0 or 1/1)
grep '^#\|0/1\|1/0' $VCFLocation/$name.vcf > $VCFLocation/$name"_no_monomorph.vcf" 

echo "Step 2 done: gtc2vcf done" >> ${PIPELINE_DIR}/Software/Code/Progress.txt
sh ${PIPELINE_DIR}/Software/Code/3_vcf2imputation.sh  $name   $token
