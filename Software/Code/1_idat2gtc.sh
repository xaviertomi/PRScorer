name=$1
PIPELINE_DIR="My_PIPELINE_DIR" # Replace with your actual pipeline directory path
idatLocation="${PIPELINE_DIR}/Data/"$name"/IDAT"
BeadChipType=$2
BeadChipFile="${PIPELINE_DIR}/Software/Cluster_Manifest/"$BeadChipTsype
bpm=$BeadChipFile"/Hg19/*.bpm"
csv=$BeadChipFile"/Hg19/*.csv"
egt=$BeadChipFile"/Cluster/*.egt"


echo "" > ${PIPELINE_DIR}/Software/Code/Progress.txt
echo "Step 1: Starting IDAT to GTC conversion..." >> ${PIPELINE_DIR}/Software/Code/Progress.txt
GTCLocation="${PIPELINE_DIR}/Data/"$name"/GTC"
mkdir -p $GTCLocation

#idat->gtc
bcftools +idat2gtc -i $idatLocation \
    -b $bpm \
    -e $egt \
    -o $GTCLocation
    #--snp-map snipmap
echo "Step1 done: idat2gtc done" >> ${PIPELINE_DIR}/Software/Code/Progress.txt

#GTC->VCF
sh ${PIPELINE_DIR}/Software/Code/2_gtc2vcf.sh   $name   $BeadChipType
