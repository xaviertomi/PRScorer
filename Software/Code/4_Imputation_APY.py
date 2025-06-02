import requests
import sys

url = 'https://imputationserver.helmholtz-munich.de/api/v2'
#https://imputationserver.helmholtz-munich.de/api/v2/jobs/submit/imputationserver@2.1.5
token = sys.argv[1]
project_name = sys.argv[2]


#=========================================PARAMETERS=========================================
data = {
    'job-name': project_name,
    'refpanel': '1000g-phase-3-v5',
    'build': 'hg19',
    'population': 'super_pop', 
    'password': 'AS4P3RS3C4R3DP5SSW0RD' # Password for the job, can be anything and must be the exact same in 5_Imputation_Download.sh
}

headers = {'X-Auth-Token' : token}
#=========================================UPLOAD=========================================

#code pour upload d'un seul fichier
#code if needed to upload a single VCF file
""" # 1 file
vcf = '/storage/home/x_tomisinecvandoren/Ksein/sandbox/All_hg19/22.non_chr_sorted.vcf'
files = {'files' : open(vcf, 'rb')}
r = requests.post(url + '/jobs/submit/imputationserver@2.1.5',
                  files=files, data=data, headers=headers) """


#code pour upload de plusieurs fichiers
#code to upload multiple VCF files
per_chr_positon="/storage/home/x_tomisinecvandoren/Ksein/FinalPipelineV3/Data/"+project_name+"/VCF/per_chr/"

# Définition des fichiers en dur (chr1 à 22)
#Hardcoded files from chr1 to chr22
vcf_files = [
    per_chr_positon+"/1.sorted.vcf.gz",
    per_chr_positon+"/2.sorted.vcf.gz",
    per_chr_positon+"/3.sorted.vcf.gz",
    per_chr_positon+"/4.sorted.vcf.gz",
    per_chr_positon+"/5.sorted.vcf.gz",
    per_chr_positon+"/6.sorted.vcf.gz",
    per_chr_positon+"/7.sorted.vcf.gz",
    per_chr_positon+"/8.sorted.vcf.gz",
    per_chr_positon+"/9.sorted.vcf.gz",
    per_chr_positon+"/10.sorted.vcf.gz",
    per_chr_positon+"/11.sorted.vcf.gz",
    per_chr_positon+"/12.sorted.vcf.gz",
    per_chr_positon+"/13.sorted.vcf.gz",
    per_chr_positon+"/14.sorted.vcf.gz",
    per_chr_positon+"/15.sorted.vcf.gz",
    per_chr_positon+"/16.sorted.vcf.gz",
    per_chr_positon+"/17.sorted.vcf.gz",
    per_chr_positon+"/18.sorted.vcf.gz",
    per_chr_positon+"/19.sorted.vcf.gz",
    per_chr_positon+"/20.sorted.vcf.gz",
    per_chr_positon+"/21.sorted.vcf.gz",
    per_chr_positon+"/22.sorted.vcf.gz"
]

# Création de la liste de tuples pour l'envois des fichiers
# Create a list of tuples for sending files
files = [('files', open(vcf, 'rb')) for vcf in vcf_files]
r = requests.post(url + '/jobs/submit/imputationserver@2.1.5',
                  files=files, data=data, headers=headers)

##=========================================POST=========================================
if r.status_code != 200:
  print(r.json()['message'])
  raise Exception('POST /jobs/submit/imputationserver@2.1.5 {}'.format(r.status_code))

# print response and job ID
print(r.json()['message'])
print(r.json()['id'])

#=========================================DOWNLOAD=========================================
#Pas moyen de DL -> voir imputationbot ou uploader manuellement le lien
#No way to download -> see imputationbot or manually upload the link

print("Veuillez vérifier vos mails dans les 30~40 min qui suivent")
