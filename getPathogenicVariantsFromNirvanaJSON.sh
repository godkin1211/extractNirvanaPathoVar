#!/bin/bash
## Parameter setting
# Input & output
inputJSONfile=$1
outputFile=$2
JSON_location=${inputJSONfile%/*}
JSON_file=${inputJSONfile##*/}
# target gene-related info: genomic regions & selected transcripts
apc_region="chr5:112707498-112846239"
atm_region="chr11:108222804-108369102"
brca1_region="chr17:43044295-43170327"
brca2_region="chr13:32315086-32400268"
cdkn2a_region="chr9:21967752-21995324"
chek2_region="chr22:28687743-28742422"
mlh1_region="chr3:36993226-37050896"
msh2_region="chr2:47403067-47709830"
msh6_region="chr2:47695530-47810101"
palb2_region="chr16:23603160-23641321"
pms2_region="chr7:5970925-6009130"
apc_transcript="NM_000038"
atm_transcript="NM_000051"
brca1_transcript="NM_007294"
brca2_transcript="NM_000059"
cdkn2a_transcript="NM_000077"
chek2_transcript="NM_007194"
mlh1_transcript="NM_000249"
msh2_transcript="NM_000251"
msh6_transcript="NM_000179"
palb2_transcipt="NM_024675"
pms2_transcript="NM_000535"
all_transcripts=($apc_transcript $atm_transcript $brca1_transcript $brca2_transcript $cdkn2a_transcript $chek2_transcript $mlh1_transcript $msh2_transcript $msh6_transcript $palb2_transcipt $pms2_transcript)
# Nirvana JSON parser commands
jasixCMD="docker run --rm -v $PWD:/work -w /work annotation/nirvana:3.14 dotnet /opt/nirvana/Jasix.dll"

# Check if the JSON file exists
[[ -e $inputJSONfile ]] || { echo "ERROR: JSON file not found"; exit 1; }

# Check if the output file exists
[[ -e $outputFile ]] && { echo -e "\033[1;33m* Warning: output file has already existed! Please change the file name of output file!\033[0m"; exit 1;}

# Executing Jasix.dll to extract called variants in targeted genes
#echo -e "\033[1;32m*** Extract variant informations ***\033[0m"
calledVariantsInTargetedGenes=$($jasixCMD -i $JSON_location/$JSON_file -q $apc_region -q $atm_region -q $brca1_region -q $brca2_region -q $cdkn2a_region -q $chek2_region -q $mlh1_region -q $msh2_region -q $msh6_region -q $palb2_region -q $pms2_region)

# Count the number of variants in targeted genes
numVariantsInTargetedGenes=$(echo $calledVariantsInTargetedGenes | jq '.positions | length')

# Check if there are any pathogenic variants in targeted genes
#echo -e "\033[1;32m* Checking $numVariantsInTargetedGenes variants in targeted genes\033[0m"
touch $outputFile && echo "Gene,Transcript,Genotype,HGVSc,HGVSp,ClinVar_Significance" >> $outputFile

for i in $(seq 0 `expr $numVariantsInTargetedGenes - 1`); do
    numVariantsFoundInSamePosition=$(echo $calledVariantsInTargetedGenes | jq ".positions[$i].variants | length")
    thisPosition=$(echo $calledVariantsInTargetedGenes | jq ".positions[$i]")
    filterResult=$(echo $thisPosition | jq ".filters[0]" | sed 's/\"//g')
    [[ "$filterResult" != "PASS" ]] && continue
    genotype=$(echo $thisPosition | jq '.samples[0].genotype' | sed 's/\"//g')
    refAllele=$(echo $thisPosition | jq ".refAllele" | sed 's/\"//g')
    altAllele=$(echo $thisPosition | jq ".altAlleles[0]" | sed 's/\"//g')
    chrom=$(echo $thisPosition | jq ".chromosome" | sed 's/\"//g')
    pos=$(echo $thisPosition | jq ".position")
    transcript="NULL"
    gene="NULL"
    hgvsc="NULL"
    hgvsp="NULL"
    for j in $(seq 0 `expr $numVariantsFoundInSamePosition - 1`); do
        thisVariant=$(echo $thisPosition | jq ".variants[$j]")
        numClinVarRecordsInThisVariant=$(echo $thisVariant | jq ".clinvar | length")
        numTranscriptsAssociatedWithThisVariant=$(echo $thisVariant | jq ".transcripts | length")
        clinvar_significance="NULL"
        for k in $(seq 0 `expr $numClinVarRecordsInThisVariant - 1`); do
            thisClinVarRecord=$(echo $thisVariant | jq ".clinvar[$k]")
            clinvar_rec_refAllele=$(echo $thisClinVarRecord | jq ".refAllele" | sed 's/\"//g')
            clinvar_rec_altAllele=$(echo $thisClinVarRecord | jq ".altAllele" | sed 's/\"//g')
            if [[ "$refAllele" == "$clinvar_rec_refAllele" && "$altAllele" == "$clinvar_rec_altAllele" ]]; then
                clinvar_significance=$(echo $thisClinVarRecord | jq ".significance[]")
            fi
        done
        for k in $(seq 0 `expr $numTranscriptsAssociatedWithThisVariant - 1`); do
            thisTranscriptRecord=$(echo $thisVariant | jq ".transcripts[$k]")
            thisTranscript=$(echo $thisTranscriptRecord | jq ".transcript")
            for trpt in ${all_transcripts[@]}; do
                if [[ "$thisTranscript" =~ "$trpt" ]]; then
                    transcript=$thisTranscript
                    hgvsc=$(echo $thisTranscriptRecord | jq ".hgvsc")
                    hgvsp=$(echo $thisTranscriptRecord | jq ".hgvsp")
                    [[ "$hgvsc" == "null" ]] && hgvsc="NULL"
                    [[ "$hgvsp" == "null" ]] && hgvsp="NULL"
                    case $trpt in 
                        NM_000038) gene="APC"
                            ;;
                        NM_000051) gene="ATM"
                            ;;
                        NM_007294) gene="BRCA1"  
                            ;;
                        NM_000059) gene="BRAC2"
                            ;;
                        NM_000077) gene="CDKN2A"
                            ;;
                        NM_007194) gene="CHEK2"
                            ;;
                        NM_000249) gene="MLH1"
                            ;;
                        NM_000251) gene="MSH2"
                            ;;
                        NM_000179) gene="MSH6"
                            ;;
                        NM_024675) gene="PALB2"
                            ;;
                        NM_000535) gene="PMS2"
                            ;;
                        *) gene="NULL"
                    esac
                    echo -e "$gene,$transcript,$genotype,$hgvsc,$hgvsp,$clinvar_significance" >> $outputFile
                    continue
                fi
            done
            [[ "$transcript" != "NULL" ]] && continue
        done
    done
done

echo -e "\033[1m***** Job Completed! *****\033[0m"
