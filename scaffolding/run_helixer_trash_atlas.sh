
### ON ATLAS




## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats

six=crefra
FASTA=Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa        
mkdir ${six}_TRASH2
     
six=cserru
FASTA=Cs-KelloggPI219580-DRAFT-PanAnd-1.0.chrSuperScaf.fa  
mkdir ${six}_TRASH2     
six=hconto
FASTA=Hc-AUB53_1-DRAFT-PanAnd-1.0.chrSuperScaf.fa   
mkdir ${six}_TRASH2          
six=ppanic
FASTA=Pi-Clark-DRAFT-PanAnd-1.0.chrSuperScaf.fa    
mkdir ${six}_TRASH2        





conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=crefra; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa  -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}_TRASH2 -p 46"
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=cserru; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/Cs-KelloggPI219580-DRAFT-PanAnd-1.0.chrSuperScaf.fa  -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}_TRASH2 -p 46"
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=hconto; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/Hc-AUB53_1-DRAFT-PanAnd-1.0.chrSuperScaf.fa  -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}_TRASH2 -p 46"
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=ppanic; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/Pi-Clark-DRAFT-PanAnd-1.0.chrSuperScaf.fa  -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}_TRASH2 -p 46"

### then helixer
cd ..
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=crefra; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path Cr-AUB069-DRAFT-PanAnd-1.0.chrSuperScaf.fa   --lineage land_plant --gff-output-path ${six}.chrSuperScaf.helixer.gff3'
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=cserru; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path Cs-KelloggPI219580-DRAFT-PanAnd-1.0.chrSuperScaf.fa  --lineage land_plant --gff-output-path ${six}.chrSuperScaf.helixer.gff3'
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=hconto; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path Hc-AUB53_1-DRAFT-PanAnd-1.0.chrSuperScaf.fa  --lineage land_plant --gff-output-path ${six}.chrSuperScaf.helixer.gff3'
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=ppanic; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path Pi-Clark-DRAFT-PanAnd-1.0.chrSuperScaf.fa  --lineage land_plant --gff-output-path ${six}.chrSuperScaf.helixer.gff3'

## not doing this yet - kind of hard?
## need to set up subphaser input on cbsu first!!!!!!



conda activate SubPhaser
cd subphaser
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}BothHapsUnfilteredaggressive_subphaserinput.txt .

## generate subphaser in put through my script from anchorwave output (need to improve usability)
six=hcontoCHR
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=hcontoCHR; subphaser -i ../Hc-AUB53_1-DRAFT-PanAnd-1.0.chrSuperScaf.fa -c ${six}_subphaserinput.txt -pre ${six}_ -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"

six=tdacn1
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=tdacn1; subphaser -i ../Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fa -c ${six}_subphaserinput.txt -pre ${six}_ -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"

six=tdacs1
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=tdacs1; subphaser -i ../Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fa -c ${six}_subphaserinput.txt -pre ${six}_ -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"



six=avirgi
FASTA=Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta  
mkdir ${six}_TRASH2        

conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=avirgi; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta   -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}_TRASH2 -p 46"
cd ..
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=avirgi; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta   --lineage land_plant --gff-output-path ${six}.helixer.gff3'







