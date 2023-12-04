
## prune bs < 10
conda activate anchorwave_new ## conda install -c bioconda newick_utils
nw_ed  paspalum_anchors_aster.tre 'i & b<=10' o > paspalum_anchors_aster-BS10.tre

## run
../../cbsu_projects/hackathon_december2022/ASTER-Linux/bin/astral-pro -t 24 -o paspalum_anchors_ASTRALout.tre paspalum_anchors_aster.tre
../../cbsu_projects/hackathon_december2022/ASTER-Linux/bin/astral-pro -t 24 -o paspalum_anchors_ASTRALout-BS10.tre paspalum_anchors_aster-BS10.tre
