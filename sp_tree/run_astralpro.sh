
## prune bs < 10
conda activate anchorwave_new ## conda install -c bioconda newick_utils
nw_ed paspalum_anchors_aster.2023-12-05.tre 'i & b<=10' o > paspalum_anchors_aster-BS10.2023-12-05.tre

## run
../../cbsu_projects/hackathon_december2022/ASTER-Linux/bin/astral-pro -t 24 -o paspalum_anchors_ASTRALout.2023-12-05.tre paspalum_anchors_aster.2023-12-05.tre
../../cbsu_projects/hackathon_december2022/ASTER-Linux/bin/astral-pro -t 24 -o paspalum_anchors_ASTRALout-BS10.2023-12-05.tre paspalum_anchors_aster-BS10.2023-12-05.tre
