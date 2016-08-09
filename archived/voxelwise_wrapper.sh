subjDataName="/import/monstrum2/Users/angelgar/nbackDprime/n404_longitudinal_dataset_nback_dprime.rds"
OutDirRoot="/import/monstrum2/Users/angelgar/nbackDprime/voxelWiseAnalysis"
namePaths="nback.path"
maskName="/import/speedy/eons/group_results_n1601/frac2back/masks/mask_n1129_FINAL_binary.nii.gz"
smooth=0
inclusionName="dimentional.include"
subjID="bblid"
covsFormula="~sex+all.dprime+s(age)"
randomFormula="~(1|bblid)"
ncores=20

Rscript ~/nbackDprime/gamm4_voxelwise.R -c $subjDataName -o $OutDirRoot -p $namePaths -m $maskName -i $inclusionName -u $subjID -f $covsFormula -r $randomFormula -n 20 -s 0



logfile="/home/agarza/logs"
errfile="/home/agarza/logs"


qsub -V -S /share/apps/R/R-3.1.1/bin/Rscript -cwd -o ${logfile} -e ${errfile} -binding linear:10 -pe unihost 10 /home/agarza/gamm4_voxelwise.R -c $cov -i $in -o $out -s $outs -a $imageName -d $imageIds -m $mask -f $covariates -n 10





subjDataName <- "/import/monstrum2/Users/angelgar/nbackDprime/n404_longitudinal_dataset_nback_dprime.rds"
OutDirRoot <- "/import/monstrum2/Users/angelgar/nbackDprime/voxelWiseAnalysis"
namePaths <- opt$imagepaths
maskName  <- opt$mask
smooth <- opt$smoothing
inclusionName <- opt$inclusion
covsFormula <- opt$formula
randomFormula <- opt$random
ncores <- opt$numbercores