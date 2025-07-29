##### Collect metadata from omnicrobe database #####

#### Use taxids to get info on known habitats, phenotypes and use for all assembled MAGs

#### VARIABLES AND WORKING DIRECTORY ####
scriptdir=".."
outdir="../../output/community_analysis"

#### PREP INPUT ####
cut -d, -f 2 $outdir/names_to_ids_filt.csv | sort | uniq > $outdir/taxids.tmp

# Remove NAs and header
sed -i '/NA/d' $outdir/taxids.tmp
sed -i '/ids/d' $outdir/taxids.tmp

#### RUN SCRIPT ####
python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "habitat"
python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "phenotype"
python $scriptdir/access_omnicrobe_db.py $outdir/taxids.tmp $outdir "use"

rm $outdir/taxids.tmp
