muscle -in data/all.consensus.wref2.fasta -out data/all.consensus.wref2.aln.fasta
python ./scripts/mask_alignment_using_vcf.py -i data/all.consensus.wref2.aln.fasta -o data/all.consensus.wref2.aln.mask.fasta -v data/problematic_sites_sarsCov2.vcf -r Wuhan-Hu-1
mkdir data/iqtree_output
iqtree -redo -nt AUTO -me 0.05 -bb 1000 -wbtl -czb -m GTR --prefix data/iqtree_output/michigan_timetree -s data/all.consensus.wref2.aln.mask.fasta -o Wuhan/WH01/2019
treetime --tree data/iqtree_output/michigan_timetree.treefile --dates data/all_sample_dates.tsv --aln data/all.consensus.wref2.aln.mask.fasta --outdir data/treetime_output --confidence --clock-filter 0 --clock-rate 0.0008 --clock-std-dev 0.0004 --coalescent skyline --keep-root
