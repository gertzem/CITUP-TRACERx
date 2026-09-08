set -euf

mkdir -p results/summary
python scripts/read_mach2_obj.py $(find data/mach2_metient -name '*.refined.tree'  | 
    grep 'T-00*\.refined.tree' | sort) > results/summary/mach2_migrations.tsv
python scripts/read_metient_obj.py $(find data/mach2_metient -name '*p.lung.pkl.gz'  | sort) > results/summary/metient_migrations.tsv