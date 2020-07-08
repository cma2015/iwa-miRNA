
set -x

script_path=$1/scripts
PATH="$PATH:${script_path}"

outxt=$4

out_dir=${outxt%.*}

if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}
fi

cd ${out_dir}

python ${script_path}/databases_mining.py --curpath $1 --species "$2" \
--databases $3 --outpath ./

Rscript -e "rmarkdown::render('Overview.Rmd')"

mv Overview.html $outxt
mkdir -p ${outxt%.*}_files
mv Overview_files ${outxt%.*}_files/Overview_files

zip Collection.zip miRBase.txt PlantsmallRNAgenes.txt PmiREN.txt sRNAanno.txt

mv Collection.zip $5
