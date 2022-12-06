# nf-kallisto

Create the test directory:
```
mkdir -p ~/nf-kallisto-test/raw_data
```

Download the demo data:
```
cd ~/nf-kallisto-test/raw_data
curl -J -O https://datashare.mpcdf.mpg.de/s/jcEaS5vqpJO0lOy/download
curl -J -O https://datashare.mpcdf.mpg.de/s/XHanbnjfvQ9rACD/download
curl -J -O https://datashare.mpcdf.mpg.de/s/sIebkRdMfMSweq2/download
curl -J -O https://datashare.mpcdf.mpg.de/s/zoNxS9vRI7jl77y/download
curl -J -O https://datashare.mpcdf.mpg.de/s/0WHGNIhjJC792lY/download
curl -J -O https://datashare.mpcdf.mpg.de/s/ZlM0lWKPh8KrP6B/download
curl -J -O https://datashare.mpcdf.mpg.de/s/o3O6BKaEXqB7TTo/download
```

Download the paramaters files:
```
cd ~/nf-kallisto-test
curl -J -O https://raw.githubusercontent.com/mpg-age-bioinformatics/nf-kallisto/main/params.json
```

Run the workflow:
```
RELEASE=1.0.0
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry images -profile local && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry get_genome -profile local && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry write_cdna -profile local && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry index -profile local && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry check_strand -profile local && \
nextflow run mpg-age-bioinformatics/nf-kallisto -r ${RELEASE} -params-file params.json -entry map_reads -profile local
```


## Contributing

Check the last tag, add a new one, push it and make a release|:
```
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push
git push origin --tags
gh release create <tag> 
```