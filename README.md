## POI: Personal Omics Interpreter

![PyPI](https://img.shields.io/pypi/v/panno?color=pink)  ![Conda](https://img.shields.io/conda/v/lyaqing/panno?color=blue&label=conda) ![AppVeyor](https://img.shields.io/appveyor/build/PreMedKB/PAnno)

POI (Personal Omics Interpreter) is an efficient and user-friendly tool that assists clinicians and researchers in resolving patients' multi-omics therapeutic biomarkers to obtain treatments supported by clinical evidence or potentially feasible.

```

### Usage
Once installed, you can use PAnno by navigating to your VCF file and entering the corresponding three-letter abbreviation of the population:

```Shell
panno -s sample_id -i germline_vcf -p population -o outdir
```

* Required arguments
```Shell
-s, --sample_id TEXT            Sample ID that will be displayed in the PAnno report.

-i, --germline_vcf TEXT         Unannotated VCF file, preferably germline variant.

-p, --population [AAC|AME|EAS|EUR|LAT|NEA|OCE|SAS|SSA]
                                The three-letter abbreviation for biogeographic groups:
                                AAC (African American/Afro-Caribbean), AME (American),
                                EAS (East Asian), EUR (European), LAT (Latino),
                                NEA (Near Eastern), OCE (Oceanian),
                                SAS (Central/South Asian), SSA (Sub-Saharan African).

-o, --outdir TEXT               Create report in the specified output path.
```

### Input data
#### 1. Germline VCF file

PAnno directly uses the NGS-derived germline VCF file as input and assumes it has undergone quality control. Therefore, if the VCF file is of poor quality, inaccurate diplotypes and inappropriate clinical recommendations may be reported.

PAnno requires the VCF file aligned to the GRCh38 reference genome given the increasing generality and the built-in diplotype definition dependency version.


#### 2. Population
There are nine biogeographic groups supported by PAnno. Please use the ***three-letter abbreviation*** as input. This is to prevent errors caused by special symbols such as spaces.

**AAC** (African American/Afro-Caribbean), **AME** (American), **EAS** (East Asian), **EUR** (European), **LAT** (Latino), **NEA** (Near Eastern), **OCE** (Oceanian), **SAS** (Central/South Asian), **SSA** (Sub-Saharan African).

More information is available at https://www.pharmgkb.org/page/biogeographicalGroups.

### Output data

The report is created in `${sample_id}.html` at the `outdir` by default.

For more detailed instructions, run `panno -h`.
