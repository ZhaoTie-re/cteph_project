params.filePath = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/VQSR.v4/'
params.samplelist = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/sample_ls/cteph_naga_lst'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/ALL'

params.TommoFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/database/Tommo_hg38/tommo-54kjpn-20230626r3-GRCh38-af-autosome.vcf.gz'
params.TommoFile_tbi = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/database/Tommo_hg38/tommo-54kjpn-20230626r3-GRCh38-af-autosome.vcf.gz.tbi'
params.Tommohdr = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/database/Tommo_hg38/hdr.txt'

Channel
    .of(1..22)
    .set{ chr_ch }

process selectPASS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/01.pass", mode: 'symlink'

    input:
    val filePath from params.filePath
    val chr from chr_ch
    path samplelist from params.samplelist

    output:
    tuple chr, file(pass_vcf), file(pass_vcf_tbi) into select_pass_ch

    script:
    vcf = filePath + 'all.VQSR3.chr' + chr + '.vcf.gz'
    pass_vcf = chr + '.pass.mac1.vcf.gz'
    pass_vcf_tbi = chr + '.pass.mac1.vcf.gz.tbi'

    """
    bcftools view --threads 2 -S $samplelist --force-samples -Ou $vcf | bcftools view --threads 2 -f"PASS" -c1 -Oz -o $pass_vcf
    bcftools index --threads 2 -t $pass_vcf
    """
}

process Variantfilter {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/02.VTfilter", mode: 'symlink'

    input:
    tuple chr, file(pass_vcf), file(pass_vcf_tbi) from select_pass_ch

    output:
    tuple chr, file(vcf), file(vcf_tbi) into GTfilter_ch

    script:
    vcf = chr + '.vft.vcf.gz'
    vcf_tbi = chr + '.vft.vcf.gz.tbi'

    """
    bcftools view --threads 2 -i 'VQSLOD > 10 & MQ > 58.75' -Oz -o $vcf $pass_vcf
    bcftools index --threads 2 -t $vcf
    """
}

process GTfilter {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/03.GTfilter1", mode: 'symlink'

    input:
    tuple chr, file(pass_vcf), file(pass_vcf_tbi) from GTfilter_ch

    output:
    tuple chr, file(vcf), file(vcf_tbi) into GT_addAF_ch

    script:
    vcf = chr + '.vft.gft1.vcf.gz'
    vcf_tbi = chr + '.vft.gft1.vcf.gz.tbi'

    """
    singularity exec -B /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project:/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/gatk_4.1.3.0.sif gatk --java-options "-Xmx1G" VariantFiltration \
    -R /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/database/refseq/hs38DH.fa \
    -V $pass_vcf \
    --genotype-filter-name "LowGQ" \
    --genotype-filter-expression "GQ<20" \
    --genotype-filter-name "LowDP" \
    --genotype-filter-expression "DP<10" \
    --set-filtered-genotype-to-no-call true \
    --create-output-variant-index true \
    -O $vcf

    echo $vcf
    """
}

process normVCF {
    executor 'slurm'
    queue 'gr10478b'
    time '7d'
    tag "${chr}"

    publishDir "${params.outdir}/04.normVCF_keepAD", mode: 'symlink'

    input:
    tuple chr, file(pass_vcf), file(pass_vcf_tbi) from GT_addAF_ch

    output:
    tuple chr, file(norm_fix_vcf), file(norm_fix_vcf_tbi) into GT2Nocall_ch

    script:
    norm_vcf = chr + '.norm.vcf.gz'
    norm_fix_vcf = chr + '.norm_fix.vcf.gz'
    norm_fix_vcf_tbi = chr + '.norm_fix.vcf.gz.tbi'
    """
    bcftools norm --threads 2 -m- -w 10000 --keep-sum AD -f /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/database/refseq/hs38DH.fa -Oz -o $norm_vcf $pass_vcf
    bcftools view $norm_vcf | sed 's/nan/NaN/g' | bgzip > $norm_fix_vcf
    bcftools index --threads 2 -t $norm_fix_vcf
    rm $norm_vcf
    """
    //bcftools norm --threads 5 -m - -w 10000 -f ~/projects/database/refseq/hs38DH.fa -Oz -o $norm_vcf $pass_vcf
}

process GT2Nocall {
    executor 'slurm'
    queue 'gr10478b'
    time '7d'
    tag "${chr}"

    publishDir "${params.outdir}/05.GT2NCall_PASS_noABB", mode: 'symlink'

    input:
    tuple chr, file(norm_fix_vcf), file(norm_fix_vcf_tbi) from GT2Nocall_ch

    output:
    tuple chr, file(GT2Nocall_vcf), file(GT2Nocall_vcf_tbi) into filterAC_ch
    script:
    GT2Nocall_vcf = chr + '.GT_ABB.vcf.gz'
    GT2Nocall_vcf_tbi = chr + '.GT_ABB.vcf.gz.tbi'

    """
    singularity exec -B /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project:/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/gatk_4.1.3.0.sif gatk --java-options "-Xmx1G" VariantFiltration \
    -V $norm_fix_vcf \
    --genotype-filter-name "NAN" \
    --genotype-filter-expression "AF == 'NaN'" \
    --set-filtered-genotype-to-no-call true \
    --create-output-variant-index true \
    -O $GT2Nocall_vcf
    """
}

process filterAC {
    executor 'slurm'
    queue 'gr10478b'
    time '7d'
    tag "${chr}"

    publishDir "${params.outdir}/06.filtAC_vcf", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from filterAC_ch

    output:
    tuple chr, file(outvcf), file(outvcf_tbi) into addAF_ch

    script:
    outvcf = chr + '.vep.filtAC1.vcf.gz'
    outvcf_tbi = chr + '.vep.filtAC1.vcf.gz.tbi'

    """
    bcftools view --threads 2 -Oz -o $outvcf -c 1 $vcf
    bcftools index --threads 2 -t $outvcf
    """
}

process addTommo_AF {
    executor 'slurm'
    queue 'gr10478b'
    time '7d'
    tag "${chr}"

    publishDir "${params.outdir}/07.addTommo_AF_vcf", mode: 'symlink'
    
    input:
    tuple chr, file(vcf), file(vcf_tbi) from addAF_ch
    path Tommo from params.TommoFile
    path Tommo_tbi from params.TommoFile_tbi
    path header from params.Tommohdr

    output:
    tuple chr, file("${chr}.tommo.vcf.gz"), file("${chr}.tommo.vcf.gz.tbi") into vcf_vep_ch

    script:
    """
    bcftools annotate --threads 2 -a $Tommo -c "TOMMO_AF:=AF" -h $header -Ou $vcf | bcftools annotate --threads 2 -I +'%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o ${chr}.tommo.vcf.gz
    bcftools index --threads 2 -t ${chr}.tommo.vcf.gz
    """
}