params.filePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/ALL/08.info_recalculate/'
params.samplelist = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_project/sample_ls/omics_lst'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Genome/OMICS'

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
    vcf = filePath + chr + '.tommo.snp.reinfo.vcf.gz'
    pass_vcf = chr + '.pass.mac1.vcf.gz'
    pass_vcf_tbi = chr + '.pass.mac1.vcf.gz.tbi'

    """
    bcftools view --threads 2 -S $samplelist --force-samples -Ou $vcf | bcftools view --threads 2 -f"PASS" -c1 -Oz -o $pass_vcf
    bcftools index --threads 2 -t $pass_vcf
    """
}

process info_recalculate {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/02.info_recalculate", mode: 'symlink'

    input:
    tuple val(chr), file(snp_vcf), file(snp_vcf_tbi) from select_pass_ch

    output:
    tuple val(chr), file(info_vcf), file(info_vcf_tbi) into bed_prepare_ch

    script:
    info_vcf = chr + '.omics.tommo.reinfo.vcf.gz'
    info_vcf_tbi = chr + '.omics.tommo.reinfo.vcf.gz.tbi'
    """
    bcftools +fill-tags ${snp_vcf} -Oz -o ${info_vcf} --threads 2 -- -t 'AF,AC,AN,DP:1=int(sum(FORMAT/DP))'
    bcftools index -t ${info_vcf} --threads 2
    """
}