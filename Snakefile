ALLFAMILIES = ['helium', 'neon', 'argon', 'krypton']
ALLSAMPLES = ['mother', 'father', 'proband']

rule all:
    input:
        expand('{family}-{sample}-reads{cor}.fq.gz',
               family=ALLFAMILIES, sample=ALLSAMPLES, cor=['', '-cor'])

rule simulategenome:
    input:
        'human.order6.mm',
    output:
        '{family}-refr.fa.gz'
    params:
        numseqs=lambda w: config[w.family]['numchr'],
        chrlen=lambda w: config[w.family]['chrlen'],
        seed=lambda w: config[w.family]['seeds']['genome']
    shell:
        'nuclmm simulate '
        '--order 6 --numseqs {params.numseqs} --seqlen {params.chrlen} '
        '--seed {params.seed} --out - {input} '
        '| gzip -n -c > {output}'

rule simulatetriovariants:
    input:
        '{family}-refr.fa.gz'
    output:
        '{family}.vcf',
        '{family}-mother.fasta.gz',
        '{family}-father.fasta.gz',
        '{family}-proband.fasta.gz',
    params:
        ninherited=lambda w: config[w.family]['ninherited'],
        ndenovo=lambda w: config[w.family]['ndenovo'],
        seed=lambda w: config[w.family]['seeds']['trio']
    shell:
        'kevlar gentrio '
        '--inherited {params.ninherited} --de-novo {params.ndenovo} '
        '--seed {params.seed} --prefix {wildcards.family} '
        '--vcf {output[0]} {input} '
        '&& gzip -n -f {wildcards.family}-{{mother,father,proband}}.fasta'

rule sequencing:
    input:
        '{family}-{sample}.fasta.gz'
    output:
        '{family}-{sample}-reads-1.fq',
        '{family}-{sample}-reads-2.fq',
    params:
        seed=lambda w: config[w.family]['seeds']['sequencing'][w.sample],
        numreads=lambda w: config[w.family]['chrlen'] \
                           * config[w.family]['numchr'] \
                           * 30 / 200  # 30x coverage with 2x100bp PE reads
    shell:
        'wgsim '
        '-e 0.005 -r 0.0 -d 450 -s 50 -N {params.numreads} '
        '-1 100 -2 100 -S {params.seed} '
        '{wildcards.family}-{wildcards.sample}.fasta.gz '
        '{wildcards.family}-{wildcards.sample}-reads-{{1,2}}.fq'

rule trustedkmers:
    input:
        expand('{{family}}-{sample}-reads-{end}.fq', sample=ALLSAMPLES,
               end=[1, 2])
    output:
        '{family}-trustedkmers-bloomfilter.lighter',
    params:
        genomesize=lambda w: config[w.family]['chrlen'] \
                             * config[w.family]['numchr']
    threads: 16
    shell:
        'lighter '
        '-r {input[0]} -r {input[1]} -r {input[2]} -r {input[3]} '
        '-r {input[4]} -r {input[5]} '
        '-K 25 {params.genomesize} -t {threads} '
        '-saveTrustedKmers {output}'

rule errorcorrect:
    input:
        '{family}-{sample}-reads-1.fq',
        '{family}-{sample}-reads-2.fq',
        '{family}-trustedkmers-bloomfilter.lighter',
    output:
        '{family}-{sample}-reads-cor-1.fq',
        '{family}-{sample}-reads-cor-2.fq',
    params:
        genomesize=lambda w: config[w.family]['chrlen'] \
                             * config[w.family]['numchr']
    shell:
        'lighter '
        '-r {input[0]} -r {input[1]} '
        '-K 25 {params.genomesize} '
        '-loadTrustedKmers {input[2]} '
        '&& mv {wildcards.family}-{wildcards.sample}-reads-{{1.cor,cor-1}}.fq '
        '&& mv {wildcards.family}-{wildcards.sample}-reads-{{2.cor,cor-2}}.fq'

rule interleave_compress:
    input:
        '{prefix}-1.fq',
        '{prefix}-2.fq',
    output:
        '{prefix}.fq.gz'
    shell:
        'paste '
        '<(paste - - - - < {input[0]}) '
        '<(paste - - - - < {input[1]}) '
        '| tr "\\t" "\\n" '
        '| gzip -n -c '
        '> {output}'
