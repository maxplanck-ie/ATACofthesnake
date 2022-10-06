rule bamIx:
  input: os.path.join(config['dirs']['bamdir'], "{sample}.bam")
  output: os.path.join(config['dirs']['bamdir'], "{sample}.bam.bai")
  conda: config['envs']['seqtools']
  threads: 5
  shell:'''
  samtools index -@ {threads} {input}
  '''
