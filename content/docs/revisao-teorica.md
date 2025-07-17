---
title: Revisão teórica
date: 2024-02-17
weight: 1
---

## Controle de Qualidade Pré-mapeamento

{{% steps %}}

### Arquivos FASTQ

Um arquivo FASTQ é o formato padrão para armazenar sequências de leitura brutas de sequenciamento de nova geração (NGS), como as geradas por plataformas Illumina. Ele contém não apenas a sequência de nucleotídeos lida, mas também uma pontuação de qualidade para cada base, permitindo a avaliação da confiabilidade dos dados.

![FASTQ](estrutura_fastq.png "Estrutura de um arquivo FASTQ.")

### Estrutura de um arquivo FASTQ

Cada leitura em um arquivo FASTQ é composta por quatro linhas:

**1. Linha de Identificação (Header):**  
Começa com `@` seguido por um identificador único da leitura (ex: ID do instrumento, número da corrida, coordenadas no *flowcell*).  
Exemplo:  
`@M02568:1:000000000-A227B:1:1101:17540:2009 1:N:0:1`

**2. Linha da Sequência:**  
Contém a sequência de nucleotídeos (A, T, C, G e N para bases incertas).  
Exemplo:  
`GTTTCCCCTCGCTGCGAAGGGCCGTCCGAGCGAAGCCGCAGGGTTTCAGCCATGGC`

**3. Linha Separadora:**  
Começa com `+` e pode opcionalmente repetir o identificador da leitura. Geralmente, é usada apenas para separar a sequência da linha de qualidade.  
Exemplo:  
`+`

**4. Linha de Qualidade (Phred Score):**  
Contém uma string de caracteres ASCII que representam a pontuação de qualidade Phred para cada base na linha de sequência correspondente. Uma pontuação mais alta indica uma maior confiança na base.  
Exemplo:  
`AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ`


[Artigo que sistematizou a estrutura do arquivo FASTQ.](https://academic.oup.com/nar/article/38/6/1767/3112533)

### FastQC

O **FastQC** é uma ferramenta essencial para a análise inicial da qualidade de arquivos FASTQ. Ele gera um relatório detalhado que sumariza diversas métricas de qualidade das leituras brutas, permitindo identificar potenciais problemas antes do mapeamento.

**Como funciona:**

O FastQC analisa um ou mais arquivos FASTQ e produz um relatório HTML com gráficos e tabelas para cada uma das seguintes métricas:

![FastQC](fastqc.png "Exemplo de relatório de qualidade gerado no FastQC a partir de um arquivo FASTQ.")

- **Qualidade por Base (Per base sequence quality):**  
  Mostra a distribuição das pontuações de qualidade Phred em cada posição das leituras. Quedas na qualidade em determinadas posições (especialmente no final das leituras) são comuns e importantes de serem observadas.

- **Qualidade por Leitura (Per sequence quality scores):**  
  Apresenta a distribuição das pontuações médias de qualidade para todas as leituras. Ajuda a identificar se uma grande proporção das leituras possui baixa qualidade geral.

- **Conteúdo de GC por Leitura (Per base GC content):**  
  Exibe a porcentagem de bases Guanina e Citosina em cada posição da leitura. Desvios significativos do esperado (com base no organismo) podem indicar contaminação ou viés de sequenciamento.

- **Conteúdo de GC por Sequência (Per sequence GC content):**  
  Mostra a distribuição do conteúdo de GC médio para todas as leituras.

- **Deduplicação de Sequências (Sequence Duplication Levels):**  
  Estima a proporção de leituras duplicadas, que podem ser um indicativo de super-amplificação por PCR.

- **Conteúdo de Adaptadores (Adapter Content):**  
  Detecta a presença de sequências de adaptadores de sequenciamento, que precisam ser removidas no pré-processamento.

- **Conteúdo N (N Content):**  
  Indica a porcentagem de bases "N" (bases não identificadas) em cada posição. Um alto percentual de Ns pode indicar problemas no sequenciamento.

Cada métrica é sinalizada com um status (**Normal**, **Aviso** ou **Falha**), facilitando a identificação de áreas problemáticas.

[Documentação detalhada do FastQC.](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### MultiQC

Embora o FastQC seja excelente para relatórios individuais, a análise de centenas ou milhares de amostras seria impraticável. É aqui que o **MultiQC** se torna indispensável.

**Como funciona:**

O MultiQC é uma ferramenta que agrega e sumariza os relatórios de controle de qualidade de múltiplas amostras, gerados por diversas ferramentas de bioinformática (incluindo FastQC, SAMtools, GATK, Trimmomatic, BWA, entre outras). Ele lê os arquivos de saída dessas ferramentas e os consolida em um único relatório HTML interativo.

![MultiQC](multiqc.gif "Exemplo de relatório de qualidade gerado no MultiQC a partir de vários arquivos FASTQ.")

**Vantagens do MultiQC:**

- **Visão Geral Consolidada:**  
  Em vez de abrir dezenas ou centenas de relatórios FastQC individualmente, o MultiQC apresenta todas as métricas em um único *dashboard*, facilitando a comparação entre amostras.

- **Detecção de Outliers:**  
  Permite identificar rapidamente amostras com problemas de qualidade que se destacam do restante do conjunto de dados.

- **Interatividade:**  
  O relatório HTML gerado pelo MultiQC é interativo, permitindo que o usuário explore gráficos, filtre amostras e visualize detalhes específicos de cada métrica.

- **Suporte a Múltiplas Ferramentas:**  
  Além do FastQC, o MultiQC pode integrar resultados de outras etapas do pipeline, como o mapeamento (BWA-MEM), o processamento de BAM (SAMtools), e a chamada de variantes (GATK), fornecendo uma visão holística do controle de qualidade de todo o fluxo de trabalho.

[Documentação detalhada do MultiQC.](https://docs.seqera.io/multiqc)

{{% /steps %}}

## Próximo passo

Agora que revisamos os fundamentos das métricas de qualidade de arquivos de sequenciamento de DNA, siga para a prática:

{{< cards >}}
  {{< card url="../guide/_index" title="Mão na Massa!" icon="document-duplicate" >}}
{{< /cards >}}
