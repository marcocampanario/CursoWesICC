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

## Explorando e Interpretando Arquivos VCF

{{% steps %}}

### Depois do FASTQ [...] Antes do VCF

Após o controle de qualidade dos arquivos FASTQ, os dados de sequenciamento são mapeados em um genoma humano de referência (GRCh38 / HG38 ou GRCh37/HG19), tornando-se arquivos SAM (*Sequence Alignment Map*) que contém as informações sobre a qualidade dos mapeamentos de leituras.

![SAM](sam.jpg "Estrutura de um arquivo SAM.")

Devido ao grande volume de dados ocupado por esse formato, transforma-se e trabalha-se com sua versão binária BAM (*Binary Alignment Map*), para todas as etapas de controle de qualidade pós mapeamento que antecedem a chamada e genotipagem de variantes genéticas.

Após a chamada de variantes genéticas e genotipagem conjunta (*joint calling*), as variações genotípicas em divergência com o genoma de referência são anotadas no padrão do formato VCF (*Variant Call Format*).

### Arquivo VCF

É o formato padrão para armazenar informações sobre variantes genéticas (SNPs, INDELs e outras).

![VCF](vcf.png "Estrutura de um arquivo VCF.")

#### Colunas Essenciais do formato VCF

- **CHROM**: Cromossomo onde a variante está localizada.
- **POS**: Posição inicial da variante no cromossomo.
- **ID**: Identificador da variante (ex: rsID do dbSNP, se disponível).
- **REF**: Alelo de referência na posição.
- **ALT**: Alelo(s) alternativo(s) observado(s).
- **QUAL**: Escore de qualidade Phred para a chamada da variante.
- **FILTER**: Indicação se a variante passou pelos filtros de qualidade.
- **INFO**: Campo de informações adicionais sobre a variante (densamente anotado).
- **FORMAT**: Define a ordem e o tipo de dados para cada amostra.
- **Colunas de Amostras**: Contêm os genótipos e outros dados específicos da amostra.

#### Interpretação dos Genótipos

- **0/0** - Homozigoto para o alelo de referência.
- **0/1** - Heterozigoto (um alelo de referência e um alelo alternativo).
- **1/1** - Homozigoto para o alelo alternativo.
- **./.** - Genótipo não chamado (leituras insuficientes ou de baixa qualidade).

#### Campo INFO

Contém as diversas anotações sobre cada variante, como:

- **QD (Quality by Depth)**: Escore de qualidade normalizado pela profundidade de cobertura.
- **FS (FisherStrand)**: Escore para viés de strand.
- **MQ (MappingQuality)**: Qualidade de mapeamento das leituras de suporte.

#### Anotadores

Anotadores, como exemplo o ANNOVAR, podem agregar informações de contexto genômico ao campo INFO do VCF:

- **Anotações funcionais**: gene, tipo de variante genética, mudança de aminoácido (caso *missense*), etc.
- **Frequências populacionais** (gnomAD, 1000 Genomes, ABraOM).
- **Classificações de patogenicidade** (ClinVar, REVEL, CADD).


{{% /steps %}}

## Visualizando dados com a linguagem R

{{% steps %}}

### O que é a Linguagem R?

R é uma linguagem de programação e um ambiente de software livre, amplamente utilizado para computação estatística e gráficos. É uma ferramenta poderosa para:

- **Processamento de Dados**: Importar, limpar, transformar e organizar grandes conjuntos de dados.
- **Análise Estatística**: Realizar uma vasta gama de análises estatísticas, desde testes simples a modelos complexos.
- **Visualização de Dados**: Criar gráficos de alta qualidade para explorar e apresentar resultados de forma eficaz.

### Por que R na Bioinformática?

No contexto da bioinformática, R é ideal para:

- Processar e analisar dados de sequenciamento (genômica, transcriptômica, etc.).
- Realizar análises de expressão diferencial.
- Interpretar e visualizar resultados de anotação de variantes genéticas.
- Desenvolver modelos preditivos para características biológicas ou doenças.
- Realizar análises estatísticas.

### O que é RStudio?

Enquanto R é a linguagem e o ambiente de execução, o RStudio é um Ambiente de Desenvolvimento Integrado (IDE) que facilita a escrita, execução e depuração de códigos R. Pense nele como uma "central de comando" que torna a experiência com R muito mais amigável e produtiva.

Principais funcionalidades do RStudio:

- **Editor de Código**: Onde você escreve e edita seus scripts R. Oferece recursos como realce de sintaxe e autocompletar.
- **Console**: Onde os comandos R são executados e os resultados são exibidos.
- **Ambiente**: Mostra todos os objetos (variáveis, dados, funções) que você criou na sessão atual.
- **Arquivos, Plots, Pacotes, Ajuda, Viewer**: Painéis dedicados para gerenciar arquivos, visualizar gráficos gerados, instalar e carregar pacotes, acessar a documentação e pré-visualizar arquivos HTML.

### Programação com Pacotes em R

Uma das maiores forças do R é sua vasta coleção de **pacotes**. Um pacote é um conjunto de funções, dados e documentação pré-escritos que estendem as capacidades básicas do R.

**Como funcionam os pacotes**:

1. **Instalação**: Antes de usar um pacote pela primeira vez, você precisa instalá-lo. Isso geralmente é feito com o comando `install.packages("nome_do_pacote")`. Pense nisso como "baixar" uma nova ferramenta para sua caixa de ferramentas R.
1. **Carregamento**: Após a instalação, para usar as funções de um pacote em uma sessão R, você precisa "carregá-lo" usando o comando `library(nome_do_pacote)`. Isso torna as funções do pacote disponíveis para uso.
1. **Uso**: Uma vez carregado, você pode chamar as funções do pacote diretamente em seu script.

**Exemplos de pacotes comuns e suas aplicações na bioinformática**:

- **`dplyr`**: Essencial para manipulação e transformação de dados (filtrar, selecionar, agrupar, etc.).
- **`ggplot2`**: O pacote padrão e mais robusto para a criação de gráficos estatísticos de alta qualidade e altamente personalizáveis.
- **`vcfR`**: Específico para bioinformática, permite ler, manipular e visualizar arquivos VCF (Variant Call Format).
- **`stringr`**: Facilita a manipulação de strings (textos), útil para extrair informações de anotações de variantes.

Ao utilizar pacotes, você não precisa "reinventar a roda" para cada tarefa. A comunidade R já desenvolveu e disponibilizou ferramentas otimizadas para diversas finalidades, permitindo que você se concentre na análise e interpretação dos seus dados.

{{% /steps %}}

## Próximo passo

Agora que revisamos os fundamentos das métricas de qualidade de arquivos de sequenciamento de DNA, siga para a prática:

{{< cards >}}
  {{< card url="../guide/" title="Mão na Massa!" icon="document-duplicate" >}}
{{< /cards >}}