---
title: Prática I
weight: 1
---

## Análise de arquivos FASTQ

Nesta seção prática, vamos utilizar as ferramentas FastQC e MultiQC para analisar a qualidade de arquivos FASTQ de exemplo.

Os arquivos FASTQ estão localizados no diretório `~/data/fastq/` no servidor.

**Objetivos**:

- Gerar relatórios de qualidade para arquivos FASTQ individuais usando FastQC.
- Consolidar múltiplos relatórios FastQC em um único relatório interativo usando MultiQC.
- Interpretar os resultados para identificar potenciais problemas de qualidade.


## Passo-a-passo

### 01: Navegar até o diretório dos arquivos FASTQ

Primeiro, acesse o diretório onde seus arquivos FASTQ estão armazenados.

```yaml
cd ~/data/fastq/
```

### 02: Executar FastQC para cada arquivo FASTQ

Para cada arquivo FASTQ, execute o FastQC individualmente. O FastQC criará um arquivo `.html` e um arquivo `.zip` para cada FASTQ no mesmo diretório.

```yaml
cd ~/data/fastq/
```

Page `example.md`:

```yaml
fastqc sample1.fastq.gz
```

```yaml
fastqc sample2.fastq.gz
```

```yaml
# Repita para todos os arquivos FASTQ que você deseja analisar!
```

### 03: Executar MultiQC para consolidar os relatórios

Após executar o FastQC para todos os seus arquivos, execute o MultiQC no diretório que contém os resultados do FastQC (os arquivos `.html` ou `.zip` gerados). O MultiQC irá procurar automaticamente pelos arquivos de saída do FastQC e gerar um relatório consolidado.

```yaml
multiqc .
```

- **Observação**: O ponto `.` indica que o MultiQC deve procurar pelos arquivos de resultados no diretório atual. Se seus arquivos FastQC estiverem em um subdiretório específico, você pode especificar o caminho (ex: `multiqc ./fastqc_results/`).

### 04: Visualizar o relatório MultiQC

O MultiQC gerará um arquivo chamado `multiqc_report.html` no diretório onde foi executado. Você precisará transferir este arquivo para o seu computador local (usando `fpt` ou FileZilla, por exemplo) e abri-lo em um navegador web para interagir com o relatório.

### 05: Interpretação dos Resultados:

Ao abrir o `multiqc_report.html` no seu navegador, você verá um dashboard interativo. Preste atenção nas seguintes seções:

- *General Statistics*: Fornece um resumo de alto nível para cada amostra.

- *FastQC: Per Base Sequence Quality*: Verifique se a qualidade das bases se mantém alta ao longo da leitura. Quedas abruptas, especialmente no final, podem indicar a necessidade de corte de reads.

- *FastQC: Per Sequence Quality Scores*: Observe a distribuição das pontuações de qualidade médias. Idealmente, a maioria das leituras deve ter alta qualidade.

- *FastQC: Adapter Content*: Garanta que a quantidade de adaptadores seja mínima ou ausente. A presença de adaptadores indica a necessidade de remoção (usando ferramentas como Trimmomatic).

- *FastQC: Sequence Duplication Levels*: Níveis muito altos de duplicação podem sugerir super-amplificação por PCR e podem impactar a profundidade de cobertura real.

- *FastQC: Per Base GC Content*: Compare o conteúdo de GC com o esperado para o organismo. Desvios podem indicar contaminação.

Os relatórios do MultiQC são interativos. Você pode clicar nas seções para expandir os gráficos, passar o mouse sobre os pontos de dados para obter mais informações e usar as opções de filtragem para visualizar subconjuntos de amostras.
