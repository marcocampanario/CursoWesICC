---
title: Prática III
weight: 3
---

## Explorando e Interpretando Arquivos VCF

Nesta seção prática, vamos criar um script R para carregar o arquivo VCF, extrair informações de anotação (especialmente o impacto funcional e a classificação de patogenicidade) e gerar gráficos.


**Objetivos**:

- Conhecer a interface gráfica básica do RStudio e princípios da linguagem R.
- Aplicar pacotes de manipulação de dados pelo R para quantificar e analisar variantes genéticas.
- Gerar gráficos para apresentação de resultados.


## Passo-a-passo

### 01: Criar script R

Abra o bloco de notas e salve o conteúdo abaixo como `analyze_vcf_annotations.R`:

```yaml
# Instalar pacotes se necessário
# install.packages("vcfR")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("stringr")

# Carregar pacotes
library(vcfR)
library(ggplot2)
library(dplyr)
library(stringr)

# --- Caminho para o arquivo VCF ---

vcf_file <- "~/data/variants/example.vcf"

# Verificar se o arquivo existe
if (!file.exists(vcf_file)) {
  stop("Arquivo VCF não encontrado. Verifique o caminho: ", vcf_file)
}

# Carregar o arquivo VCF
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Extrair informações do campo INFO (onde ANNOVAR e outras anotações estão)
# Nota: O formato exato das anotações no campo INFO pode variar.
# Para ANNOVAR, as anotações são frequentemente concatenadas e separadas por ponto e vírgula.
# Este é um exemplo genérico; você pode precisar ajustar baseado no seu VCF.

# Converta o VCF para um data.frame para facilitar a manipulação
# Isso é uma forma simplificada, em um cenário real, você pode precisar parsing mais robusto do INFO

vcf_df <- vcfR2tidy(vcf) %>%

  unnest_wider(c(Fix, Info)) %>% # Expande as colunas Fix e Info

  select(CHROM, POS, ID, REF, ALT, FILTER, starts_with("ANNOVAR"), # Inclua colunas de ANNOVAR

         starts_with("ClinVar"), # Inclua colunas de ClinVar

         starts_with("gnomAD"), # Inclua colunas de gnomAD

         starts_with("REVEL")) # Inclua colunas de REVEL


# --- Exemplos de Análise e Plotagem ---

# 1. Distribuição de Variantes por Impacto Funcional (ANNOVAR)

# Supondo que "ANNOVAR_Func.refGene" e "ANNOVAR_ExonicFunc.refGene" existam.

# Você pode precisar ajustar os nomes das colunas de acordo com o seu VCF.

if ("ANNOVAR_ExonicFunc.refGene" %in% colnames(vcf_df)) {

  # Limpar e agrupar as funções exonicas para plotagem

  functional_impact <- vcf_df %>%

    filter(FILTER == "PASS") %>% # Apenas variantes que passaram nos filtros

    mutate(ExonicFunc = case_when(

      grepl("frameshift", ANNOVAR_ExonicFunc.refGene, ignore.case = TRUE) ~ "Frameshift",

      grepl("stopgain", ANNO
```

### 02: Executar o software RStudio na interface gráfica X2GO

…

Descrever a execução por linha no RStúdio (Ctrl+Enter), pwd e setpwd
 Salvar plot e afins.
