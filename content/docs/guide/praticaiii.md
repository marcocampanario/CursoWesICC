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

### 01: Introdução ao R

Abra o RStudio, copie e cole o conteúdo abaixo e salve como `intro_R` (a extensão .R é automática).

```yaml
# Linhas iniciadas por # são comentários

# O R funciona como uma calculadora
2 + 3
10 / 2

# Mostrar o diretório de trabalho atual
getwd()

# Criar um objeto
numero_variantes <- 6

# Consultar seu conteúdo
numero_variantes

# Criar um vetor com regiões funcionais
regioes <- c("exonic", "intronic", "exonic",
             "splicing", "intronic", "exonic")

# Contar os elementos do vetor
length(regioes)

# Acessar o primeiro elemento
regioes[1]

# Criando um dataframe
vcf_exemplo <- data.frame(
  Func.refGene = regioes,
  CLNSIG = c("Pathogenic", "Benign", "Uncertain_significance",
             "Likely_pathogenic", "Benign", ".")
)

# Mostrar as primeiras linhas
head(vcf_exemplo)

# Mostrar os nomes das colunas
names(vcf_exemplo)

# Contar as linhas
nrow(vcf_exemplo)

# Acessar uma coluna pelo nome
vcf_exemplo$Func.refGene

# Quais linhas possuem variantes exônicas?
vcf_exemplo$Func.refGene == "exonic"

# Criar uma tabela contendo somente variantes exônicas
exonicas_exemplo <- subset(
  vcf_exemplo,
  Func.refGene == "exonic"
)

exonicas_exemplo

# Quais linhas possuem uma destas classificações?
vcf_exemplo$CLNSIG %in% c("Pathogenic", "Likely_pathogenic")

# Contar variantes por região
contagem_regiao_exemplo <- table(vcf_exemplo$Func.refGene)

contagem_regiao_exemplo

# Transformar a contagem em um data.frame
resumo_regiao_exemplo <- as.data.frame(contagem_regiao_exemplo)

resumo_regiao_exemplo

# Contar variantes por classificação
table(vcf_exemplo$CLNSIG)


#####

install.packages("writexl")
install.packages("dplyr")
install.packages("stringr")
install.packages("tidyr")
install.packages("glue")
install.packages("ggplot2")
install.packages("scales")
install.packages("data.table")


library(writexl)
library(dplyr)
library(stringr)
library(tidyr)
library(glue)
library(ggplot2)
library(scales)

```

### 01: Criar script R

Abra o RStudio, copie e cole o conteúdo abaixo e salve como `praticaiii` (a extensão .R é automática).

```yaml
#install.packages("dplyr")
#install.packages("writexl")

library(writexl)
library(dplyr)
library(stringr)
library(tidyr)
library(glue)
library(ggplot2)
library(scales)
library(data.table)

### Ctrl+Shift+H
### Selecione a pasta com os arquivos da prática como seu diretório de trabalho
getwd()

### PARTE 1: ORGANIZAÇÃO DOS DADOS

# Arquivo ANNOVAR .txt ---------------------------------------------------------

vcf <- data.table::fread("mgp-hg38.txt", sep = "\t", quote = "")

# Contar os tipos de variantes quanto à região ---------------------------------

tipo_var_regiao <- as.data.frame(table(vcf$Func.refGene))

# Selecionar variantes quanto à região -----------------------------------------

filt1 <- filter(vcf, Func.refGene == "exonic")
filt2 <- filter(vcf, Func.refGene == "intronic")
filt3 <- filter(vcf, Func.refGene == "splicing")

# Contar os tipos de variantes quanto à patogenicidade -------------------------

tipo_var_patogenicidade <- as.data.frame(table(vcf$CLNSIG))

# Selecionar variantes quanto à patogenicidade ---------------------------------

filt4 <- filter(vcf, CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Pathogenic,_low_penetrance|other", "Likely_pathogenic"))
filt5 <- filter(vcf, CLNSIG %in% c("Benign", "Benign/Likely_benign", "Likely_benign"))
filt6 <- filter(vcf, CLNSIG %in% c("Uncertain_significance"))

### PLOTS ----------------------------------------------------------------------

### PARTE 2: VISUALIZAÇÃO DOS DADOS

### PLOT - REGIAO FUNCIONAL
# Garantir que os dados estão organizados
tipo_var_regiao <- tipo_var_regiao %>%
  rename(Regiao = Var1, Contagem = Freq) %>%
  arrange(desc(Contagem))

# Criar gráfico de barras
ggplot(tipo_var_regiao, aes(x = reorder(Regiao, -Contagem), y = Contagem)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = format(Contagem, big.mark = ".", decimal.mark = ",")), 
            vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribuição de variantes por região funcional",
    x = "Região",
    y = "Número de variantes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ","))
### DÊ ZOOM E SALVE O GRÁFICO COMO "fig1"

### PLOT - PATOGENICIDADE
# Garantir que os dados estão organizados
tipo_var_patogenicidade <- tipo_var_patogenicidade %>%
  rename(Patogenicidade = Var1, Contagem = Freq) %>%
  arrange(desc(Contagem))

# Criar gráfico de barras
ggplot(tipo_var_patogenicidade, aes(x = reorder(Patogenicidade, -Contagem), y = Contagem)) +
  geom_bar(stat = "identity", fill = "tomato") +
  geom_text(aes(label = format(Contagem, big.mark = ".", decimal.mark = ",")), 
            vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribuição de variantes por patogenicidade (CLNSIG)",
    x = "Classificação",
    y = "Número de variantes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ","))
### DÊ ZOOM E SALVE O GRÁFICO COMO "fig2"

### PLOT - PATOGENICIDADE (excluindo não classificadas)
# Garantir que os dados estão organizados
tipo_var_patogenicidade_filtrado <- tipo_var_patogenicidade %>%
  filter(!is.na(Patogenicidade) & Patogenicidade != ".") %>%
  arrange(desc(Contagem))

# Criar gráfico de barras
ggplot(tipo_var_patogenicidade_filtrado, aes(x = reorder(Patogenicidade, -Contagem), y = Contagem)) +
  geom_bar(stat = "identity", fill = "tomato") +
  geom_text(aes(label = format(Contagem, big.mark = ".", decimal.mark = ",")), 
            vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribuição de variantes por patogenicidade (CLNSIG - somente classificadas)",
    x = "Classificação",
    y = "Número de variantes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::label_number(big.mark = ".", decimal.mark = ","))
### DÊ ZOOM E SALVE O GRÁFICO COMO "fig3"

### EXCEL ----------------------------------------------------------------------

### PARTE 3: EXPLORAÇÃO DOS BANCOS DE DADOS PÚBLICOS

# Salvar o objeto filt4 como tabela .xlsx
writexl::write_xlsx(filt4, "variantes_patogenicas.xlsx")

# Criar objeto com variantes comuns
variantes_comuns <- vcf %>%
  mutate(AF = as.numeric(na_if(as.character(AF), "."))) %>%
  filter(!is.na(AF), AF > 0.05) %>%
  filter(Func.refGene == "exonic")

# Conferir quantas variantes atendem ao critério
nrow(variantes_comuns)

# Fixar a semente para reproduzir o mesmo sorteio
set.seed(555)

# Sortear 100 linhas sem reposição, mantendo todas as colunas
variantes_comuns_100 <- variantes_comuns %>%
  slice_sample(n = 100, replace = FALSE)

# Exportar para Excel no diretório de trabalho atual
writexl::write_xlsx(
  variantes_comuns_100,
  "variantes_comuns.xlsx"
)

```

### 02: Executar o software RStudio na interface gráfica X2GO

Agora podemos trabalhar com o script em R.

#### Executar linha por linha
- Para executar **uma linha ou uma seleção** do script:
  - Clique na linha desejada (ou selecione o trecho).
  - Pressione **Ctrl + Enter** (ou **Cmd + Enter** no Mac).
  - O comando será executado no **Console**, que fica no painel inferior esquerdo.

## Atividade 03

- Considerando a classificação de posição funcional das variantes genéticas no genoma, qual o tipo de variante genética mais abundante entre as amostras sequenciadas?

- E quanto à classificação de patogenicidade, qual o tipo de variante genética mais abundante?

- Escolha uma variante genética classificada como patogênica e explore suas características no [Ensembl](https://www.ensembl.org/index.html) e no [UCSC Genome Browser](https://genome.ucsc.edu/). Lembre-se, todas as variantes genéticas conhecidas possuem rsID no [dbSNP](https://www.ncbi.nlm.nih.gov/snp/).