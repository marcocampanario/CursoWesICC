---
title: Prática II
weight: 2
---

## Explorando e Interpretando Arquivos VCF

Nesta seção prática, vamos utilizar estratégias de bioniformática para analisar arquivos VCF.

Para isso, utilizaremos um arquivo VCF pequeno e anotado, contendo algumas variantes patogênicas e de significado incerto para fins de demonstração. Este arquivo está localizado em `~/data/variants/example.vcf.gz` no servidor.


**Objetivos**:

- Compreender a estrutura de um arquivo VCF (Variant Call Format).
- Interpretar as informações contidas nas diferentes colunas de um VCF, incluindo genótipos e anotações.
- Utilizar comandos de linha e scripts R para inspecionar e filtrar variantes.

## Passo-a-passo

### 01: Descomprimir um Arquivo VCF e Visualizar seu Cabeçalho

Antes de tudo, vamos descompactar o arquivo VCF de exemplo e inspecionar seu cabeçalho. O cabeçalho VCF contém metadados importantes sobre o pipeline de chamada de variantes e as ferramentas de anotação usadas.

{{% steps %}}

```yaml
# Navegar até o diretório das variantes

cd ~/data/variants/
```

```yaml
# Descompactar o arquivo VCF de exemplo (se for .gz)

gunzip -k example.vcf.gz 
```

```yaml
# Visualizar as primeiras linhas do arquivo VCF, incluindo o cabeçalho

# As linhas com '##' são metadados, a linha com '#' define as colunas.

head -n 50 example.vcf 
```

{{% /steps %}}

### 02: Inspecionar o Conteúdo do VCF com `grep` e `less`

Vamos usar `less` para navegar pelo arquivo VCF e `grep` para buscar por variantes específicas ou filtrar por critérios simples.

{{% steps %}}

```yaml
# Visualizar o arquivo VCF completo usando less (pressione 'q' para sair)

less example.vcf
```

```yaml
# Procurar por variantes em um cromossomo específico (ex: chr1)

grep -E '^#|chr1\t' example.vcf | less
```

```yaml
# Procurar por variantes que passaram em todos os filtros (campo FILTER = PASS)

grep -E '^#|PASS\t' example.vcf | less
```

```yaml
# Procurar por uma variante com um ID específico (se houver, ex: rs12345)

grep -E '^#|rs12345\t' example.vcf | less
```

{{% /steps %}}

### 03: Extrair e Filtrar Variantes com `bcftools`

O pacote de ferramentas `bcftools` é um poderoso aliado na exploração de arquivos VCF. Vamos usá-lo para extrair informações específicas e aplicar filtros. Lembre-se que ele já está instalado no servidor de aula, mas deve ser instalado em seu ambiente Linux para uso.

[Sobre instalação do BCFtools.](https://samtools.github.io/bcftools/howtos/install.html)

{{% steps %}}

```yaml
# Contar o número total de variantes (excluindo linhas de cabeçalho)

bcftools view -H example.vcf | wc -l
```

```yaml
# Filtrar variantes que são classificadas como "Pathogenic" ou "Likely_pathogenic" nas colunas de anotação 
# (exemplo de anotação ClinVar)
# NOTA: O nome exato da anotação no campo INFO ou em coluna específica pode variar dependendo da ferramenta de anotação 
# (ex: ANN, ClinVar_CLNSIG)
# Para este exemplo, usaremos uma anotação 'ClinVar_Significance'
# -i (ou --include): Indica que você quer incluir apenas as variantes que satisfazem a condição especificada
# 

bcftools view -i 'INFO/ClinVar_Significance="Pathogenic" | INFO/ClinVar_Significance="Likely_pathogenic"' example.vcf | less
```

```yaml
# Filtrar por SNPs e INDELs

# SNPs

bcftools view -v snps example.vcf | less

# INDELs

bcftools view -v indels example.vcf | less 

```

{{% /steps %}}

