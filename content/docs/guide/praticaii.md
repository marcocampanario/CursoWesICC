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

