---
title: Prática Linux
weight: 4
---

## Linux para Bioinformática

Nesta aula prática, vamos explorar comandos básicos do Linux aplicados a arquivos comuns em Bioinformática.

## Passo-a-passo

### 00: Preparação

{{% steps %}}

```bash
unzip aula_linux.zip
cd aula_linux
pwd
ls
```

{{% /steps %}}

### 01: Navegação

{{% steps %}}

```bash
pwd                      # onde estou?
ls                       # o que há aqui?
ls -l                    # formato detalhado (opção)
ls -lh dados             # tamanhos legíveis (-h) e argumento (pasta)
cd dados                 # entrar (caminho relativo)
pwd
cd ..                    # subir um nível
cd ~                     # ir para a pasta pessoal
cd aula_linux/dados      # voltar (tente usar Tab!)
```

{{% /steps %}}

**Pontos-chave:** caminho absoluto vs. relativo, `.` e `..`, Tab para autocompletar, ↑ para repetir.

### 02: Organizar arquivos

{{% steps %}}

```bash
cd ~/aula_linux          # ajuste se descompactou em outro lugar
mkdir resultados         # cria a pasta onde salvaremos as saídas
cp dados/genes.tsv resultados/
ls resultados
mv resultados/genes.tsv resultados/genes_copia.tsv   # mv também renomeia
touch teste.txt          # cria arquivo vazio
rm teste.txt             # apaga SEM lixeira
```

{{% /steps %}}

**Ponto-chave:** demonstrar o `rm` só em arquivo de teste. Mencionar `rm -r` (pastas) como perigoso.

### 03: Olhar o conteúdo

{{% steps %}}

```bash
cd dados
cat log_pipeline.txt         # arquivo pequeno: mostra tudo
head sequencias.fasta        # 10 primeiras linhas
head -n 4 amostras.fastq     # 1 read = 4 linhas
tail -n 3 log_pipeline.txt   # últimas linhas (onde o log termina)
less sequencias.fasta        # navegar com setas; q para sair
wc -l amostras.fastq         # 96 linhas
wc -l *                      # o * representa todos os arquivos
```

{{% /steps %}}

**Ponto-chave:** por que `cat` não serve para arquivos de GBs, e `head`/`less` sim.

### 04: Buscar com `grep`

{{% steps %}}

```bash
grep ">" sequencias.fasta                   # só os cabeçalhos
grep -c ">" sequencias.fasta                # -c conta: 10 sequências
grep "Escherichia" metadados.csv            # filtra linhas
grep -i "error" log_pipeline.txt            # -i ignora maiúsc./minúsc.
grep -n ERROR log_pipeline.txt              # -n mostra o nº da linha
grep -v ">" sequencias.fasta | head -3      # -v inverte: só as sequências
grep -A1 "qualidade=baixa" sequencias.fasta # -A1: linha seguinte também
grep -c GAATTC sequencias.fasta             # busca de motivo (sítio EcoRI): 5 linhas
grep -c "^@READ" amostras.fastq             # ^ = início da linha: 24 reads
```

{{% /steps %}}

### 05: Combinar e salvar: pipes e redirecionamento

{{% steps %}}

```bash
grep ">" sequencias.fasta | wc -l                           # 10
cut -d"," -f2 metadados.csv | tail -n +2 | sort | uniq -c   # contagem por organismo
cut -f2,6 genes.tsv | head -5                               # colunas separadas por TAB
grep ERROR log_pipeline.txt > ../resultados/problemas.txt   # > cria/sobrescreve
grep WARNING log_pipeline.txt >> ../resultados/problemas.txt # >> acrescenta
cat ../resultados/problemas.txt
```

{{% /steps %}}

**Ponto-chave:** comandos simples + pipe = análise. `cut -d` define o separador; `-f` escolhe a coluna.

## Desafio

1. Quantas sequências há em `sequencias.fasta`?
2. Quantos reads há em `amostras.fastq`? (duas formas: `wc -l` ÷ 4, ou `grep -c "^@READ"`)
3. Quantas sequências têm `qualidade=baixa`?
4. Quantas linhas `ERROR` há no log?
5. Quantos genes estão na fita negativa (`-`)? Dica: `grep -c -P "\t-\t" genes.tsv`
6. Quantas amostras de *Escherichia* foram aprovadas? Dica: dois `grep` em pipe.
7. Salve em `resultados/` a lista de genes do plasmídeo (`grep plasmideo1 genes.tsv`).
