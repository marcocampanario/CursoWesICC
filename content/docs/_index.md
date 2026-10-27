---
linkTitle: Práticas
title: Sobre o curso
---

👋 Bem-vindes ao nosso curso de genômica humana focado em análise de dados de sequenciamento de exoma completo!

<!--more-->

## Próximo passo

{{< cards >}}
  {{< card url="revisao-teorica" title="Revisão Teórica" icon="document-text" subtitle="O que você precisa saber antes de começar as práticas!" >}}
{{< /cards >}}

## Mão na massa!

{{< cards >}}
  {{< card url="https://marcocampanario.github.io/CursoWesICC/docs/guide/praticai" title="Prática I" icon="document-duplicate" >}}
  {{< card url="https://marcocampanario.github.io/CursoWesICC/docs/guide/praticaii" title="Prática II" icon="document-duplicate" >}}
  {{< card url="https://marcocampanario.github.io/CursoWesICC/docs/guide/praticaiii" title="Prática III" icon="document-duplicate" >}}
{{< /cards >}}

## Prática Linux

{{< cards >}}
  {{< card url="https://learn.microsoft.com/pt-br/windows/wsl/install" title="Instalando WSL" icon="command-line" subtitle="Guia oficial da Microsoft para instalar o Subsistema do Windows para Linux (WSL)." >}}
  {{< card url="https://marcocampanario.github.io/CursoWesICC/aula_linux.zip" title="Arquivos da aula" icon="arrow-down-tray" subtitle="Baixe os arquivos utilizados na prática de Linux." >}}
{{< /cards >}}

# Roteiro da aula prática: Linux para Bioinformática

## Etapa 0. Preparação

\`\`\`bash
unzip aula_linux.zip
cd aula_linux
pwd
ls
\`\`\`

## Etapa 1. Navegação

\`\`\`bash
pwd                      # onde estou?
ls                       # o que há aqui?
ls -l                    # formato detalhado (opção)
ls -lh dados              # tamanhos legíveis (-h) e argumento (pasta)
cd dados                  # entrar (caminho relativo)
pwd
cd ..                     # subir um nível
cd ~                      # ir para a pasta pessoal
cd aula_linux/dados       # voltar (tente usar Tab!)
\`\`\`

Pontos-chave: caminho absoluto vs. relativo, \`.\` e \`..\`, Tab para autocompletar, ↑ para repetir.

## Etapa 2. Organizar arquivos

\`\`\`bash
cd ~/aula_linux          # ajuste se descompactou em outro lugar
mkdir resultados         # cria a pasta onde salvaremos as saídas
cp dados/genes.tsv resultados/
ls resultados
mv resultados/genes.tsv resultados/genes_copia.tsv   # mv também renomeia
touch teste.txt          # cria arquivo vazio
rm teste.txt             # apaga SEM lixeira
\`\`\`

Ponto-chave: demonstrar o \`rm\` só em arquivo de teste. Mencionar \`rm -r\` (pastas) como perigoso.

## Etapa 3. Olhar o conteúdo

\`\`\`bash
cd dados
cat log_pipeline.txt         # arquivo pequeno: mostra tudo
head sequencias.fasta        # 10 primeiras linhas
head -n 4 amostras.fastq     # 1 read = 4 linhas
tail -n 3 log_pipeline.txt   # últimas linhas (onde o log termina)
less sequencias.fasta        # navegar com setas; q para sair
wc -l amostras.fastq         # 96 linhas
wc -l *                      # o * representa todos os arquivos
\`\`\`

Ponto-chave: por que \`cat\` não serve para arquivos de GBs, e \`head\`/\`less\` sim.

## Etapa 4. Buscar com grep

\`\`\`bash
grep ">" sequencias.fasta                    # só os cabeçalhos
grep -c ">" sequencias.fasta                 # -c conta: 10 sequências
grep "Escherichia" metadados.csv             # filtra linhas
grep -i "error" log_pipeline.txt              # -i ignora maiúsc./minúsc.
grep -n ERROR log_pipeline.txt                # -n mostra o nº da linha
grep -v ">" sequencias.fasta | head -3        # -v inverte: só as sequências
grep -A1 "qualidade=baixa" sequencias.fasta   # -A1: linha seguinte também
grep -c GAATTC sequencias.fasta               # busca de motivo (sítio EcoRI): 5 linhas
grep -c "^@READ" amostras.fastq               # ^ = início da linha: 24 reads
\`\`\`

## Etapa 5. Combinar e salvar: pipes e redirecionamento

\`\`\`bash
grep ">" sequencias.fasta | wc -l                            # 10
cut -d"," -f2 metadados.csv | tail -n +2 | sort | uniq -c   # contagem por organismo
cut -f2,6 genes.tsv | head -5                                # colunas separadas por TAB
grep ERROR log_pipeline.txt > ../resultados/problemas.txt    # > cria/sobrescreve
grep WARNING log_pipeline.txt >> ../resultados/problemas.txt # >> acrescenta
cat ../resultados/problemas.txt
\`\`\`

Ponto-chave: comandos simples + pipe = análise. \`cut -d\` define o separador; \`-f\` escolhe a coluna.

## Etapa 6. Desafio

1. Quantas sequências há em \`sequencias.fasta\`?
2. Quantos reads há em \`amostras.fastq\`? (duas formas: \`wc -l\` ÷ 4, ou \`grep -c "^@READ"\`)
3. Quantas sequências têm \`qualidade=baixa\`?
4. Quantas linhas \`ERROR\` há no log?
5. Quantos genes estão na fita negativa (\`-\`)? Dica: \`grep -c -P "\t-\t" genes.tsv\`
6. Quantas amostras de *Escherichia* foram aprovadas? Dica: dois \`grep\` em pipe.
7. Salve em \`resultados/\` a lista de genes do plasmídeo (\`grep plasmideo1 genes.tsv\`).
